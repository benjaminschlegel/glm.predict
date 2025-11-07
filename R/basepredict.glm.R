basepredict.glm = function(model, values, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, 
                           type = c("any", "simulation", "bootstrap"), summary = TRUE){
  # check inputs
  if(sum("glm" %in% class(model)) == 0){
    stop("model has to be of type glm()")
  }
  if(length(values) != length(coef(model))){
    stop("the length of values is not identical to the number of coefficient of the model")
  }
  if(!is.numeric(sim.count) | round(sim.count) != sim.count){
    stop("sim.count has to be whole number")
  }
  if(!is.numeric(conf.int)){
    stop("conf.int has to be numeric")
  }
  if(!is.null(set.seed) & !is.numeric(set.seed)){
    stop("set.seed must be numeric")
  }
  
  type = match.arg(type)
  
  if(type == "bootstrap" && "svyglm" %in% class(model)){
    warning("Boostrap not supported for survey()-models, using simulations instead.")
    type = "simulation"
  }

  # model type
  model.type = family(model)
  link = model.type[2]
  
  if(type == "any"){
    if("svyglm" %in% class(model)){
      type = "simulation"
    }else if(nrow(model$data) < 500){
      type = "bootstrap"
      message("Type not specified: Using bootstrap as n < 500")
    }else{
      type = "simulation"
      message("Type not specified: Using simulation as n >= 500")
    }
  }
  
  if(type == "simulation"){
    if(is.null(sigma)){
      sigma = stats::vcov(model)
    }
    if(nrow(sigma) != length(values)){
      warning("sigma and values do not match, ignoring the specified sigma")
      sigma = stats::vcov(model)
    }
    if(!is.null(set.seed)){
      set.seed(set.seed)
    }
    betas_sim = MASS::mvrnorm(sim.count, coef(model), sigma)
    # get the predicted probabilities/values with the inverse link function
    pred = calculate_glm_pred(betas_sim, values, link)
  }else{ # bootstrap
    boot = function(x, model){
      data = model$data
      sample_data = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      if("(weights)" %in% colnames(data)){
        w <- sample_data[["(weights)"]]
        coef(update(model, data = sample_data, weights = w))
      }else{
        coef(update(model, data = sample_data))
      }
    }
    betas_boot = do.call('rbind', lapply(seq_len(sim.count), boot, model))
    # get the predicted probabilities/values with the inverse link function
    pred = calculate_glm_pred(betas_boot, values, link)
  }
  
  # return all simulated / bootstrapped values if summary is FALSE
  if(!summary){
    return(pred)
  }
  
  # calculate mean and confident interval
  confint_lower = (1 - conf.int) / 2 
  result = t(as.matrix(c(mean(pred, na.rm = TRUE),
                          quantile(pred, c(confint_lower, 1 - confint_lower), na.rm = TRUE))))
  # name the output matrix
  colnames(result) = c("Mean", 
                        paste0(100 * confint_lower,"%"), 
                        paste0(100 * (1 - confint_lower),"%"))
  result
}

calculate_glm_pred = function(betas, x, link){
  yhat = betas %*% x
  
  # the inverse link functions
  if(link == "logit"){
    return(exp(yhat) / (1 + exp(yhat)))
  }
  if(link == "log"){
    return(exp(yhat))
  }
  if(link == "identity"){
    return(yhat)
  }
  if(link == "probit"){ 
    return(pnorm(yhat))
  }
  if(link == "cauchit"){
    return(tan(pi * (yhat - 0.5)))
  }
  if(link == "cloglog"){
    return(exp(-exp(yhat)) * (-1 + exp(exp(yhat))))
  }
  if(link == "sqrt"){
    return(yhat * yhat)
  }
  if(link == "1/mu^2"){
    return(1 / sqrt(yhat))
  }
  if(link == "inverse"){
    return(1 / yhat)
  }  
}
