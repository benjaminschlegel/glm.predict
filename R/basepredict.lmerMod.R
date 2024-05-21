basepredict.lmerMod = function(model, values, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, 
                           type = c("any", "simulation", "bootstrap"), summary = TRUE){
  # check inputs
  if(sum("lmerMod" %in% class(model)) == 0){
    stop("model has to be of type lmer()")
  }
  if(length(values) != length(fixef(model))){
    stop("the length of values is not identical to the number of coefficient of the model")
  }
  if(!is.numeric(sim.count) | round(sim.count) != sim.count){
    stop("sim.count has to be a whole number")
  }
  if(!is.numeric(conf.int)){
    stop("conf.int has to be numeric")
  }
  if(!is.null(set.seed) & !is.numeric(set.seed)){
    stop("set.seed must be numeric")
  }
  
  type = match.arg(type)
  
  if(type == "any"){
    if(nrow(model.frame(model)) < 500){
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
    betas_sim = MASS::mvrnorm(sim.count, fixef(model), sigma)
    # get the predicted probabilities/values with the inverse link function
    pred = betas_sim %*% values
  }else{ # bootstrap
    boot = function(x, model){
      data = model.frame(model)
      sample_data = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      if("(weights)" %in% colnames(data)){
        fixef(update(model, data = sample_data, weights = `(weights)`))
      }else{
        fixef(update(model, data = sample_data))
      }
    }
    betas_boot = do.call('rbind', lapply(seq_len(sim.count), boot, model))
    # get the predicted probabilities/values with the inverse link function
    pred = betas_boot %*% values
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
