basepredict.tobit = function(model, values, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL,
                             type = c("any", "simulation", "bootstrap"), summary = TRUE){
  # check inputs
  if(sum("tobit" %in% class(model)) == 0){
    stop("model has to be of type lm()")
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
  
  if(type == "bootstrap"){
    warning("Boostrap not supported for tobit()-models, using simulations instead.")
    type = "simulation"
  }
  
  if(type %in% c("simulation", "any")){
    if(is.null(sigma)){
      sigma = stats::vcov(model)
      sigma = sigma[-nrow(sigma), ]
      sigma = sigma[, -ncol(sigma)]
    }
    if(nrow(sigma) != length(values)){
      warning("sigma and values do not match, ignoring the specified sigma")
      sigma = stats::vcov(model)
      sigma = sigma[-nrow(sigma), ]
      sigma = sigma[, -ncol(sigma)]
    }
    if(!is.null(set.seed)){
      set.seed(set.seed)
    }
    betas_sim = MASS::mvrnorm(sim.count, coef(model), sigma)
    # get the predicted probabilities/values
    yhat = betas_sim %*% values
  }
    
    pred = yhat * pnorm(yhat / model$scale) 
    
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
