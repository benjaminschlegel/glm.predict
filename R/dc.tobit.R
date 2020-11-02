dc.tobit = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
                 type = c("any", "simulation", "bootstrap")){
  # check inputs
  if(is.null(values) && (is.null(values1) || is.null(values2))){
    stop("Either values1 and values2 or values has to be specified!")
  }
  if(!is.null(values)){
    l = length(values)
    values1 = values[1 : (l/2)]
    values2 = values[(l/2 + 1) : l]
  }
  if(sum("tobit" %in% class(model)) == 0){
    stop("model has to be of type lm()")
  }
  if(length(values1) != length(coef(model))){
    stop("the length of values1 is not identical to the number of coefficient of the model")
  }
  if(length(values2) != length(coef(model))){
    stop("the length of values2 is not identical to the number of coefficient of the model")
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
    if(nrow(sigma) != length(values1)){
      warning("sigma and values do not match, ignoring the specified sigma")
      sigma = vcov(model)
      sigma = sigma[-nrow(sigma), ]
      sigma = sigma[, -ncol(sigma)]
    }
    if(!is.null(set.seed)){
      set.seed(set.seed)
    }
    betas_sim = MASS::mvrnorm(sim.count, coef(model), sigma)
    yhat1 = betas_sim %*% values1
    yhat2 = betas_sim %*% values2
  }
  
  pred1 = yhat1 * pnorm(yhat1 / model$scale) 
  pred2 = yhat2 * pnorm(yhat2 / model$scale) 
  
  diff = pred1 - pred2
  
  all = cbind(pred1, pred2, diff)
  
  
  confint_lower = (1 - conf.int) / 2 
  result = apply(all, 2, quantile, probs = c(confint_lower, 1 - confint_lower))
  result = t(rbind(apply(all, 2, mean), result))
  
  colnames(result) = c("Mean", 
                       paste0(100 * confint_lower,"%"), 
                       paste0(100 * (1 - confint_lower),"%"))
  rownames(result) = c("Case 1", "Case 2", "Difference")
  result
}