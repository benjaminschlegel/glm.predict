dc.lm = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
                  type = c("any", "simulation", "bootstrap"), summary = TRUE){
  # check inputs
  if(is.null(values) && (is.null(values1) || is.null(values2))){
    stop("Either values1 and values2 or values has to be specified!")
  }
  if(!is.null(values)){
    l = length(values)
    values1 = values[1 : (l/2)]
    values2 = values[(l/2 + 1) : l]
  }
  if(sum("lm" %in% class(model)) == 0){
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
    if(nrow(sigma) != length(values1)){
      warning("sigma and values do not match, ignoring the specified sigma")
      sigma = vcov(model)
    }
    if(!is.null(set.seed)){
      set.seed(set.seed)
    }
    betas_sim = MASS::mvrnorm(sim.count, coef(model), sigma)
    pred1 = betas_sim %*% values1
    pred2 = betas_sim %*% values2
  }else{ # bootstrap
    boot = function(x, model){
      data = model.frame(model)
      sample_data = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      if("(weights)" %in% colnames(data)){
        coef(update(model, data = sample_data, weights = `(weights)`))
      }else{
        coef(update(model, data = sample_data))
      }
    }
    betas_boot = do.call('rbind', lapply(seq_len(sim.count), boot, model))
    pred1 = betas_boot %*% values1
    pred2 = betas_boot %*% values2
  }
  
  diff = pred1 - pred2
  
  all = cbind(pred1, pred2, diff)
  
  # return all simulated / bootstrapped values if summary is FALSE
  if(!summary){
    return(all)
  }
  
  
  confint_lower = (1 - conf.int) / 2 
  result = apply(all, 2, quantile, probs = c(confint_lower, 1 - confint_lower))
  result = t(rbind(apply(all, 2, mean), result))
  
  colnames(result) = c("Mean", 
                       paste0(100 * confint_lower,"%"), 
                       paste0(100 * (1 - confint_lower),"%"))
  rownames(result) = c("Case 1", "Case 2", "Difference")
  result
}

