basepredict.mlogit = function(model,values,sim.count=1000,conf.int=0.95,sigma=NULL,set.seed=NULL,
                                type = c("any", "simulation", "bootstrap")){
  
  # check inputs
  if(sum("mlogit" %in% class(model)) == 0){
    stop("model has to be of type mlogit()")
  }
  choices = names(model$freq)
  beta_names = names(coef(model))
  n_multinomial = do.call(sum, lapply(choices, grepl, beta_names))
  n_conditional = length(beta_names) - n_multinomial
  if(length(values) != n_multinomial / (length(choices) - 1) + n_conditional){
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
    warning("Bootstrap not supported for mlogit. Using simulation instead.")
  }
  type = "simulation"
  
  
  betas = coef(model)
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  
  # simulation
  betas_draw = MASS::mvrnorm(sim.count, betas, sigma)


  pred = matrix(nrow = length(choices), ncol = sim.count)
  

  
  for(i in 1:sim.count){
    current_betas = betas_draw[i,]
    n = length(choices) - 1
    pos_conditional_vars = (n+1):(n+n_conditional)
    pos_multinomial_vars = max(pos_conditional_vars) + 
      (seq_len((length(betas) - max(pos_conditional_vars)) / n) - 1) * n 
    sim_temp = current_betas[seq_len(n)]
    for(pos in pos_conditional_vars){
      sim_temp = cbind(sim_temp, current_betas[pos])
    }
    
    for(j in pos_multinomial_vars){
      sim_temp = cbind(sim_temp, current_betas[j + seq_len(n)])
    }
    
    yhat = c(0, sim_temp %*% values)
    e = exp(yhat)
    
    for(j in seq_along(choices)){
      pred[j, i] = e[j] / sum(e)
    }
  }
  
  confint_lower = (1 - conf.int) / 2
  result = matrix(nrow = n + 1, ncol = 3)
  colnames(result) = c("mean", 
                       paste0(100 * confint_lower, "%"), 
                       paste0(100 * (1 - confint_lower), "%"))
  rownames(result) = choices
  for(j in 1:(n+1)){
    result[j,] = c(mean(pred[j,]), quantile(pred[j,],probs = c(confint_lower, 1 - confint_lower)))
  }
  result
}
