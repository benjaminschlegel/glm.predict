dc.mlogit = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
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
  if(sum("mlogit" %in% class(model)) == 0){
    stop("model has to be of type mlogit()")
  }
  
  choices = names(model$freq)
  beta_names = names(coef(model))
  n_multinomial = do.call(sum, lapply(choices, grepl, beta_names))
  n_conditional = length(beta_names) - n_multinomial
  if(length(values1) != n_multinomial / (length(choices) - 1) + n_conditional){
    stop("the length of values1 is not identical to the number of coefficient of the model")
  }
  if(length(values2) != n_multinomial / (length(choices) - 1) + n_conditional){
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
  
  for(v in values){
    if(length(v) != length(choices)-1 & length(v) != 1){
      stop("values need to be of length 1 or length of choices-1")
    }
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
  
  x = matrix(ncol = 2, nrow = length(choices))
  pred1 = pred2 = matrix(nrow = length(choices), ncol = sim.count)
  
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
    
    yhat1 = c(0, diag(sim_temp %*% t(as.matrix(expand.grid(values1)))))
    yhat2 = c(0, diag(sim_temp %*% t(as.matrix(expand.grid(values2)))))
    e1 = exp(yhat1)
    e2 = exp(yhat2)
    
    for(j in seq_along(choices)){
      pred1[j, i] = e1[j] / sum(e1)
      pred2[j, i] = e2[j] / sum(e2)
    }
  }
  
  
  
  
  diff = pred1 - pred2
  
  # return all simulated / bootstrapped values if summary is FALSE
  if(!summary){
    return(list(pred1 = pred1, pred2 = pred2, dc = diff))
  }
  
  confint_lower = (1 - conf.int)/2
  result = matrix(nrow=(n+1),ncol=9)
  colnames(result) = c("Mean1", paste0("1:", 100 * confint_lower, "%"), paste0("1:", 100 * (1 - confint_lower), "%"), 
                       "Mean2", paste0("2:", 100 * confint_lower, "%"), paste0("2:", 100 * (1 - confint_lower), "%"), 
                       "Mean.Diff", paste0("diff:", 100 * confint_lower, "%"), paste0("diff:", 100 * (1 - confint_lower), "%"))
  rownames(result) = choices
  
  for(j in 1:(n+1)){
    result[j,] = c(mean(pred1[j,]), quantile(pred1[j,], probs = c(confint_lower, 1 - confint_lower)), 
                   mean(pred2[j,]), quantile(pred2[j,], probs = c(confint_lower, 1 - confint_lower)),
                   mean(diff[j,]), quantile(diff[j,], probs = c(confint_lower, 1 - confint_lower)))
  }
  result
}
