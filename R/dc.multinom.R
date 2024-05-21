dc.multinom = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
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
  if(sum("multinom" %in% class(model)) == 0){
    stop("model has to be of type multinom()")
  }
  if(length(values1) != ncol(coef(model))){
    stop("the length of values1 is not identical to the number of coefficient of the model")
  }
  if(length(values2) != ncol(coef(model))){
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
  
  betas = coef(model)
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  n.coefs = ncol(betas)
  n = nrow(betas)
  
  if(type == "simulation"){
    betas_draw = matrix(nrow = sim.count, ncol = n.coefs * n)
    
    for(j in 1:n){
      from = (n.coefs*(j-1)+1)
      to = (n.coefs*j)
      if(!is.null(set.seed)){
        set.seed(set.seed)
      }
      betas_draw[, from:to] = MASS::mvrnorm(sim.count, betas[j, ], sigma[from:to, from:to])
    }
  }else{ #
    boot = function(x, model){
      data = model.frame(model)
      sample_data = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      if("(weights)" %in% colnames(data)){
        unlist(as.list(t(coef(update(model, data = sample_data, weights = `(weights)`)))))
      }else{
        unlist(as.list(t(coef(update(model, data = sample_data)))))
      }
      
      
    }
    betas_draw = do.call('rbind', lapply(seq_len(sim.count), boot, model))
  }
  
  x = matrix(ncol = 2, nrow = n + 1)
  pred1 = pred2 = matrix(nrow = n + 1, ncol = sim.count)
  
  for(i in 1:sim.count){
    sim.temp = NULL
    for(j in 1:n){
      from = n.coefs * (j-1) + 1
      to = n.coefs * j
      if(is.null(sim.temp)){
        sim.temp = betas_draw[i, from:to]
      }else{
        sim.temp = rbind(sim.temp, betas_draw[i,from:to])
      }
      
    }
    x[,1] = c(0, sim.temp %*% values1)
    x[,2] = c(0, sim.temp %*% values2)
    
    e = exp(x)
    
    for(j in 1:(n+1)){
      pred1[j,i] = e[j, 1] / colSums(e)[1]
      pred2[j,i] = e[j, 2] / colSums(e)[2]
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
  rownames(result) = model$lev
  
  for(j in 1:(n+1)){
    result[j,] = c(mean(pred1[j,]), quantile(pred1[j,], probs = c(confint_lower, 1 - confint_lower)), 
                   mean(pred2[j,]), quantile(pred2[j,], probs = c(confint_lower, 1 - confint_lower)),
                   mean(diff[j,]), quantile(diff[j,], probs = c(confint_lower, 1 - confint_lower)))
  }
  result
}
