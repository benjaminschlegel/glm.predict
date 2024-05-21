basepredict.multinom = function(model,values,sim.count=1000,conf.int=0.95,sigma=NULL,set.seed=NULL,
                                type = c("any", "simulation", "bootstrap"), summary = TRUE){
  
  # check inputs
  if(sum("multinom" %in% class(model)) == 0){
    stop("model has to be of type multinom()")
  }
  if(length(values) != ncol(coef(model))){
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

  pred = matrix(nrow = n + 1, ncol = sim.count)
  for(i in 1:sim.count){
    sim.temp = NULL
    for(j in 1:n){
      from = (n.coefs * (j-1) + 1)
      to = n.coefs * j
      if(is.null(sim.temp)){
        sim.temp = betas_draw[i, from:to]
      }else{
        sim.temp = rbind(sim.temp, betas_draw[i, from:to])
      }
    }
    
    yhat = c(0, sim.temp %*% values)
    e = exp(yhat)
    
    for(j in 1:(n+1)){
      pred[j, i] = e[j] / sum(e)
    }
  }
  
  # return all simulated / bootstrapped values if summary is FALSE
  if(!summary){
    return(pred)
  }
  
  confint_lower = (1 - conf.int) / 2
  result = matrix(nrow = n + 1, ncol = 3)
  colnames(result) = c("mean", 
                       paste0(100 * confint_lower, "%"), 
                       paste0(100 * (1 - confint_lower), "%"))
  rownames(result) = model$lev
  for(j in 1:(n+1)){
    result[j,] = c(mean(pred[j,]), quantile(pred[j,],probs = c(confint_lower, 1 - confint_lower)))
  }
  result
}
