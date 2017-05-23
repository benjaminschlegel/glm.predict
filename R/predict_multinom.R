multinom.predict = function(model,values,sim.count=1000,conf.int=0.95,sigma=NULL,set.seed=NULL){
  mu = coef(model)
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  n.coefs = length(mu[1,])
  n = length(mu[,1])
  sim = matrix(nrow=sim.count, ncol = n.coefs * n)
  
  ev = matrix(nrow=(n+1),ncol=sim.count)
  
  for(j in 1:n){
    from = (n.coefs*(j-1)+1)
    to = (n.coefs*j)
    if(!is.null(set.seed)){
      set.seed(set.seed)
    }
    sim[,from:to] = MASS::mvrnorm(sim.count,mu=mu[j,],Sigma=sigma[from:to,from:to])
  }
  
  for(i in 1:sim.count){
    
    sim.temp = NULL
    for(j in 1:n){
      from = (n.coefs*(j-1)+1)
      to = (n.coefs*j)
      if(is.null(sim.temp)){
        sim.temp = sim[i,from:to]
      }else{
        sim.temp = rbind(sim.temp,sim[i,from:to])
      }
      
    }
    
    x = c(0,sim.temp %*% values)
    
    e = exp(x)
    
    for(j in 1:(n+1)){
      ev[j,i] = e[j]/sum(e)
    }
    
  }
  
  upper = conf.int + (1 - conf.int)/2
  lower = (1 - conf.int)/2
  result = matrix(nrow=(n+1),ncol=3)
  colnames(result) = c("mean", paste0(100 * lower, "%"), 
                       paste0(100 * upper, "%"))
  rownames(result) = model$lev
  
  for(j in 1:(n+1)){
    result[j,] = c(mean(ev[j,]), quantile(ev[j,],probs = c(lower,upper)))
  }
  
  return(result)
}
