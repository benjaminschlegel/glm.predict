multinom.predict = function(model,values,sim.count=1000,conf.int=0.95,sigma=NULL){
  mu = coef(model)
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  n.coefs = length(mu[1,])
  n = length(mu[,1])
  sim = matrix(ncol=n.coefs,nrow=n)
  
  ev = matrix(nrow=(n+1),ncol=sim.count)
  
  for(i in 1:sim.count){
    for(j in 1:n){
      sim[j,] = MASS::mvrnorm(mu=mu[j,],Sigma=sigma[(n.coefs*(j-1)+1):(n.coefs*j),(n.coefs*(j-1)+1):(n.coefs*j)])
    }
    
    x = c(0,sim %*% values)
    
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