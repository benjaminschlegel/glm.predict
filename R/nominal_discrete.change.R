nominal_discrete.change = function(model,v1,v2,sim.count=1000,conf.int=0.95,sigma=NULL){
  mu = coef(model)
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  n.coefs = length(mu[1,])
  n = length(mu[,1])
  sim = matrix(ncol=n.coefs,nrow=n)
  
  x = matrix(ncol=2,nrow=n+1)
  
  ev1 = matrix(nrow=(n+1),ncol=sim.count)
  ev2 = matrix(nrow=(n+1),ncol=sim.count)
  
  for(i in 1:sim.count){
    for(j in 1:n){
      
      sim[j,] = MASS::mvrnorm(mu=mu[j,],Sigma=sigma[(n.coefs*(j-1)+1):(n.coefs*j),(n.coefs*(j-1)+1):(n.coefs*j)])
    }
    
    x[,1] = c(0,sim %*% v1)
    x[,2] = c(0,sim %*% v2)
    
    e = exp(x)
    
    for(j in 1:(n+1)){
      ev1[j,i] = e[j,1]/colSums(e)[1]
      ev2[j,i] = e[j,2]/colSums(e)[2]
    }
    
  }
  
  diff = ev1-ev2
  
  upper = conf.int + (1 - conf.int)/2
  lower = (1 - conf.int)/2
  result = matrix(nrow=(n+1),ncol=9)
  colnames(result) = c("Mean1", paste0("1:", 100 * lower, "%"), 
                       paste0("1:", 100 * upper, "%"), "Mean2", paste0("2:", 
                                                                       100 * lower, "%"), paste0("2:", 100 * upper, "%"), 
                       "Mean.Diff", paste0("diff:", 100 * lower, "%"), paste0("diff:", 
                                                                              100 * upper, "%"))
  rownames(result) = model$lev
  
  for(j in 1:(n+1)){
    result[j,] = c(mean(ev1[j,]), quantile(ev1[j,],probs = c(lower,upper)),mean(ev2[j,]), quantile(ev2[j,],probs = c(lower,upper)),
                   mean(diff[j,]),quantile(diff[j,],probs = c(lower,upper)))
  }
  
  return(result)
}