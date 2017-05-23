nominal_discrete.change = function(model,values1,values2,sim.count=1000,conf.int=0.95,sigma=NULL,set.seed=NULL){
  mu = coef(model)
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  n.coefs = length(mu[1,])
  n = length(mu[,1])
  sim = matrix(nrow=sim.count, ncol = n.coefs * n)
  
  x = matrix(ncol=2,nrow=n+1)
  
  ev1 = matrix(nrow=(n+1),ncol=sim.count)
  ev2 = matrix(nrow=(n+1),ncol=sim.count)
  
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
    x[,1] = c(0,sim.temp %*% values1)
    x[,2] = c(0,sim.temp %*% values2)
    
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
