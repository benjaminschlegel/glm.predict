discrete.change <-
function(model,values1,values2,sim.count=1000,conf.int=0.95){
  
  # model type
  model.type = family(model)
  link = model.type[2]  
  
  n = sim.count
  mu = coef(model)
  sigma = vcov(model)
  sim = MASS::mvrnorm(n, mu, sigma)
  size = length(values1)
  
  v = matrix(1:size,size)
  v = v[,-1]
  
  v = cbind(values1,values2)
  ev = cbind(rep(NA,n),rep(NA,n))
  
  
  for(i in 1:n){
    x = c(NA,NA)
    for(j in 1:2){
      x[j] = sum(sim[i,]%*%v[,j])

      # the inverse link functions
      if(link == "logit"){
        ev[i, j] = exp(x[j])/(1+exp(x[j]))
      }
      if(link == "log"){
        ev[i, j] = exp(x[j])
      }
      if(link == "identity"){
        ev[i, j] = x[j]
      }
      if(link == "probit"){ 
        ev[i, j] = pnorm(x[j])
      }
      if(link == "cauchit"){
        ev[i, j] = tan(pi*(x[j]-0.5))
      }
      if(link == "cloglog"){
        ev[i, j] = exp(-exp(x[j]))*(-1+exp(exp(x[j])))
      }
      if(link == "sqrt"){
        ev[i, j] = x[j]*x[j]
      }
      if(link == "1/mu^2"){
        ev[i, j] = 1/sqrt(x[j])
      }
      if(link == "inverse"){
        ev[i,j] = 1/x[j]
      }
    }
  }
  
  diff = matrix(1:n,n)
  diff = ev[,1]-ev[,2]
  
  all = cbind(ev,diff)
  
  results = matrix(1:length(all[1,]),length(all[1,]))
  results = results[,-1]
  
  for(i in 1:length(all[1,])){
    results = cbind(results,c(mean(all[,i]),quantile(all[,i],(1-conf.int)/2),quantile(all[,i],conf.int+(1-conf.int)/2)))
  }
  
  results = t(results[1:3,])
  colnames(results) = c("Mean",paste0(100*((1-conf.int)/2),"%"),paste0(100*(conf.int+(1-conf.int)/2),"%"))
  rownames(results) = c("Case 1","Case 2","Difference")
  
  return(results)
}
