discrete.change <-
function(model,values1,values2,sim.count=1000,conf.int=0.95,sigma=NULL,set.seed=NULL){
  
  # model type
  model.type = family(model)
  link = model.type[2]  
  
  n = sim.count
  mu = coef(model)
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  if(!is.null(set.seed)){
    set.seed(set.seed)
  }
  sim = MASS::mvrnorm(n, mu, sigma)
  size = length(values1)
  
  v = matrix(1:size,size)
  v = v[,-1]
  
  v = cbind(values1,values2)
  ev = cbind(rep(NA,n),rep(NA,n))

  ev[,1] = sapply(1:n,simu.glm,j=1,sim=sim,v=v,link=link)
  ev[,2] = sapply(1:n,simu.glm,j=2,sim=sim,v=v,link=link)

  diff = matrix(1:n,n)
  diff = ev[,1]-ev[,2]
  
  all = cbind(ev,diff)
  
  results = matrix(1:length(all[1,]),length(all[1,]))
  results = results[,-1]
  
  for(i in 1:length(all[1,])){
    results = cbind(results,c(mean(all[,i],na.rm=T),quantile(all[,i],(1-conf.int)/2,na.rm=T),quantile(all[,i],conf.int+(1-conf.int)/2,na.rm=T)))
  }
  
  results = t(results[1:3,])
  colnames(results) = c("Mean",paste0(100*((1-conf.int)/2),"%"),paste0(100*(conf.int+(1-conf.int)/2),"%"))
  rownames(results) = c("Case 1","Case 2","Difference")
  
  return(results)
}

simu.glm = function(i,j,sim,v,link){
  x = c(NA,NA)
  x[j] = sum(sim[i,]%*%v[,j])
  
  # the inverse link functions
  if(link == "logit"){
    return(exp(x[j])/(1+exp(x[j])))
  }
  if(link == "log"){
    return(exp(x[j]))
  }
  if(link == "identity"){
    return(x[j])
  }
  if(link == "probit"){ 
    return(pnorm(x[j]))
  }
  if(link == "cauchit"){
    return(tan(pi*(x[j]-0.5)))
  }
  if(link == "cloglog"){
    return(exp(-exp(x[j]))*(-1+exp(exp(x[j]))))
  }
  if(link == "sqrt"){
    return(x[j]*x[j])
  }
  if(link == "1/mu^2"){
    return(1/sqrt(x[j]))
  }
  if(link == "inverse"){
    return(1/x[j])
  }  
}
