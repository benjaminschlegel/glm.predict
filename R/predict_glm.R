glm.predict_ <-
function(model,values,sim.count=1000,conf.int=0.95,sigma=NULL){
  
  # model type
  model.type = family(model)
  link = model.type[2]  
  
  n = sim.count
  mu = coef(model)
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  sim = MASS::mvrnorm(n, mu, sigma)
  size = length(values)
  
  v = values
  ev = rep(NA,n)
  
  
  for(i in 1:n){
      x = sum(sim[i,]%*%v)

      # the inverse link functions
      if(link == "logit"){
        ev[i] = exp(x)/(1+exp(x))
      }
      if(link == "log"){
        ev[i] = exp(x)
      }
      if(link == "identity"){
        ev[i] = x
      }
      if(link == "probit"){ 
        ev[i] = pnorm(x)
      }
      if(link == "cauchit"){
        ev[i] = tan(pi*(x-0.5))
      }
      if(link == "cloglog"){
        ev[i] = exp(-exp(x))*(-1+exp(exp(x)))
      }
      if(link == "sqrt"){
        ev[i] = x*x
      }
      if(link == "1/mu^2"){
        ev[i] = 1/sqrt(x)
      }
      if(link == "inverse"){
        ev[i] = 1/x
      }
  }
  

  results = t(as.matrix(c(mean(ev,na.rm=T),quantile(ev,(1-conf.int)/2,na.rm=T),quantile(ev,conf.int+(1-conf.int)/2,na.rm=T))))

  colnames(results) = c("Mean",paste0(100*((1-conf.int)/2),"%"),paste0(100*(conf.int+(1-conf.int)/2),"%"))
  
  return(results)
}
