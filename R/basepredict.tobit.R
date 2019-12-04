basepredict.tobit = function(model, values, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL){
    
    n = sim.count
    mu = coef(model)
    if(is.null(sigma)){
      sigma = vcov(model)
    }
    sigma = sigma[-nrow(sigma), -nrow(sigma)]
    if(!is.null(set.seed)){
      set.seed(set.seed)
    }
    sim = MASS::mvrnorm(n, mu, sigma)
    size = length(values)
    
    v = values
    ev = rep(NA,n)
    
    
    for(i in 1:n){
      ev[i] = sum(sim[i,]%*%v)
    }
    
    ev = ev * pnorm(ev / model$scale) 
    
    results = t(as.matrix(c(mean(ev,na.rm=T),quantile(ev,(1-conf.int)/2,na.rm=T),quantile(ev,conf.int+(1-conf.int)/2,na.rm=T))))
    
    colnames(results) = c("Mean",paste0(100*((1-conf.int)/2),"%"),paste0(100*(conf.int+(1-conf.int)/2),"%"))
    
    return(results)
}
