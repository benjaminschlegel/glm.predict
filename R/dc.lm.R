dc.lm = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL){

    if(is.null(values) && (is.null(values1) || is.null(values2))){
      stop("Either values1 and values2 or values has to be specified!")
    }
    if(!is.null(values)){
      l = length(values)
      values1 = values[1 : (l/2)]
      values2 = values[(l/2 + 1) : l]
    }
    
    # model type
    model.type = family(model)
    link = model.type[2]  
    
    n = sim.count
    mu = coef(model)
    if(is.null(sigma)){
      sigma = stats::vcov(model)
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
    
    ev[,1] = sapply(1:n,simu.lm,j=1,sim=sim,v=v,link=link)
    ev[,2] = sapply(1:n,simu.lm,j=2,sim=sim,v=v,link=link)
    
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

simu.lm = function(i,j,sim,v,link){
  x = c(NA,NA)
  x[j] = sum(sim[i,]%*%v[,j])
  return(x[j])
}
