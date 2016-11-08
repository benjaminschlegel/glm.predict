polr.predict = function(model,values,sim.count = 1000, conf.int = 0.95, sigma=NULL){
  
  # check for correct input
  if(length(values) != length(coefficients(model))){
    stop("values has different length then coef")
  }
  
  # initialize variables
  l = length(values)
  n = sim.count
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  level.count = length(model$lev)
  kappa.count = level.count - 1
  
  x = matrix(values,nrow=l,ncol=1,byrow=T) 
  
  draw = matrix(NA,nrow=n,ncol=l+kappa.count,byrow=T)
  beta = matrix(coef(model),nrow=1,ncol=l,byrow=T)
  zeta = matrix(model$zeta,nrow=1,ncol=kappa.count,byrow=T)
  estim = cbind(beta,zeta)
  b = matrix(NA,nrow=n,ncol=l,byrow=T)

  kappa = list()
  for(i in 1:kappa.count){
    kappa[[length(kappa)+1]] = matrix(NA,nrow=n,ncol=1,byrow=TRUE)
  }  
  ev = matrix(NA,nrow=n,ncol=level.count,byrow=TRUE)
  
  # simulation
  draw[,] = MASS::mvrnorm(n,estim,sigma)
  b[,]<-draw[,1:l]
  for(i in 1:kappa.count){
    kappa[[i]][,] = draw[,l+i]
  }
  
  # calculate the discrete changes
  for (i in 1:n)
  {
    for(j in 1:level.count){
        if(j == 1){
          ev[i,j] = exp(kappa[[j]][i,]-b[i,]%*%x)/(1+exp(kappa[[j]][i,]-b[i,]%*%x))
        }else if(j == level.count){
          ev[i,j] = 1/(1+exp(kappa[[j-1]][i,]-b[i,]%*%x))
        }else{
          ev[i,j] = exp(kappa[[j]][i,]-b[i,]%*%x)/(1+exp(kappa[[j]][i,]-b[i,]%*%x)) -
            exp(kappa[[j-1]][i,]-b[i,]%*%x)/(1+exp(kappa[[j-1]][i,]-b[i,]%*%x))
        }
    }
  }
  
  # prepare the results
  upper = conf.int + (1 - conf.int)/2
  lower = (1 - conf.int)/2
  result = matrix(NA,nrow=level.count,ncol=3,byrow=T)
  for(i in 1:level.count){
    result[i,] = c(mean(ev[,i]),quantile(ev[,i],prob=c(lower,upper)))
  }
  colnames(result) = c("mean",paste0(100*lower,"%"),paste0(100*upper,"%"))
  rownames(result) = model$lev
  
  return(result)
}