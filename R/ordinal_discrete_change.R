ordinal_discrete.change = function(model,v1,v2,sim.count = 1000, conf.int = 0.95, sigma=NULL){
  
  # check for correct imput
  if(length(v1) != length(v2) & length(v2) != length(coefficients(model))){
    stop("v1 or v2 have different length then coef")
  }
  
  # initialize variables
  l = length(v1)
  n = sim.count
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  level.count = length(model$lev)
  kappa.count = level.count - 1
  
  x = list()
  x[[1]] = matrix(v1,nrow=l,ncol=1,byrow=T) 
  x[[2]] = matrix(v2,nrow=l,ncol=1,byrow=T)
  
  draw = matrix(NA,nrow=n,ncol=l+kappa.count,byrow=T)
  beta = matrix(coef(model),nrow=1,ncol=l,byrow=T)
  zeta = matrix(model$zeta,nrow=1,ncol=kappa.count,byrow=T)
  estim = cbind(beta,zeta)
  b = matrix(NA,nrow=n,ncol=l,byrow=T)

  kappa = list()
  for(i in 1:kappa.count){
    kappa[[length(kappa)+1]] = matrix(NA,nrow=n,ncol=1,byrow=TRUE)
  }  
  delta = matrix(NA,nrow=n,ncol=level.count,byrow=TRUE)
  ev = matrix(NA,nrow=n,ncol=2*level.count,byrow=TRUE)
  
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
      for(k in 1:2){
        if(j == 1){
          ev[i,j+(k-1)*level.count] = exp(kappa[[j]][i,]-b[i,]%*%x[[k]])/(1+exp(kappa[[j]][i,]-b[i,]%*%x[[k]]))
        }else if(j == level.count){
          ev[i,j+(k-1)*level.count] = 1/exp(kappa[[j-1]][i,]-b[i,]%*%x[[k]])
        }else{
          ev[i,j+(k-1)*level.count] = exp(kappa[[j]][i,]-b[i,]%*%x[[k]])/(1+exp(kappa[[j]][i,]-b[i,]%*%x[[k]])) -
            exp(kappa[[j-1]][i,]-b[i,]%*%x[[k]])/(1+exp(kappa[[j-1]][i,]-b[i,]%*%x[[k]]))
        }
      }
      delta[i,j] = ev[i,j] - ev[i,j+level.count]
    }
  }
  
  # prepare the results
  upper = conf.int + (1 - conf.int)/2
  lower = (1 - conf.int)/2
  result = matrix(NA,nrow=level.count,ncol=9,byrow=T)
  for(i in 1:level.count){
    result[i,] = c(mean(ev[,i]),quantile(ev[,i],prob=c(lower,upper)),
                   mean(ev[,i+level.count]),quantile(ev[,i+level.count],prob=c(lower,upper)),
                   mean(delta[,i]),quantile(delta[,i],prob=c(lower,upper)))
  }
  colnames(result) = c("Mean1",paste0("1:",100*lower,"%"),paste0("1:",100*upper,"%"),"Mean2",paste0("2:",100*lower,"%"),paste0("2:",100*upper,"%"),"Mean.Diff",paste0("diff:",100*lower,"%"),paste0("diff:",100*upper,"%"))
  rownames(result) = model$lev
  
  return(result)
}