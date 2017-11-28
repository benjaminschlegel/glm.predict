getFactorCombinations = function(n){
  result = matrix(NA,nrow=n*(n-1)/2,ncol=2)
  result[,1] = rep(1:(n-1),(n-1):1)
  result[,2] = unlist(sapply(2:n,seq,to=n))
  return(result)
}