getDummies = function(n){
  result = matrix(0, ncol = n-1, nrow=n)
  result[2:n,] = diag(n-1)
  return(result)
}

getCombinations = function(n){
  result = matrix(NA,nrow=n*(n-1)/2,ncol=2)
  result[,1] = rep(1:(n-1),(n-1):1)
  result[,2] = unlist(sapply(2:n,seq,to=n))
  return(result)
}
