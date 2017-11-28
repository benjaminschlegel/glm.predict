getDummies = function(n){
  result = matrix(0, ncol = n-1, nrow=n)
  result[2:n,] = diag(n-1)
  return(result)
}