getProducts = function(value, position = NULL){
  values.list = value[[1]]
  is.factor = value[[2]]
  n = length(values.list)
  
  result = 1
  results = c()
  
  for(i in 1:n){
    variable = values.list[[i]]
    if(is.list(variable)){
      variable.length = length(variable[[which(sapply(variable, length) > 1)]])
    }else if(is.factor[i]){
      variable.length = length(variable[,1])
    }else{
      variable.length = length(variable)
    }
    if(!is.null(position) && i == position && !is.factor[i]){
      result = result * (variable.length-1)
    }else if(!is.null(position) && i == position & is.factor[i]){
      gauss = (variable.length)*(variable.length-1)/2
      result = result * gauss
    }
    else{
      result = result * variable.length
    }
    results = c(results,result)
  }
  return(results)
}