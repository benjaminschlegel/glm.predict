getBaseCombinations = function(data, matrix, values, model, dv_levels = NULL, position = NULL){
  result = initialize_data.frame(data, position = position)
  log.pos = result$log.pos # at which positions are logarithms?
  result = result$result
  
  # get values from character
  value = getValues(values, data)
  values.list = value[[1]]
  is.factor = value[[2]]
  
  # get products for combination callculation
  products = getProducts(value, position)
  
  # initialize variables
  n = length(values.list)
  rows = products[length(products)]
  
  # get base combinations and labels for result
  cnames = colnames(matrix)
  value.names = grep("^[^(][^:\\^]*$",cnames, value = T)
  base.combinations = matrix(NA, nrow = rows, ncol = length(value.names))
  colnames(base.combinations) = value.names
  if(!is.null(position)){
    base.combinations_1 = base.combinations_2 = base.combinations
  }
  for(r in 1:rows){
    if(is.null(position)){
      data.frame.position = length(result) - n + 1
    }else{
      data.frame.position = length(result) - n
    }
    
    c = 1
    for(i in 1:n){
      current.product = products[i]
      current.values = values.list[[i]]
      if(i == 1){
        preproduct = 1
      }else{
        preproduct = products[i-1]
      }
      v1 = floor((r - 1) %% current.product / preproduct) + 1
      if(!is.null(position) && i == position && is.factor[i]){
        combinations = getFactorCombinations(length(current.values[,1]))
        f.v1 = combinations[v1,1]
        f.v2 = combinations[v1,2]
        base.combinations_1[r,c:(c + length(current.values[f.v1,]) - 1)] = current.values[f.v1,]
        base.combinations_2[r,c:(c + length(current.values[f.v1,]) - 1)] = current.values[f.v2,]
        
        # labels
        result[r, data.frame.position] = getLabel(data, i, f.v1)
        result[r, data.frame.position + 1] = getLabel(data, i, f.v2)
        data.frame.position = data.frame.position + 2
        c = c + length(current.values[f.v1,])
      }else if(!is.null(position) && i == position){
        v2 = v1 + 1
        base.combinations_1[r,c:(c + length(current.values[v1]) - 1)] = current.values[v1]
        base.combinations_2[r,c:(c + length(current.values[v1]) - 1)] = current.values[v2]
        
        # labels
        result[r,data.frame.position] = current.values[v1]
        result[r,data.frame.position+1] = current.values[v2]
        data.frame.position = data.frame.position + 2
        c = c + length(current.values[v1])
      }else if(is.factor[i]){
        if(is.null(position)){
          base.combinations[r,c:(c + length(current.values[v1,]) - 1)] = current.values[v1,]
        }else{
          base.combinations_1[r,c:(c + length(current.values[v1,]) - 1)] = current.values[v1,]
          base.combinations_2[r,c:(c + length(current.values[v1,]) - 1)] = current.values[v1,]
        }
        c = c + length(current.values[v1,])
        # labels
        pos = 1
        for(p in 1:length(current.values[v1,])){
          if(current.values[v1,][p]==1){
            pos = p + 1
          }
        }
        result[r, data.frame.position] = getLabel(data,i,pos)
        data.frame.position = data.frame.position + 1
      }else{
        if(is.null(position)){
          base.combinations[r, c] = current.values[v1]
        }else{
          base.combinations_1[r, c] = current.values[v1]
          base.combinations_2[r, c] = current.values[v1]
        }
        
        # labels
        if(data.frame.position %in% log.pos){
          result[r, data.frame.position] = exp(current.values[v1])
        }else{
          result[r, data.frame.position] = current.values[v1]
        }
        
        data.frame.position = data.frame.position + 1
        c = c + 1
      }
    }
  }
  if(!is.null(dv_levels)){
    dv_levels.vector = rep(dv_levels, length(result[,1]))
    result = result[rep(row.names(result),length(dv_levels)),]
    result = result[order(rownames(result)),]
    rownames(result) = 1:length(result[,1])
    result$level = dv_levels.vector
  }
  if(is.null(position)){
    return(list(result=result, base.combinations = base.combinations))
  }else{
    return(list(result=result, base.combinations_1 = base.combinations_1, base.combinations_2 = base.combinations_2))
  }
  
}