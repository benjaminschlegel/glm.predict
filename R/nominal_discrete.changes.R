nominal_discrete.changes = function(model, values, data, position=1, sim.count=1000, conf.int=0.95, sigma=NULL){
  if(!is.character(values)){
    stop("values must be given as character!")
  }
  
  values = gsub("\\s","",values)
  
  # check if any interaction
  formula = getFormulas(model) # variable names [[1]] and interaction positions[[2]]
  
  value = getValues_nominal(model,values,formula[[1]],data) # values as list [[1]] and positions of factors [[2]]
  
  products = getProducts(value,position)
  
  values.list = value[[1]]
  is.factor = value[[2]]
  n = length(values.list)
  rows = products[length(products)]
  
  result = getNames(formula[[1]],position)
  result = cbind(result,level=NA)
  
  n.levels = length(model$lev)
  levels = model$lev
  
  for(r in 1:rows){
    row.values1 = c(1)
    row.values2 = c(1)
    data.frame.position = 10
    for(i in 1:n){
      current.product = products[i]
      current.values = values.list[[i]]
      if(i == 1){
        preproduct = 1
      }else{
        preproduct = products[i-1]
      }
      v1 = floor((r-1)%%current.product/preproduct)+1
      
      row.from = (r-1)*n.levels + 1
      row.to = r*n.levels
      
      if(i == position & is.factor[i]){
        combinations = getCombinations(length(current.values[,1]))
        f.v1 = combinations[v1,1]
        f.v2 = combinations[v1,2]
        row.values1 = c(row.values1,current.values[f.v1,])
        row.values2 = c(row.values2,current.values[f.v2,])
        
        # labels
        result[row.from:row.to,data.frame.position] = getLabel_nominal(model,formula[[1]][i],f.v1,data)
        result[row.from:row.to,data.frame.position+1] = getLabel_nominal(model,formula[[1]][i],f.v2,data)
        data.frame.position = data.frame.position + 2
      }else if(i==position){
        v2 = v1+1
        row.values1 = c(row.values1,current.values[v1])
        row.values2 = c(row.values2,current.values[v2])
        
        # labels
        result[row.from:row.to,data.frame.position] = current.values[v1]
        result[row.from:row.to,data.frame.position+1] = current.values[v2]
        data.frame.position = data.frame.position + 2
      }else if(is.factor[i]){
        row.values1 = c(row.values1,current.values[v1,])
        row.values2 = c(row.values2,current.values[v1,])
        
        # labels
        pos = 1
        for(p in 1:length(current.values[v1,])){
          if(current.values[v1,][p]==1){
            pos = p+1
          }
        }
        result[row.from:row.to,data.frame.position] = getLabel_nominal(model,formula[[1]][i],pos,data)
        data.frame.position = data.frame.position + 1
      }else{
        row.values1 = c(row.values1,current.values[v1])
        row.values2 = c(row.values2,current.values[v1])
        
        # labels
        result[row.from:row.to,data.frame.position] = current.values[v1]
        data.frame.position = data.frame.position + 1
      }
    }
    
    # interactions
    pos = 1
    pos.row.values = 2
    for(ia in formula[[2]]){
      if(ia==0){ # +
        if(is.factor[pos]){
          n.dummies = length(values.list[[pos]][1,])
          pos.row.values = pos.row.values + n.dummies
        }else{
          pos.row.values = pos.row.values + 1
        }
        pos = pos + 1
      }
      if(ia==1){ # *
        if(is.factor[pos]){
          n.dummies = length(values.list[[pos]][1,])
          term.part1.v1 = row.values1[pos.row.values:(pos.row.values+n.dummies-1)]
          term.part1.v2 = row.values2[pos.row.values:(pos.row.values+n.dummies-1)]
          pos.row.values = pos.row.values + n.dummies
        }else{
          term.part1.v1 = row.values1[pos.row.values]
          term.part1.v2 = row.values2[pos.row.values]
          pos.row.values = pos.row.values + 1
        }
        if(is.factor[pos+1]){
          n.dummies = length(values.list[[pos+1]][1,])
          term.part2.v1 = row.values1[pos.row.values:(pos.row.values+n.dummies-1)]
          term.part2.v2 = row.values2[pos.row.values:(pos.row.values+n.dummies-1)]
          pos.row.values = pos.row.values + n.dummies
        }else{
          term.part2.v1 = row.values1[pos.row.values]
          term.part2.v2 = row.values2[pos.row.values]
          pos.row.values = pos.row.values + 1
        }
        row.values1 = c(row.values1,term.part1.v1*term.part2.v1)
        row.values2 = c(row.values2,term.part1.v2*term.part2.v2)
        pos = pos + 2
      }
    }
    subresult = nominal_discrete.change(model,row.values1,row.values2,sim.count,conf.int,sigma)
    
      for(k in 1:n.levels){
        current.row = (r-1)*n.levels + k
        result[current.row,]$mean1 = subresult[k,1]
        result[current.row,]$mean2 = subresult[k,4]
        result[current.row,]$lower1 = subresult[k,2]
        result[current.row,]$upper1 = subresult[k,3]
        result[current.row,]$lower2 = subresult[k,5]
        result[current.row,]$upper2 = subresult[k,6]
        result[current.row,]$mean.diff = subresult[k,7]
        result[current.row,]$lower.diff = subresult[k,8]
        result[current.row,]$upper.diff = subresult[k,9]
        result[current.row,]$level = levels[k]
      }
    
    #     print(paste0(r,"-v1: ",toString(row.values1)))
    #     print(paste0(r,"-v2: ",toString(row.values2)))
  }
  return(result)
}

getFormulas = function(model){
  formula = formula(model)
  dvs = unlist(strsplit(as.character(formula)[3]," + ", fixed=T))
  ia = c()
  for(dv in dvs){
    if(grepl(":",dv)){
      ia = c(ia,2)
      stop("Formula contains :, works only with * for interaction")
    }else if(grepl("\\*",dv)){
      ia = c(ia,1)
    }else{
      ia = c(ia, 0)
    }
  }
  
  formula = gsub(" \\* ", " + ", as.character(formula)[3])
  formula = gsub(" : ", " + ", formula)
  formula = gsub("\\*", " + ", formula)
  formula = gsub(":", " + ", formula)
  
  # remove dublicates
  temp.formula = unlist(strsplit(formula," + ", fixed=T))
  temp.formula = unique(temp.formula)
  
  return(list(temp.formula,ia))
}

getValues_nominal = function(model,values,formula,data){
  
  result = list()
  pos = 1
  current.values = NA
  use.data.stop.error = length(data)
  
  values.vector = unlist(strsplit(values,";"))
  is.factor = rep(F,length(values.vector))
  for(value in values.vector){
    if(grepl("^mode$",value,ignore.case = TRUE)){ # Mode
      varName = formula[pos]
      if(!is.null(data)){
        data.v = data[,grep(varName,colnames(data),value=T)[1]]
        mode = Mode(data.v,na.rm=T)
        if(is.numeric(mode)){
          current.values = mode
        }else{
          n = length(levels(data.v))
          dummies = getDummies(n)
          current.values = matrix(dummies[which(levels(data.v)==mode),],nrow=1)
          is.factor[pos] = T
        }
      }
    } # mode
    else if(grepl("^mean$",value,ignore.case = TRUE)){ # mean
      varName = formula[pos]
      data.v = data[,varName]
      if(!is.numeric(data.v)){
        stop("Cannot calculate the mean of a non numeric variable")
      }
      current.values = mean(data.v, na.rm=T)
    } # mean
    else if(grepl("^median$",value,ignore.case = TRUE)){ # median
      varName = formula[pos]
      data.v = data[,varName]
      if(!is.numeric(data)){
        stop("Cannot calculate the median of a non numeric variable")
      }
      current.values = median(data.v, na.rm=T)
    } # median
    else if(grepl("^Q[0-9]+$",value,ignore.case = TRUE)){ # quantile
      n.quantile = as.numeric(unlist(strsplit(value,"[Q\\]")))[2]
      varName = formula[pos]
      data.v = data[,varName]
      if(!is.numeric(data.v)){
        stop("Cannot calculate the quantiles of a non numeric variable")
      }
      current.values = quantile(data.v,probs=seq(from=0,to=1,length.out =n.quantile+1),na.rm = T)
    } # quantile
    else if(grepl("^min$",value,ignore.case = TRUE)){ # min
      varName = formula[pos]
      data.v = data[,varName]
      if(!is.numeric(data.v)){
        stop("Cannot calculate the minimum of a non numeric variable")
      }
      current.values = min(data.v,na.rm = T)
    } # min
    else if(grepl("^max$",value,ignore.case = TRUE)){ # max
      varName = formula[pos]
      data.v = data[,varName]
      if(!is.numeric(data.v)){
        stop("Cannot calculate the maximum of a non numeric variable")
      }
      current.values = max(data.v,na.rm = T)
    } # max
	else if(grepl("^F[0-9]+\\([0-9]+\\)$",value,ignore.case = TRUE)){ # single factor (old version)
      components = as.numeric(unlist(strsplit(value,"[F\\(\\)]")))
      n = components[2]
      x = components[3]
      dummies = getDummies(n)
      current.values = matrix(dummies[x,],nrow=1)
      is.factor[pos] = T
    } # single factor value (old version)
    else if(grepl("^F\\([0-9]+\\)$",value,ignore.case = TRUE)){ # single factor
      components = as.numeric(unlist(strsplit(value,"[F\\(\\)]")))
      varName = formula[pos]
      data.v = data[,varName]
      n = length(levels(data.v))
      x = components[3]
      dummies = getDummies(n)
      current.values = matrix(dummies[x,],nrow=1)
      is.factor[pos] = T
    } # single factor value
    else if(grepl("^F[0-9]+$",value,ignore.case = TRUE)){ # factor (old version)
      n = as.numeric(unlist(strsplit(value,"F")))[2]
      current.values = getDummies(n)
      is.factor[pos] = T
    } # factor (old version)
    else if(grepl("^F$",value,ignore.case = TRUE)){ # factor
      varName = formula[pos]
      data.v = data[,varName]
      n = length(levels(data.v))
      current.values = getDummies(n)
      is.factor[pos] = T
    } # factor
    else if(grepl("^(-?[0-9]+(\\.[0-9]+)?)-(-?[0-9]+(\\.[0-9]+)?),(-?[0-9]+(\\.[0-9]+)?)$",value)){ # from-to,by
      components = as.numeric(unlist(strsplit(value,"[-,]")))
		i.container = c()
		for(i in 1:length(components)){
		  if(components[i]==""){
			components[i+1] = paste0("-",components[i+1])
			i.container = c(i.container,i)
		  }
		}
		if(length(i.container)>0){
		  components = components[-i.container]
		}
      current.values = seq(from=components[1],to=components[2],by=components[3])
    } # from-to,by
    else if(grepl("^(-?[0-9]+(\\.[0-9]+)?)-(-?[0-9]+(\\.[0-9]+)?)$",value)){ # from-to
      components = unlist(strsplit(value,"-"))
		i.container = c()
		for(i in 1:length(components)){
		  if(components[i]==""){
			components[i+1] = paste0("-",components[i+1])
			i.container = c(i.container,i)
		  }
		}
		if(length(i.container)>0){
		  components = components[-i.container]
		}
      current.values = components[1]:components[2]
    } # from-to
    else if(grepl("^(-?[0-9]+(\\.[0-9]+)?)(,-?[0-9]+(\\.[0-9]+)?)*$",value)){ # value1[, value2 [, ...]]
      current.values = as.numeric(unlist(strsplit(value,",")))
    } # value1[, value2 [, ...]]
    else { # invalid syntax
      print(value)
      stop("values has invalid syntax!")
    } # invalid syntax
    result = c(result,list(current.values))
    pos = pos + 1
  }
  return(list(result,is.factor))
}

getProducts = function(value,position){
  values.list = value[[1]]
  is.factor = value[[2]]
  n = length(values.list)
  
  result = 1
  results = c()
  
  for(i in 1:n){
    variable = values.list[[i]]
    if(is.factor[i]){
      variable.length = length(variable[,1])
    }else{
      variable.length = length(variable)
    }
    if(i == position & !is.factor[i]){
      result = result * (variable.length-1)
    }else if(i == position & is.factor[i]){
      gauss = (variable.length)*(variable.length-1)/2
      result = result * gauss
    }else{
      result = result * variable.length
    }
    results = c(results,result)
  }
  return(results)
}

getDummies = function(n){
  if(n == 2){
    return(matrix(0:1,nrow=2,ncol=1))
  }
  n = n-1
  I = iterpc::iterpc(2, n, label=c(0,1), order=T, replace=T)
  grid = iterpc::getall(I)
  grid = grid[order(rowSums(grid)),]
  grid = subset(grid,rowSums(grid)<=1)
  grid = grid[,n:1]
  return(grid)
}

getCombinations = function(n){
  I = iterpc::iterpc(n, 2)
  grid = iterpc::getall(I)
  return(grid)
}

getNames = function(names,position){
  new.names = c("mean1","mean2","lower1","upper1","lower2","upper2","mean.diff","lower.diff","upper.diff")
  for(i in 1:length(names)){
    if(i != position){
      new.names = c(new.names,names[i])
    }else{
      new.names = c(new.names,paste0(names[i],".1"),paste0(names[i],".2"))
    }
  }
  result = data.frame(t(rep(NA,length(new.names))))
  colnames(result) = new.names
  return(result)
}

getLabel_nominal = function(model,varName,pos,data){
  data.v = data[,grep(varName,colnames(data),value=T)[1]]
  labels = levels(data.v)
  return(labels[pos])
}

Mode = function(x,na.rm=FALSE) {
  if(na.rm){
    x = na.omit(x)
  }
  ux = unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}
