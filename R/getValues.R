getValues = function(values, data){
  
  result = list()
  pos = 1
  current.values = NA
  
  values.vector = unlist(strsplit(values,";"))
  is.factor = rep(F,length(values.vector))
  for(value in values.vector){
    varName = colnames(data)[pos+1] # +1 because of y
    var = data[[varName]]
    if(grepl("^mode$",value,ignore.case = TRUE)){ # Mode
      mode = Mode(var,na.rm=T)
      if(is.numeric(mode)){
        current.values = mode
      }else{
        n = length(levels(var))
        dummies = getDummies(n)
        current.values = matrix(dummies[which(levels(var)==mode),],nrow=1)
        is.factor[pos] = T
      }
    } # mode
    else if(grepl("^mean$",value,ignore.case = TRUE)){ # mean
      if(!is.numeric(var)){
        stop("Cannot calculate the mean of a non numeric variable")
      }
      current.values = mean(var, na.rm=T)
    } # mean
    else if(grepl("^median$",value,ignore.case = TRUE)){ # median
      if(!is.numeric(var)){
        stop("Cannot calculate the median of a non numeric variable")
      }
      current.values = median(var, na.rm=T)
    } # median
    else if(grepl("^Q[0-9]+$",value,ignore.case = TRUE)){ # quantile
      n.quantile = as.numeric(unlist(strsplit(value,"[Q\\]")))[2]
      if(!is.numeric(var)){
        stop("Cannot calculate the quantiles of a non numeric variable")
      }
      current.values = quantile(var,probs=seq(from=0,to=1,length.out =n.quantile+1),na.rm = T)
    } # quantile
    else if(grepl("^min$",value,ignore.case = TRUE)){ # min
      if(!is.numeric(var)){
        stop("Cannot calculate the minimum of a non numeric variable")
      }
      current.values = min(var,na.rm = T)
    } # min
    else if(grepl("^max$",value,ignore.case = TRUE)){ # max
      if(!is.numeric(var)){
        stop("Cannot calculate the maximum of a non numeric variable")
      }
      current.values = max(var,na.rm = T)
    } # max
    else if(grepl("^F\\([0-9]+(,[0-9]+)*\\)$",value,ignore.case = TRUE)){ # single factor
      components = as.numeric(unlist(strsplit(value,"[F\\(,\\)]")))
      n = length(levels(var))
      x = components[c(-1,-2)]
      dummies = getDummies(n)
      current.values = matrix(dummies[x,], nrow = length(x))
      is.factor[pos] = T
    } # get specific factor levels
    else if(grepl("^F$",value,ignore.case = TRUE)){ # factor
      n = length(levels(var))
      current.values = getDummies(n)
      is.factor[pos] = T
    } # factor
    else if(grepl("^(-?[0-9]+(\\.[0-9]+)?)-(-?[0-9]+(\\.[0-9]+)?),([0-9]+(\\.[0-9]+)?)$",value)){ # from-to,by
      components = as.numeric(unlist(strsplit(value,"[-,]")))
      i.container = c()
      for(i in 1:length(components)){
        if(is.na(components[i]) || components[i]==""){
          components[i+1] = paste0("-",components[i+1])
          i.container = c(i.container,i)
        }
      }
      if(length(i.container)>0){
        components = components[-i.container]
      }
      components = as.numeric(components)
      current.values = seq(from=components[1],to=components[2],by=components[3])
    } # from-to,by
    else if(grepl("^(-?[0-9]+(\\.[0-9]+)?)-(-?[0-9]+(\\.[0-9]+)?)$",value)){ # from-to
      components = unlist(strsplit(value,"-"))
      i.container = c()
      for(i in 1:length(components)){
        if(is.na(components[i]) || components[i]==""){
          components[i+1] = paste0("-",components[i+1])
          i.container = c(i.container,i)
        }
      }
      if(length(i.container)>0){
        components = components[-i.container]
      }
      components = as.numeric(components)
      current.values = components[1]:components[2]
    } # from-to
    else if(grepl("^(-?[0-9]+(\\.[0-9]+)?)(,-?[0-9]+(\\.[0-9]+)?)*$",value)){ # value1[, value2 [, ...]]
      current.values = as.numeric(unlist(strsplit(value,",")))
    } # value1[, value2 [, ...]
    else if(grepl("^log\\((-?[0-9]+(\\.[0-9]+)?)-(-?[0-9]+(\\.[0-9]+)?),(-?[0-9]+(\\.[0-9]+)?)\\)$",value)){ # from-to,by
      value = gsub("log\\(", "", value)
      value = gsub("\\)", "", value)
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
      current.values = log(seq(from=components[1],to=components[2],by=components[3]))
    } # log(from-to,by)
    else if(grepl("^log\\((-?[0-9]+(\\.[0-9]+)?)-(-?[0-9]+(\\.[0-9]+)?)\\)$",value)){ # from-to
      value = gsub("log\\(", "", value)
      value = gsub("\\)", "", value)
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
      current.values = log(components[1]:components[2])
    } # log(from-to)
    else if(grepl("^log\\((-?[0-9]+(\\.[0-9]+)?)(,-?[0-9]+(\\.[0-9]+)?)*\\)$",value)){ # value1[, value2 [, ...]]
      value = gsub("log\\(", "", value)
      value = gsub("\\)", "", value)
      current.values = log(as.numeric(unlist(strsplit(value,","))))
    } # log(value1[, value2 [, ...]])
    else { # invalid syntax
      message(value)
      stop("values has invalid syntax!")
    } # invalid syntax
    result = c(result,list(current.values))
    pos = pos + 1
  }
  return(list(result,is.factor))
}
