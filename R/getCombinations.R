getCombinations = function(matrix, base.combinations, model, dv_levels){
  c = 1
  ncol = ncol(matrix)
  if(inherits(model,"mlogit")){
    cnames = colnames(matrix)
    for(choice in dv_levels){
      cnames = gsub(paste0(":", choice), "", cnames)
    }
    ncol = length(unique(cnames))
  }
  
  if(inherits(model,"vglm")){
    cnames = colnames(matrix)
    cnames = gsub(":[0-9]+", "", cnames)
    ncol = length(unique(cnames)) - 1
  }
  
  if(!inherits(model,"polr") & !inherits(model,"vglm")){
    combinations = matrix(NA, nrow = nrow(base.combinations), ncol = ncol)
    combinations[,1] = 1
    c = c + 1
  }else if(inherits(model,"vglm")){
    combinations = matrix(NA, nrow = nrow(base.combinations), ncol = ncol)
  }else{
    combinations = matrix(NA, nrow = nrow(base.combinations), ncol = ncol - 1)
  }
  combinations[,c : (c + ncol(base.combinations) - 1)] = base.combinations
  c = c + ncol(base.combinations)
  
  # add interactions and polygons
  if(ncol > ncol(base.combinations) + 1){ # add interactions
    cnames = unique(gsub(":[0-9]+", "",colnames(matrix)))
    for(name in cnames){
      if(grepl(":", name)){
        parts = unlist(strsplit(name, ":"))
        if(dim(base.combinations)[1]>1){
          combinations[,c] = apply(base.combinations[,parts], 1, prod)
        }else{
          combinations[,c] = prod(base.combinations[,parts])
        }
        c = c + 1
      }else if(grepl("\\^[0-9]", name)){
        temp.name = gsub("I\\(","",name)
        temp.name = gsub("\\)","",temp.name)
        parts = unlist(strsplit(temp.name, "\\^"))
        combinations[,c] = base.combinations[,parts[1]] ^ as.numeric(parts[2])
        c = c + 1
      }
    }
  }
  combinations
}