predicts = function(model, values, position=NULL, sim.count=1000, conf.int=0.95, sigma=NULL, set.seed=NULL, doPar = TRUE){
  if(!is.character(values)){
    stop("values must be given as character!")
  }
  
  if(inherits(model, "multinom")){
    doPar = F
  }
  
  # remove any empty space in values
  values = gsub("\\s","",values)
  
  # get data
  full_data = stats::model.frame(model)
  if(inherits(model,"polr") || inherits(model,"multinom")){
    if(!is.null(levels(full_data[,1]))){
      dv_levels = levels(full_data[,1])
    }else{
      dv_levels = unique(full_data[,1])
    }
  }else{
    dv_levels = NULL
  }
  #data = full_data[,-1]  # data without y
  matrix = stats::model.matrix(model)
  
  # get base combinations
  char_pos = which(sapply(full_data, is.character))
  for(i in char_pos){
    full_data[,i] = as.factor(full_data[,i])
  }
  temp = getBaseCombinations(full_data, matrix, values, model, dv_levels, position)
  result = temp[["result"]]
  if(is.null(position)){
    base.combinations = temp[["base.combinations"]]
  }else{
    base.combinations_1 = temp[["base.combinations_1"]]
    base.combinations_2 = temp[["base.combinations_2"]]
  }
 
  
  # add other things to base combinations
  if(is.null(position)){
    combinations = getCombinations(matrix, base.combinations, model)
  }else{
    combinations_1 = getCombinations(matrix, base.combinations_1, model)
    combinations_2 = getCombinations(matrix, base.combinations_2, model)
  }
  
  
  cores = parallel::detectCores() - 1
  if(doPar && cores > 1){
    # set up parallel cluster
    cl = parallel::makeCluster(cores)
    
    if(is.null(position)){
      parallel::clusterExport(cl, varlist = c("basepredict.lm","basepredict.glm","basepredict.polr","basepredict.multinom","basepredict.tobit"), envir=environment())
      parallel::clusterEvalQ(cl, library("MASS"))
      parallel::clusterEvalQ(cl, library("nnet"))
      
      # simulate
      if(is.null(dv_levels)){
        result[, 1:3] = t(parallel::parApply(cl, combinations, 1, basepredict, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed))
      }else{
        temp = parallel::parApply(cl, combinations, 1, basepredict, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed)
        result[, 1:3] = t(do.call(rbind,lapply(1:3, getResultMatrix, result_matrix = temp, levels = length(dv_levels), base.combinations = base.combinations)))
      }
    }else{
      parallel::clusterExport(cl, varlist = c("dc.lm", "dc.glm","dc.polr","dc.multinom", "simu.glm", "dc.tobit"), envir=environment())
      parallel::clusterEvalQ(cl, library("MASS"))
      parallel::clusterEvalQ(cl, library("nnet"))
      
      # simulate
      combinations = cbind(combinations_1,combinations_2)
      if(is.null(dv_levels)){
        result[, 1:9] = t(parallel::parApply(cl, combinations, 1, dc, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed))
        result[,  c("val1_lower", "val1_upper", "val2_mean", "val2_upper", "dc_mean", "dc_lower")] = 
          result[,  c("val2_mean", "dc_mean", "val1_lower", "dc_lower", "val1_upper", "val2_upper")]
      }else{
        temp = parallel::parApply(cl, combinations, 1, dc,  model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed)
        result[, 1:9] = t(do.call(rbind,lapply(1:9, getResultMatrix, result_matrix = temp, levels = length(dv_levels), base.combinations = base.combinations_1)))
      }
    }

    # stopp parallel cluster
    parallel::stopCluster(cl)
  }else{
    # simulate
    if(is.null(position)){
      if(is.null(dv_levels)){
        result[, 1:3] = t(apply(combinations, 1, basepredict, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed))
      }else{
        temp = apply(combinations, 1, basepredict, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed)
        result[, 1:3] = t(do.call(rbind,lapply(1:3, getResultMatrix, result_matrix = temp, levels = length(dv_levels), base.combinations = base.combinations)))
      }
      
    }else{
      combinations = cbind(combinations_1,combinations_2)
      if(is.null(dv_levels)){
        result[, 1:9] = t(apply(combinations, 1, dc, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed))
        result[,  c("val1_lower", "val1_upper", "val2_mean", "val2_upper", "dc_mean", "dc_lower")] = 
          result[,  c("val2_mean", "dc_mean", "val1_lower", "dc_lower", "val1_upper", "val2_upper")]
      }else{
        temp = apply(combinations, 1, dc, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed)
        result[, 1:9] = t(do.call(rbind,lapply(1:9, getResultMatrix, result_matrix = temp, levels = length(dv_levels), base.combinations = base.combinations_1))) 
      }
    }
  }
  
  # return result data.frame
  result
}


