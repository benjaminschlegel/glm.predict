predicts = function(model, values, position = NULL, sim.count = 1000, conf.int=0.95, sigma = NULL, set.seed = NULL, doPar = FALSE,
                    type = c("any", "simulation", "bootstrap")){
  if(!is.character(values)){
    stop("values must be given as character!")
  }
  if("vglm" %in% class(model)){
    full_data = VGAM::model.frame(model)
  }else{
    full_data =  stats::model.frame(model)
  }
                     
                     
  if(any(c("lmerMod", "glmerMod") %in% class(model))){
    full_data = full_data[,-which(colnames(full_data) %in% names(ranef(model)))]
  }
  
  # collapse values to one character, if given as vector
  if(length(values) > 1){
    values = paste(values, collapse = ";")
  }
  
  if("tobit" %in% class(model)){
    colnames(full_data)[1] = "y"
  }
  
  # reshape mlogit data
  if("dfidx" %in% class(full_data)){ 
    choices = levels(full_data$idx[[2]])
    full_data = as.data.frame(full_data)
    pos_idx = which(colnames(full_data) == "idx")
    full_data = full_data[, -(pos_idx:ncol(full_data))]
    full_data[,1] = as.factor(choices)
    
  }
  
  # remove weights column
  if("(weights)" %in% colnames(full_data)){ 
    full_data = full_data[,-which(colnames(full_data) == "(weights)")]
  }
  
  # remove polynomial values
  full_data = full_data[, grep("^[^(][^:\\^]*$", colnames(full_data), value = TRUE)]
  if(length(unlist(strsplit(values, ";"))) != ncol(full_data) - 1){
    stop("The length of values does not match the number of independend variables.")
  }
  
  if(!is.null(position) && (!is.numeric(position) || position != round(position))){
    stop("position must be a whole number or NULL.")
  }
  
  if(inherits(model, "multinom") && doPar){
    doPar = FALSE
    warning("Parallel version not supported for multinom() models. Setting doPar to FALSE.")
  }
  
  type = match.arg(type)
  
  if(type == "any"){
    if(nrow(full_data) < 500){
      type = "bootstrap"
      message("Type not specified: Using bootstrap as n < 500")
    }else{
      type = "simulation"
      message("Type not specified: Using simulation as n >= 500")
    }
  }
  
  # remove any empty space in values
  values = gsub("\\s","",values)
  
  # get data
  if(inherits(model,"polr") || inherits(model,"multinom") || inherits(model, "mlogit")){
    if(!is.null(levels(full_data[,1]))){
      dv_levels = levels(full_data[,1])
    }else{
      dv_levels = levels(as.factor(full_data[, 1]))
    }
  }else if(inherits(model,"vglm")){
    dv_levels = model@extra$colnames.y
  }else{
    dv_levels = NULL
  }
  # data = full_data[,-1]  # data without y
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
    combinations = getCombinations(matrix, base.combinations, model, dv_levels)
  }else{
    combinations_1 = getCombinations(matrix, base.combinations_1, model, dv_levels)
    combinations_2 = getCombinations(matrix, base.combinations_2, model, dv_levels)
  }

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cores <- 2L
  } else {
    # use all cores in devtools::test()
    cores <- parallel::detectCores()
  }
  
  if(doPar && cores > 1){
    # set up parallel cluster
    cl = parallel::makeCluster(cores)
    
    if(is.null(position)){
      parallel::clusterExport(cl, varlist = c("basepredict.lm","basepredict.glm","basepredict.polr","basepredict.multinom","basepredict.tobit", "calculate_glm_pred", "basepredict.mlogit"), envir=environment())
      parallel::clusterEvalQ(cl, library("MASS"))
      parallel::clusterEvalQ(cl, library("nnet"))
      parallel::clusterEvalQ(cl, library("mlogit"))
      parallel::clusterEvalQ(cl, library("dfidx"))
      
      # simulate
      if(is.null(dv_levels)){
        result[, 1:3] = t(parallel::parApply(cl, combinations, 1, basepredict, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, type = type))
      }else{
        temp = parallel::parApply(cl, combinations, 1, basepredict, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, type = type)
        result[, 1:3] = t(do.call(rbind,lapply(1:3, getResultMatrix, result_matrix = temp, levels = length(dv_levels), base.combinations = base.combinations)))
      }
    }else{
      parallel::clusterExport(cl, varlist = c("dc.lm", "dc.glm","dc.polr","dc.multinom", "calculate_glm_pred", "dc.tobit", "dc.mlogit"), envir=environment())
      parallel::clusterEvalQ(cl, library("MASS"))
      parallel::clusterEvalQ(cl, library("nnet"))
      parallel::clusterEvalQ(cl, library("mlogit"))
      parallel::clusterEvalQ(cl, library("dfidx"))
      
      # simulate
      combinations = cbind(combinations_1,combinations_2)
      if(is.null(dv_levels)){
        result[, 1:9] = t(parallel::parApply(cl, combinations, 1, dc, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, type = type))
        result[,  c("val1_lower", "val1_upper", "val2_mean", "val2_upper", "dc_mean", "dc_lower")] = 
          result[,  c("val2_mean", "dc_mean", "val1_lower", "dc_lower", "val1_upper", "val2_upper")]
      }else{
        temp = parallel::parApply(cl, combinations, 1, dc,  model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, type = type)
        result[, 1:9] = t(do.call(rbind,lapply(1:9, getResultMatrix, result_matrix = temp, levels = length(dv_levels), base.combinations = base.combinations_1)))
      }
    }

    # stop parallel cluster
    parallel::stopCluster(cl)
  }else{
    # simulate
    if(is.null(position)){
      if(is.null(dv_levels)){
        result[, 1:3] = t(apply(combinations, 1, basepredict, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, type = type))
      }else{
        temp = apply(combinations, 1, basepredict, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, type = type)
        result[, 1:3] = t(do.call(rbind,lapply(1:3, getResultMatrix, result_matrix = temp, levels = length(dv_levels), base.combinations = base.combinations)))
      }
      
    }else{
      combinations = cbind(combinations_1,combinations_2)
      if(is.null(dv_levels)){
        result[, 1:9] = t(apply(combinations, 1, dc, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, type = type))
        result[,  c("val1_lower", "val1_upper", "val2_mean", "val2_upper", "dc_mean", "dc_lower")] = 
          result[,  c("val2_mean", "dc_mean", "val1_lower", "dc_lower", "val1_upper", "val2_upper")]
      }else{
        temp = apply(combinations, 1, dc, model = model, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, type = type)
        result[, 1:9] = t(do.call(rbind,lapply(1:9, getResultMatrix, result_matrix = temp, levels = length(dv_levels), base.combinations = base.combinations_1))) 
      }
    }
  }
  
  # return result data.frame
  result
}


