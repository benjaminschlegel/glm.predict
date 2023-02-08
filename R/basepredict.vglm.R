basepredict.vglm = function(model, values, sim.count = 1000, conf.int = 0.95, sigma=NULL, set.seed=NULL,
                            type = c("any", "simulation", "bootstrap"), summary = TRUE){
  
  # check inputs
  if(sum("vglm" %in% class(model)) == 0){
    stop("model has to be of type vglm()")
  }
  if(!("cumulative" %in% model@family@vfamily)){
    stop("only family cumulative is currently supported by glm.predict for vglm() models")
  }
  if(!is.numeric(sim.count) | round(sim.count) != sim.count){
    stop("sim.count has to be whole number")
  }
  if(!is.numeric(conf.int)){
    stop("conf.int has to be numeric")
  }
  if(!is.null(set.seed) & !is.numeric(set.seed)){
    stop("set.seed must be numeric")
  }
  
  type = match.arg(type)
  
  if(type == "any"){
    if(nrow(model.frame(model)) < 500){
      type = "bootstrap"
      message("Type not specified: Using bootstrap as n < 500")
    }else{
      type = "simulation"
      message("Type not specified: Using simulation as n >= 500")
    }
  }
  
  # initialize variables
  l = length(values)
  if(is.null(sigma)){
    sigma = stats::vcov(model)
  }
  if(nrow(sigma) != length(model@coefficients)){
    warning("sigma and coef/zeta do not match, ignoring the specified sigma")
    sigma = stats::vcov(model)
  }
  level.count = length(model@extra$colnames.y)
  kappa.count = level.count - 1
  
  x = values
  estim = model@coefficients
  
  kappa = list()
  betaByKappa = list()
  for(i in 1:kappa.count){
    kappa[[length(kappa)+1]] = matrix(NA, nrow=sim.count, ncol=1, byrow=TRUE)
  }  
  pred = matrix(NA,nrow = sim.count, ncol = level.count,byrow=TRUE)
  
  # simulation
  if(!is.null(set.seed)){
    set.seed(set.seed)
  }
  
  if(type == "simulation"){
    estim_draw = MASS::mvrnorm(sim.count, estim, sigma)
  }else{ # bootstrap
    boot = function(x, model){
      data = model.frame(model)
      colnames(data)[1] = gsub("ordered\\(", "", colnames(data)[1])
      colnames(data)[1] = gsub("\\)", "", colnames(data)[1])
      sample_data = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      model_updated = update(model, data = sample_data)
      model_updated@coefficients
    }
    estim_draw = do.call('rbind', lapply(seq_len(sim.count), boot, model))
  }
  beta_draw = estim_draw[,level.count:ncol(estim_draw)]
  for(i in 1:kappa.count){
    byLevelCols = grep(":[1-9]", colnames(beta_draw))
    cols = sort(c(byLevelCols[1:length(byLevelCols)-kappa.count %% kappa.count == i],
                seq_along(colnames(beta_draw))[-byLevelCols]))
    betaByKappa[[i]] = beta_draw[,cols]
    kappa[[i]][,] = estim_draw[,i]
  }
  
  if(length(values) != ncol(betaByKappa[[1]])){
    print(values)
    print(ncol(betaByKappa[[1]]))
    stop("the length of values is not identical to the number of coefficient of the model")
  }
  
  if(is.null(dim(beta_draw))){
    beta_draw = as.matrix(beta_draw)
  }
  for(j in 1:level.count){
    if(j == 1){
      pred[,j] = exp(kappa[[j]] - betaByKappa[[j]] %*% x) / (1 + exp(kappa[[j]]  - betaByKappa[[j]]  %*% x))
    }else if(j == level.count){
      pred[,j] = 1 / (1 + exp(kappa[[j-1]] - betaByKappa[[j-1]] %*% x))
    }else{
      pred[,j] = exp(kappa[[j]] - betaByKappa[[j]] %*% x) / (1 + exp(kappa[[j]] - betaByKappa[[j]] %*% x)) -
        exp(kappa[[j-1]] - betaByKappa[[j-1]] %*% x) / (1 + exp(kappa[[j-1]] - betaByKappa[[j-1]] %*% x))
    }
  }
  
  # return all simulated / bootstrapped values if summary is FALSE
  if(!summary){
    return(pred)
  }
  
  # prepare the results
  confint_lower = (1 - conf.int) / 2
  result = matrix(NA, nrow = level.count, ncol=3)
  for(i in 1:level.count){
    result[i,] = c(mean(pred[,i]), quantile(pred[,i], prob = c(confint_lower, 1 - confint_lower)))
  }
  colnames(result) = c("mean",paste0(100 * confint_lower,"%"),paste0(100 * (1 - confint_lower),"%"))
  rownames(result) = model@extra$colnames.y
  
  return(result)
}
