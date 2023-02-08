dc.vglm = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
                   type = c("any", "simulation", "bootstrap"), summary = TRUE){
  
  # check inputs
  if(is.null(values) && (is.null(values1) || is.null(values2))){
    stop("Either values1 and values2 or values has to be specified!")
  }
  if(!is.null(values)){
    l = length(values)
    values1 = values[1 : (l/2)]
    values2 = values[(l/2 + 1) : l]
  }
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
    if(nrow(VGAM::model.frame(model)) < 500){
      type = "bootstrap"
      message("Type not specified: Using bootstrap as n < 500")
    }else{
      type = "simulation"
      message("Type not specified: Using simulation as n >= 500")
    }
  }

  
  # initialize variables
  l = length(values1)
  n = sim.count
  if(is.null(sigma)){
    sigma = stats::vcov(model)
  }
  if(nrow(sigma) != length(model@coefficients)){
    warning("sigma and coef/zeta do not match, ignoring the specified sigma")
    sigma = stats::vcov(model)
  }
  level.count = length(model@extra$colnames.y)
  kappa.count = level.count - 1
  
  x = list()
  x[[1]] = matrix(values1,nrow=l,ncol=1) 
  x[[2]] = matrix(values2,nrow=l,ncol=1)
  

  estim = model@coefficients

  kappa = list()
  betaByKappa = list()
  for(i in 1:kappa.count){
    kappa[[length(kappa)+1]] = matrix(NA, nrow = sim.count, ncol=1)
  }  
  delta = matrix(NA, nrow = sim.count, ncol=level.count)
  pred = matrix(NA, nrow = sim.count, ncol= 2 * level.count)
  
  # simulation
  if(!is.null(set.seed)){
    set.seed(set.seed)
  }
  
  
  if(type == "simulation"){
    estim_draw = MASS::mvrnorm(sim.count, estim, sigma)
  }else{ # bootstrap
    boot = function(x, model){
      data = VGAM::model.frame(model)
      colnames(data)[1] = gsub("ordered\\(", "", colnames(data)[1])
      colnames(data)[1] = gsub("\\)", "", colnames(data)[1])
      sample_data = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      model_updated = update(model, data = sample_data)
      model_updated@coefficients
    }
    estim_draw = do.call('rbind', lapply(seq_len(sim.count), boot, model))
  }
  
  beta_draw = estim_draw[,level.count:ncol(estim_draw)]
  beta_draw = beta_draw * -1 #vglm has the wron sign compared to polr
  for(i in 1:kappa.count){
    byLevelCols = grep(":[1-9]", colnames(beta_draw))
    cols = sort(c(byLevelCols[1:length(byLevelCols)-kappa.count %% kappa.count == i],
                  seq_along(colnames(beta_draw))[-byLevelCols]))
    betaByKappa[[i]] = beta_draw[,cols]
    kappa[[i]][,] = estim_draw[,i]
  }
  
  #if(length(values1) != ncol(betaByKappa[[1]])){
  #  stop("the length of values1 is not identical to the number of coefficient of the model")
  #}
  #if(length(values2) != ncol(betaByKappa[[1]])){
  #  stop("the length of values2 is not identical to the number of coefficient of the model")
  #}
  
  # calculate the discrete changes
  for(j in 1:level.count){
    for(k in 1:2){
      if(j == 1){
        pred[, j+(k-1)*level.count] = exp(kappa[[j]] - betaByKappa[[j]] %*% x[[k]]) / (1 + exp(kappa[[j]] - betaByKappa[[j]] %*% x[[k]]))
      }else if(j == level.count){
        pred[, j+(k-1)*level.count] = 1 / (1 + exp(kappa[[j-1]] - betaByKappa[[j-1]] %*% x[[k]]))
      }else{
        pred[, j+(k-1)*level.count] = exp(kappa[[j]] - betaByKappa[[j]] %*% x[[k]]) / (1 + exp(kappa[[j]] - betaByKappa[[j]] %*% x[[k]])) -
          exp(kappa[[j-1]] - betaByKappa[[j-1]] %*% x[[k]]) / (1 + exp(kappa[[j-1]] - betaByKappa[[j-1]] %*% x[[k]]))
      }
    }
    delta[,j] = pred[,j] - pred[,j+level.count]
  }
  
  # prepare the results
  confint_lower = (1 - conf.int)/2
  result = matrix(NA,nrow=level.count,ncol=9)
  for(i in 1:level.count){
    result[i,] = c(mean(pred[,i]), quantile(pred[,i],prob=c(confint_lower, 1 - confint_lower)),
                   mean(pred[,i+level.count]),quantile(pred[,i+level.count],prob=c(confint_lower, 1 - confint_lower)),
                   mean(delta[,i]),quantile(delta[,i],prob=c(confint_lower, 1 - confint_lower)))
  }
  colnames(result) = c("mean1",paste0("val1:", 100 * confint_lower, "%"),paste0("val1:", 100 * (1 - confint_lower), "%"),
                       "mean2",paste0("val2:", 100 * confint_lower, "%"),paste0("val2:", 100 * (1 - confint_lower), "%"),
                       "dc_mean",paste0("dc:", 100 * confint_lower, "%"),paste0("dc:", 100 * (1 - confint_lower), "%"))
  rownames(result) = model@extra$colnames.y
  result
}