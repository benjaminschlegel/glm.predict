dc.polr = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
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
  if(sum("polr" %in% class(model)) == 0){
    stop("model has to be of type polr()")
  }
  if(length(values1) != length(coef(model))){
    stop("the length of values1 is not identical to the number of coefficient of the model")
  }
  if(length(values2) != length(coef(model))){
    stop("the length of values2 is not identical to the number of coefficient of the model")
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
  l = length(values1)
  n = sim.count
  if(is.null(sigma)){
    sigma = vcov(model)
  }
  level.count = length(model$lev)
  kappa.count = level.count - 1
  
  x = list()
  x[[1]] = matrix(values1,nrow=l,ncol=1) 
  x[[2]] = matrix(values2,nrow=l,ncol=1)
  

  estim = c(model$coefficients, model$zeta)

  kappa = list()
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
      data = model.frame(model)
      sample_data = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      model_updated = update(model, data = sample_data)
      c(model_updated$coefficients, model_updated$zeta)
    }
    estim_draw = do.call('rbind', lapply(seq_len(sim.count), boot, model))
  }
  
  beta_draw = estim_draw[,1:l]
  for(i in 1:kappa.count){
    kappa[[i]][,] = estim_draw[,l+i]
  }
  
  # calculate the discrete changes
  for(j in 1:level.count){
    for(k in 1:2){
      if(j == 1){
        pred[, j+(k-1)*level.count] = exp(kappa[[j]] - beta_draw %*% x[[k]]) / (1 + exp(kappa[[j]] - beta_draw %*% x[[k]]))
      }else if(j == level.count){
        pred[, j+(k-1)*level.count] = 1 / (1 + exp(kappa[[j-1]] - beta_draw %*% x[[k]]))
      }else{
        pred[, j+(k-1)*level.count] = exp(kappa[[j]] - beta_draw %*% x[[k]]) / (1 + exp(kappa[[j]] - beta_draw %*% x[[k]])) -
          exp(kappa[[j-1]] - beta_draw %*% x[[k]]) / (1 + exp(kappa[[j-1]] - beta_draw %*% x[[k]]))
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
  colnames(result) = c("Mean1",paste0("1:", 100 * confint_lower, "%"),paste0("1:", 100 * (1 - confint_lower), "%"),
                       "Mean2",paste0("2:", 100 * confint_lower, "%"),paste0("2:", 100 * (1 - confint_lower), "%"),
                       "Mean.Diff",paste0("diff:", 100 * confint_lower, "%"),paste0("diff:", 100 * (1 - confint_lower), "%"))
  rownames(result) = model$lev
  result
}