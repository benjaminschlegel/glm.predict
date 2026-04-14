basepredict.multinom = function(model,values,sim.count=1000,conf.int=0.95,sigma=NULL,set.seed=NULL,
                                type = c("any", "simulation", "bootstrap"), summary = TRUE){

  # check inputs
  if(sum("multinom" %in% class(model)) == 0){
    stop("model has to be of type multinom()")
  }
  if(length(values) != ncol(coef(model))){
    stop("the length of values is not identical to the number of coefficient of the model")
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


  betas = coef(model)
  if(is.null(sigma)){
    sigma = stats::vcov(model)
  }
  n.coefs = ncol(betas)
  n = nrow(betas)

  if(type == "simulation"){
    # set seed once, not per category (previously reset every iteration,
    # which destroyed independence across response categories)
    if(!is.null(set.seed)){
      set.seed(set.seed)
    }
    betas_draw = matrix(nrow = sim.count, ncol = n.coefs * n)
    for(j in seq_len(n)){
      from = n.coefs * (j - 1) + 1
      to   = n.coefs * j
      betas_draw[, from:to] = MASS::mvrnorm(sim.count, betas[j, ], sigma[from:to, from:to])
    }
  }else{
    boot = function(x, model){
      data = model.frame(model)
      sample_data = data[sample(seq_len(nrow(data)), replace = TRUE), ]
      if("(weights)" %in% colnames(data)){
        w <- sample_data[["(weights)"]]
        unlist(as.list(t(coef(update(model, data = sample_data, weights = w)))))
      }else{
        unlist(as.list(t(coef(update(model, data = sample_data)))))
      }
    }
    betas_draw = do.call('rbind', lapply(seq_len(sim.count), boot, model))
  }

  # Vectorized prediction: one matrix multiply per category instead of
  # a nested loop over sim.count x n with rbind()
  yhat_mat = matrix(0, nrow = n + 1, ncol = sim.count)  # row 1 = reference category (yhat = 0)
  for(j in seq_len(n)){
    from = n.coefs * (j - 1) + 1
    to   = n.coefs * j
    yhat_mat[j + 1, ] = betas_draw[, from:to] %*% values
  }
  e = exp(yhat_mat)
  pred = e / rep(colSums(e), each = n + 1)

  # return all simulated / bootstrapped values if summary is FALSE
  if(!summary){
    return(pred)
  }

  confint_lower = (1 - conf.int) / 2
  quantiles = t(apply(pred, 1, stats::quantile, probs = c(confint_lower, 1 - confint_lower)))
  result = cbind(rowMeans(pred), quantiles)
  colnames(result) = c("mean",
                       paste0(100 * confint_lower, "%"),
                       paste0(100 * (1 - confint_lower), "%"))
  rownames(result) = model$lev
  result
}
