dc.multinom = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
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
  if(sum("multinom" %in% class(model)) == 0){
    stop("model has to be of type multinom()")
  }
  if(length(values1) != ncol(coef(model))){
    stop("the length of values1 is not identical to the number of coefficient of the model")
  }
  if(length(values2) != ncol(coef(model))){
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

  # Vectorized prediction for both scenarios (row 1 = reference category, yhat = 0)
  yhat1 = matrix(0, nrow = n + 1, ncol = sim.count)
  yhat2 = matrix(0, nrow = n + 1, ncol = sim.count)
  for(j in seq_len(n)){
    from = n.coefs * (j - 1) + 1
    to   = n.coefs * j
    block = betas_draw[, from:to]
    yhat1[j + 1, ] = block %*% values1
    yhat2[j + 1, ] = block %*% values2
  }
  e1 = exp(yhat1)
  e2 = exp(yhat2)
  pred1 = e1 / rep(colSums(e1), each = n + 1)
  pred2 = e2 / rep(colSums(e2), each = n + 1)
  diff  = pred1 - pred2

  # return all simulated / bootstrapped values if summary is FALSE
  if(!summary){
    return(list(pred1 = pred1, pred2 = pred2, dc = diff))
  }

  confint_lower = (1 - conf.int) / 2
  q1 = t(apply(pred1, 1, stats::quantile, probs = c(confint_lower, 1 - confint_lower)))
  q2 = t(apply(pred2, 1, stats::quantile, probs = c(confint_lower, 1 - confint_lower)))
  qd = t(apply(diff,  1, stats::quantile, probs = c(confint_lower, 1 - confint_lower)))
  result = cbind(rowMeans(pred1), q1, rowMeans(pred2), q2, rowMeans(diff), qd)
  colnames(result) = c("Mean1", paste0("1:", 100 * confint_lower, "%"), paste0("1:", 100 * (1 - confint_lower), "%"),
                       "Mean2", paste0("2:", 100 * confint_lower, "%"), paste0("2:", 100 * (1 - confint_lower), "%"),
                       "Mean.Diff", paste0("diff:", 100 * confint_lower, "%"), paste0("diff:", 100 * (1 - confint_lower), "%"))
  rownames(result) = model$lev
  result
}
