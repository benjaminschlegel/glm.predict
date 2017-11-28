polr.predict = function(model,values,sim.count = 1000, conf.int = 0.95, sigma=NULL, set.seed=NULL){
  .Deprecated("basepredict", msg="Use the generic function 'basepredict' instead")
  basepredict(model, values, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed)
}