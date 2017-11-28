basepredict = function(model, values, sim.count = 1000, conf.int = 0.95, sigma=NULL, set.seed=NULL){
  UseMethod("basepredict")
}