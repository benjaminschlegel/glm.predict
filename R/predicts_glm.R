glm.predicts = function(model, values, sim.count=1000, conf.int=0.95, sigma=NULL, set.seed=NULL){
  .Deprecated("predicts", msg="Use 'predicts' instead")
  predicts(model, values, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, doPar = FALSE)
}