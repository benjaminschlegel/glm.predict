discrete.changes = function(model, values, position=1, sim.count=1000, conf.int=0.95, sigma=NULL,set.seed=NULL){
  .Deprecated("predicts", msg=paste0("Use 'predicts' with position = ",position, " instead"))
  predicts(model, values, position = position, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, doPar = FALSE)
}
