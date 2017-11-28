ordinal_discrete.change = function(model,values1,values2,sim.count = 1000, conf.int = 0.95, sigma=NULL,set.seed=NULL){
  .Deprecated("dc", msg="Use the generic function 'dc' instead")
  dc(model, values = NULL, sim.count = sim.count, conf.int = conf.int, sigma = sigma, set.seed = set.seed, values1 = values1, values2 = values1)
}