dc = function(model, values = NULL, sim.count = 1000, conf.int = 0.95, sigma = NULL, set.seed = NULL, values1 = NULL, values2 = NULL,
              type = c("any", "simulation", "bootstrap"), summary = TRUE){
  UseMethod("dc")
}
