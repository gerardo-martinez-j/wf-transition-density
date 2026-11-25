source("~/Documentos/Material/Papers/Wright-Fisher/functions/wf_discrete.R")

N <- 1e3
p <- 0.2
final_generation <- 1e3
generations <- c(10, 20, 30, 1e3)

plot_trajectories_wf(
  n_sim = 10, 
  N = N, 
  p = p, 
  final_generation = final_generation)