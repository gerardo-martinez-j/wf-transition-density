#===============================================================================
# Simulate Wright-Fisher allele frequency trajectories
#===============================================================================
#
# Purpose:
#   Generate Wright-Fisher allele frequency trajectories for different
#   combinations of effective population size and initial allele frequency.
#
# Outputs:
#   CSV files containing simulated allele frequencies at specified
#   generations, saved in the results/ directory.
#
#===============================================================================

library(here)

#===============================================================================
# Declaring file paths
#===============================================================================

path_wf_functions <- here("functions", "wf_discrete.R")
path_results_folder <- here("results")

# Create results folder if it doesn't exist
dir.create(
  path_results_folder,
  showWarnings = FALSE,
  recursive = TRUE
)

#===============================================================================
# Loading function to simulate discrete W-F process
#===============================================================================

source(path_wf_functions)

#===============================================================================
# Setting the parameters
#===============================================================================

params <- expand.grid(
  N  = c(100, 1000),
  x0 = c(0.1, 0.3, 0.5)
)

n_sim <- 1e5

#===============================================================================
# Running the simulations
#===============================================================================

walk(seq_len(nrow(params)), function(i) {
  
  N  <- params$N[i]
  x0 <- params$x0[i]
  
  generations <- c(0.05, 0.1, 0.5) * N
  
  df <- wf_trajectories(
    n_sim = n_sim,
    seed = 17,
    N = N,
    x0 = x0,
    generations = generations
  )
  
  write.csv(
    df,
    file.path(
      path_results_folder,
      paste0("trajectory_", x0, "_", N, ".csv")
    ),
    row.names = FALSE
  )
})