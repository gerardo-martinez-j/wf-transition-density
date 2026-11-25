#===============================================================================
# Calculate density distances for a specified initial allele frequency
#===============================================================================
#
# Purpose:
#   Simulate Wright-Fisher allele frequency trajectories for a specified
#   initial allele frequency and calculate distances between the empirical
#   distributions and the four proposed models across a range of generations.
#
#   The initial allele frequency is supplied as a command-line argument,
#   allowing the script to be run in parallel for different initial
#   frequencies.
#
# Outputs:
#   A CSV file (density_distances_x0_XXX.csv) containing the calculated
#   distances for the specified initial allele frequency across all
#   generations, saved in the results/ directory.
#
#===============================================================================

#===============================================================================
# Loading relevant packages
#===============================================================================

library(here)
library(tidyverse)

#===============================================================================
# Declaring file paths 
#===============================================================================

path_density_models <- here("functions", "density_models.R")
path_density_distances <- here("functions", "density_distances.R")
path_wf_discrete <- here("functions", "wf_discrete.R")
path_results <- here("results")

#===============================================================================
# Loading function files
#===============================================================================

source(path_density_models)
source(path_density_distances)
source(path_wf_discrete)

#===============================================================================
# Command-line argument for x0
#===============================================================================

args <- commandArgs(trailingOnly = TRUE)
x0 <- as.numeric(args[1])
cat("Running simulation for x0 =", x0, "\n")

#===============================================================================
# Simulation parameters
#===============================================================================

#selected_frequencies <- round(seq(0.2, 0.98, by = 0.02), 3)

N <- 1000
n_sim <- 1e5
seed <- 17

generations <- round(seq(0.05, 1, by = 0.05), 2) * N

#===============================================================================
# Grid for density evaluation
#===============================================================================

delta <- 0.01
x <- seq(0.01, 0.99, by = delta)

results_list <- list()

#===============================================================================
# Simulation + Distance Computation
#===============================================================================

initial_time <- Sys.time()

for (x0 in selected_frequencies) {
  
  cat("Simulating x0 =", x0, "\n")
  
  # 1) Simulate Wright–Fisher trajectories
  traj_matrix <- wf_trajectories(
    n_sim = n_sim,
    seed = seed,
    N = N,
    x0 = x0,
    generations = generations
  )
  
  # Convert to long format
  traj_df <- as.data.frame(traj_matrix)
  colnames(traj_df) <- paste0("t_", generations)
  
  trajectories_long <- traj_df %>%
    pivot_longer(
      cols = everything(),
      names_to = "generation",
      names_prefix = "t_",
      values_to = "value"
    ) %>%
    mutate(
      generation = as.numeric(generation),
      initial_frequency = x0
    )
  
  # 2) Compute distances
  for (t_gen in generations) {
    cat("Computing t =", t_gen, "\n")
    sample_vec <- trajectories_long %>%
      filter(generation == t_gen) %>%
      pull(value)
    
    dist_tbl <- compute_distances(
      sample_vec = sample_vec,
      x0 = x0,
      t_gen = t_gen,
      N = N,
      x_grid = x,
      delta = delta
    )
    
    dist_tbl$initial_frequency <- x0
    dist_tbl$generation <- t_gen
    
    results_list[[length(results_list) + 1]] <- dist_tbl
  }
}

print(Sys.time()-initial_time)

#===============================================================================
# Combine results
#===============================================================================

results <- bind_rows(results_list)

#===============================================================================
# Save results
#===============================================================================

output_file <- file.path(
  path_results,
  paste0("density_distances_x0_", sprintf("%.3f", x0), ".csv"))

write.csv(
  results,
  output_file,
  row.names = FALSE
  )
cat("Saved results to", output_file, "\n")
