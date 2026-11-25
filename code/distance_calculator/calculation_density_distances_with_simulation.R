#===============================================================================
# Calculate distances between simulated and exact Wright-Fisher transition densities
#===============================================================================
#
# Purpose:
#   Simulate Wright-Fisher allele frequency trajectories for a range of
#   initial allele frequencies and generations, and compare the resulting
#   empirical distributions with the four models defined in the text 
#   of the article.
#
# Outputs:
#   A CSV file (density_distances.csv) containing the calculated distances
#   between simulated and exact transition densities for each initial
#   allele frequency and generation, saved in the results/ directory.
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

path_results_folder <- here("results")

#===============================================================================
# Loading function files
#===============================================================================

source(path_density_models)
source(path_density_distances)

#===============================================================================
# Simulation parameters
#===============================================================================

selected_frequencies <- seq(0.02, 0.98, by = 0.02)

N <- 1000
n_sim <- 1e5
seed <- 17

generations <- seq(0.05, 1, by = 0.05) * N

#===============================================================================
# Grid for density evaluation
#===============================================================================

delta <- 0.01
x <- seq(0.01, 0.99, by = delta)

results_list <- list()

#===============================================================================
# Simulation + Distance Computation
#===============================================================================

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

#===============================================================================
# Combine results
#===============================================================================

results <- bind_rows(results_list)

#===============================================================================
# Saving file
#===============================================================================

file_name <- file.path(path_results_folder, "density_distances.csv")

write.csv(results, file_name, row.names = FALSE)
