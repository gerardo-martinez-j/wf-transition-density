#===============================================================================
# Calculate exact Wright-Fisher transition densities
#===============================================================================
#
# Purpose:
#   Calculate the exact transition density of the Wright-Fisher model for
#   a given initial allele frequency and a set of specified generations.
#   The density is evaluated over a discretized allele frequency grid.
#
# Outputs:
#   A CSV file containing the calculated transition densities for each
#   specified generation, saved in the results/ directory.
#
#===============================================================================

library(here)

#===============================================================================
# Declaring file paths 
#===============================================================================

path_density_models <- here("functions", "density_models.R")
path_results_folder <- here("results")

# Create results folder if it doesn't exist
dir.create(
  path_results_folder,
  showWarnings = FALSE,
  recursive = TRUE
)

#===============================================================================
# Loading functions to calculate the exact transition density
#===============================================================================

source(path_density_models)

#===============================================================================
# Calculating exact density
#===============================================================================

gens_to_plot <- c(0.05, 0.1, 0.5)

exact_density_df <- lapply(gens_to_plot, function(t_gen) {
  exact_transition_density_grid(
    x0 = 0.5,
    t = t_gen,
    n_samples = 1e4,
    delta = 0.001,
    n_grid = 500,
    seed = 17
  )
}) |> 
  bind_rows()

#===============================================================================
# Saving file
#===============================================================================

file_name <- file.path(path_results_folder, "exact_density.csv")

write.csv(
  exact_density_df,
  file_name,
  row.names = FALSE
  )