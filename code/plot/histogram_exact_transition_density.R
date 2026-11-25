#===============================================================================
# Compare simulated Wright-Fisher distributions with exact transition densities
#===============================================================================
#
# Purpose:
#   Compare the empirical allele frequency distributions obtained from
#   Wright-Fisher simulations with the corresponding exact transition
#   densities, for different population sizes and generations.
#
#   Histograms represent the simulated allele frequency distributions,
#   while overlaid curves represent the exact transition densities.
#
# Prerequisites:
#   The simulation trajectories and exact transition densities must be
#   calculated before running this script. Run the following scripts first:
#
#     1. code/simulation/simulation_discrete_wf.R
#     2. code/simulation_exact_transition_density.R
#
#   These scripts generate the CSV files required by the present script
#   and save them in the results/ directory.
#
# Outputs:
#   Plots comparing simulated allele frequency distributions with exact
#   transition densities for population sizes N = 100 and N = 1000.
#
#===============================================================================

#===============================================================================
# Checking prerequisites
#===============================================================================

required_files <- c(
  "trajectory_0.5_100.csv",
  "trajectory_0.5_1000.csv",
  "exact_density.csv"
)

missing_files <- required_files[
  !file.exists(here("results", required_files))
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Required files are missing from the results/ directory:\n",
      paste0("  - ", missing_files, collapse = "\n"),
      "\n\nPlease run the following scripts before proceeding:\n",
      "  1. code/simulation/simulation_discrete_wf.R\n",
      "  2. code/simulation_exact_transition_density.R"
    ),
    call. = FALSE
  )
}

#===============================================================================
# Loading relevant packages
#===============================================================================

library(here)
library(tidyverse)

#===============================================================================
# Declaring file paths 
#===============================================================================

path_density_models <- here("functions", "density_models.R")
path_results_folder <- here("results")

#===============================================================================
# Loading density functions
#===============================================================================

source(path_density_models)

#===============================================================================
# Loading simulation file
#===============================================================================

trajectory_file_names <- c("trajectory_0.5_100.csv","trajectory_0.5_1000.csv")

trajectory <- map_dfr(
  trajectory_file_names,
  function(f) {
    
    x0 <- str_match(f, "^trajectory_([^_]+)_")[,2] |> as.numeric()
    N  <- str_match(f, "_([^_]+)\\.csv$")[,2]  |> as.numeric()
    
    read.csv(file.path(path_results_folder, f)) |>
      mutate(
        x0 = x0,
        pop_size = N
      )
  }
) |>
  pivot_longer(
    cols = starts_with("t_"),
    names_to = "generation",
    names_prefix = "t_",
    values_to = "value"
  ) |>
  drop_na() |>
  mutate(
    generation = as.numeric(generation),
    generation = generation / pop_size )

gens_to_plot <- c(0.05, 0.1, 0.5)

trajectory <- trajectory |>
  filter(generation %in% gens_to_plot)

#===============================================================================
# Loading pre-calculated exact density
#===============================================================================

exact_density_df <- read.csv(
  file.path(path_results_folder, "exact_density.csv")
) |>
  filter(!(x %in% c(0.000001000, 0.9999990)))

#===============================================================================
# Plotting histogram with overlaid exact transition density
#===============================================================================

#===============================================================================
# pop_size = 100
#===============================================================================

trajectory |>
  filter(pop_size == 100) |>
  mutate(
    generation = factor(
      generation,
      levels = sort(unique(generation))
    )
  ) |>
  ggplot(aes(x = value)) +
  
  ## Histogram from simulations
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 1 / 50,
    fill = "#FFFBFF",
    color = "black"
  ) +
  
  ## Exact transition density
  geom_line(
    data = exact_density_df,
    aes(x = x, y = density, color = "Exact transition density"),
    linewidth = 1.1
  ) +
  
  scale_color_manual(
    values = c("Exact transition density" = "#006E90")
  ) +
  
  labs(
    x = expression(x[0]),
    y = "Density",
    color = NULL
  ) +
  
  facet_grid(generation ~ ., scales = "free") +
  theme_classic(base_size = 12) +
  theme(
    legend.key.width = unit(1.5, "lines"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

#===============================================================================
# pop_size = 1000
#===============================================================================

levels_gen <- sort(unique(trajectory$generation))

trajectory <- trajectory |>
  mutate(generation = factor(generation, levels = levels_gen))

exact_density_df <- exact_density_df |>
  mutate(generation = factor(generation, levels = levels_gen))

plot_exact_density <- trajectory |>
  filter(pop_size == 1000) |>
  mutate(
    generation = factor(
      generation,
      levels = sort(unique(generation))
    )
  ) |>
  ggplot(aes(x = value)) +
  
  ## Histogram from simulations
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 1 / 25,
    fill = "#EFE7DA",
    color = "black"
  ) +
  
  ## Exact transition density
  geom_line(
    data = exact_density_df,
    aes(x = x, y = density, color = "Exact transition density"),
    linewidth = 1.1
  ) +
  
  scale_color_manual(
    values = c("Exact transition density" = "#006E90")
  ) +
  labs(
    x = expression(x[0]),
    y = "Density",
    color = NULL
  ) +
  facet_wrap(~ generation, scales = "free", ncol = 1) +
  theme_classic(base_size = 12) +
  theme(
    legend.key.width = unit(1.5, "lines"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

