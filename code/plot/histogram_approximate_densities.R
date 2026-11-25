#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading relevant packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Declaring file paths 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

main_folder <- "/home/gerardo/Dropbox/FW-simulations"


path_density_models <- file.path(main_folder, "functions/density_models.R")

path_simulations <- file.path(main_folder, "results")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading function files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(path_density_models)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Listing simulation files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# List all CSV files starting with "trajectory"
files <- list.files(
  path = path_simulations, 
  pattern = "^trajectory.*\\.csv$", 
  full.names = TRUE
  )

metadata <- tibble(file = files) %>%
  mutate(
    fname = basename(file),
    x0 = as.numeric(str_match(fname, "^trajectory_([^_]+)_")[,2]),
    N  = as.numeric(str_match(fname, "_([^_]+)\\.csv$")[,2])
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation and plotting parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

selected_frequencies <- c(0.1, 0.3, 0.5)

selected_population_size <- 1000

selected_generations <- c(0.03, 0.05, 0.1)*selected_population_size

params_grid <- expand.grid(
  initial_frequency = selected_frequencies,
  generation = selected_generations
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading simulation files and shaping into long format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_traj <- metadata %>%
  filter(x0 %in% selected_frequencies) %>%
  filter(N == selected_population_size) %>%
  mutate(data = map(file, read.csv)) %>%
  select(-c(file, fname)) %>%
  unnest(data)

trajectories_long <- all_traj %>%
  pivot_longer(
    cols = starts_with("t_"),
    names_to = "generation",
    names_prefix = "t_",
    values_to = "value"
  ) %>%
  rename(initial_frequency = x0) %>%
  select(-N)

rm(all_traj)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adaptive beta kernel density estimator
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

adaptive_beta_df <- lapply(seq_len(nrow(params_grid)), function(i) {
  print(i)
  f <- params_grid$initial_frequency[i]
  t <- params_grid$generation[i]
  
  x <- trajectories_long %>%
    filter(initial_frequency == f,
           generation == t) %>%
    pull(value)
  
  dens_obj <- adaptive_beta_kernel_estimator(
    x = x,
    delta = 0.01,
    C1 = 1,
  )
  
  cat("Analyzing x0 =", f, "at time t =", t, "\n")
  
  df <- data.frame(
    x = dens_obj$x,
    density = dens_obj$y,
    initial_frequency = f,
    generation = t
  )
  
  df
}) %>%
  bind_rows # Combine into a single data frame

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Small-time approximation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute asymptotic expansions for all combinations
sta_df <- lapply(seq_len(nrow(params_grid)), function(i) {
  f <- params_grid$initial_frequency[i]
  t <- params_grid$generation[i]
  
  dens_obj <- st_approximation(
    x0 = f,
    N  = selected_population_size,
    t  = t,
    eps = 0.0001
  )
  
  df <- data.frame(
    x = dens_obj$grid,
    density = dens_obj$density,
    initial_frequency = f,
    generation = t
  )

  df
}) %>%
  bind_rows # Combine into a single data frame

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Small-time Gaussian approximation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute asymptotic expansions for all combinations
stga_df <- lapply(seq_len(nrow(params_grid)), function(i) {
  f <- params_grid$initial_frequency[i]
  t <- params_grid$generation[i]

  df <- st_gaussian_approximation_grid(
    x0 = f,
    N = selected_population_size,
    t = t
  )

  # Add columns for frequency and generation
  df$initial_frequency <- f
  df$generation <- t
  df
}) %>% 
  bind_rows()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Normal approximation calculation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute normal expansions for all combinations
normal_df <- lapply(seq_len(nrow(params_grid)), function(i) {
  f <- params_grid$initial_frequency[i]
  t <- params_grid$generation[i]

  df <- normal_wf_grid(
    x0 = f,
    N = selected_population_size,
    t = t
  )

  # Add columns for frequency and generation
  df$initial_frequency <- f
  df$generation <- t
  df
}) %>%
  bind_rows() # Combine into a single data frame

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta distribution calculation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute normal expansions for all combinations
beta_df <- lapply(seq_len(nrow(params_grid)), function(i) {
  f <- params_grid$initial_frequency[i]
  t <- params_grid$generation[i]

  df <- beta_approximation_grid(
    x0 = f,
    N = selected_population_size,
    t = t
  )

  # Add columns for frequency and generation
  df$initial_frequency <- f
  df$generation <- t
  df
}) %>%
  bind_rows # Combine into a single data frame

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p_hist <- trajectories_long %>%
  filter(initial_frequency %in% selected_frequencies,
         generation %in% selected_generations) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1/100)

y_max <- ggplot_build(p_hist)$data[[1]]$density %>%
  max(na.rm = TRUE)

colors <- c(
  "ABKE" = "black",                   
  "STA" = "#E64B35FF",   
  "STGA" = "#3C5488FF", 
  "Gaussian approximation" = "#00A087FF",    
  "Beta approximation" = "#F39B7FFF"    
)

plot_densities <- trajectories_long %>%
  filter(initial_frequency %in% selected_frequencies,
         generation %in% selected_generations) %>%
  mutate(
    generation = as.numeric(generation),
    generation = factor(
      generation,
      levels = sort(unique(generation))
    )
  ) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1/100,
                 fill = "#EFE7DA",
                 color = "black") +
  geom_line(data = adaptive_beta_df,
            aes(x = x, y = density, color = "ADE"),
            linewidth = 1) +
  geom_line(data = sta_df,
            aes(x = x, y = density, color = "STA"),
            linewidth = 1) +
  geom_line(data = stga_df,
            aes(x = x, y = density, color = "STGA"),
            linewidth = 1) +
  geom_line(data = normal_df,
            aes(x = x, y = density, color = "Gaussian approximation"),
            linewidth = 1) +
  geom_line(data = beta_df,
            aes(x = x, y = density, color = "Beta approximation"),
            linewidth = 1) +
  scale_color_manual(
    values = colors,
    breaks = c(
      "ADE",
      "STA",
      "STGA",
      "Beta approximation",
      "Gaussian approximation"
    ),
    guide = guide_legend(nrow = 2, byrow = T)
  ) + 
  labs(
    x = "x",
    y = "Density",
    color = NULL
  ) +
  facet_grid(
    generation ~ initial_frequency,
    labeller = labeller(
      generation = function(x) {
        formatC(as.numeric(x) / (selected_population_size), format = "f", digits = 2)
      }
    )
  ) +
  coord_cartesian(ylim = c(0, 10)) +
  theme_classic(base_size = 12) +
  theme(
    legend.key.width = unit(1.5, "lines"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    plot.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(t = -10),   # reduce space above legend
    axis.title.x = element_text(margin = margin(t = 5))
  )

trajectories_long %>%
  filter(initial_frequency == 0.1,
         generation == 30) %>%
  mutate(
    generation = as.numeric(generation),
    generation = factor(
      generation,
      levels = sort(unique(generation))
    )
  ) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1/100,
                 fill = "#EFE7DA",
                 color = "black") +
  geom_line(data = adaptive_beta_df %>%
              filter(initial_frequency == 0.1,
                     generation == 30),
            aes(x = x, y = density, color = "ADE"),
            linewidth = 1) +
  theme_classic()

