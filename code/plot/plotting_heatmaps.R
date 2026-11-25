#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading relevant packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Declaring file paths 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

main_folder <- "/home/gerardo/Dropbox/FW-simulations"

path_results <- file.path(main_folder, "results")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

density_values <- read.csv(
  file.path(path_results, "density_distances.csv")
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Transforming data for plotting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

density_values_long <- density_values %>%
  pivot_longer(
    cols = starts_with(c("l2_", "hell_")),
    names_to = "metric_model",
    values_to = "distance"
  ) %>%
  separate(metric_model, into = c("metric", "model"), sep = "_", extra = "merge") %>%
  mutate(
    generation = as.numeric(generation),
    generation = generation/1000,
    model = recode(model,
                   asymptotic = "STA",
                   beta = "Beta Approximation",
                   normal = "Gaussian Approximation",
                   gaussian = "STGA")
  ) %>%
  mutate(
    distance = -log10(distance),
    model = factor(model,
                   levels = c("STA", "STGA", "Beta Approximation", "Gaussian Approximation"))
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================================================================
# Hellinger distance
#===============================================================================

heatmap_hell <- density_values_long  %>%
  filter(metric == "hell") %>%
  ggplot(
    aes(x = initial_frequency,  # x-axis = x0
        y = generation,         # y-axis = generations
        fill = distance)
  ) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = expression(-log[10](Hellinger))) +
  facet_wrap(~model, scales = "free") +
  labs(
    x = "Initial Frequency (x0)",
    y = "Generation"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0))

#===============================================================================
# L2 distance
#===============================================================================

heatmap_l2 <- density_values_long  %>%
  filter(metric == "l2") %>%
  ggplot(
    aes(x = initial_frequency,  # x-axis = x0
        y = generation,         # y-axis = generations
        fill = distance)
  ) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = expression(-log[10](L^2))) +
  facet_wrap(~model, scales = "free") +
  labs(
    x = "Initial Frequency (x0)",
    y = "Generation"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0))
