#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading relevant packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation and plotting parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N <- 1000             # population size
t_max <- 1            # max continuous time
n_points <- 1000      # grid points
x0_vals <- c(0.1, 0.3, 0.5)  # initial allele frequencies

# Time grid
t <- seq(0, t_max, length.out = n_points)
n <- t * N  # convert continuous time to discrete generations

plot_data <- data.frame()

for (x0 in x0_vals) {
  
  var_stga <- t * x0 * (1 - x0)                 # small-time approximation
  var_diff <- (1 - exp(-t)) * x0 * (1 - x0)    # diffusion approximation
  var_true <- x0 * (1 - x0) * (1 - (1 - 1/N)^n) # exact discrete variance
  
  df <- data.frame(
    t = rep(t, 3),
    variance = c(var_stga, var_diff, var_true),
    model = rep(c("STGA", "Diffusion approx", "Exact discrete"), each = length(t)),
    x0 = x0
  )
  
  plot_data <- rbind(plot_data, df)
}

# Plot
ggplot(plot_data, aes(x = t, y = variance, color = model, linetype = model)) +
  geom_line(size = 1) +
  facet_wrap(~x0, ncol = 1, labeller = labeller(x0 = function(x) paste0("x0 = ", x))) +
  labs(
    x = "Time t",
    y = expression(Var(X[t])),
    color = "Model",
    linetype = "Model"
  ) +
  theme_classic() +
  scale_color_brewer(palette = "Set1")
