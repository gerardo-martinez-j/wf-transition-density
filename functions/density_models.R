#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading relevant packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(e1071)     # rbridge
library(Bolstad2)  # sintegral
library(pracma) #cot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exact transition density functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Local variance transformation
#'
#' Function to obtain a constant local variance term in the Wright-Fisher
#' diffusion.
#'
#' @param x A numeric vector of allele frequencies in \eqn{[0,1]}.
#'
#' @return The transformed value for the transformed stochastic process.
#' @export

f_trans <- function(x){
  return(2*asin(sqrt(x)))
}

#' Heat kernel for Brownian motion
#'
#' Evaluates the Gaussian transition density of Brownian motion.
#'
#' @param y Evaluation point.
#' @param y0 Initial point.
#' @param t Time parameter.
#'
#' @return Numeric value of the heat kernel.
#' @export
qt <- function(y, y0, t) {
  result <- 1/sqrt(2*pi*t)*exp(-(y0-y)^2/(2*t))
  return(result)
}

#' Monte Carlo approximation of the expectation
#' of the Wright-Fisher diffusion transition density.
#'
#' Approximates the expectation appearing in the exact transition
#' density of the Wright--Fisher diffusion using Brownian bridge
#' simulation.
#'
#' @param x Target allele frequency.
#' @param x0 Initial allele frequency.
#' @param t Scaled time.
#' @param n_samples Number of Monte Carlo samples.
#' @param delta Step size for numerical integration along the path.
#' @param seed Optional integer. If provided, sets the random seed
#'   for the Monte Carlo simulation.
#'
#' @return A numeric scalar approximating the expectation.
#'
#' @details
#' The Brownian bridge is simulated on a fixed grid, and the
#' exponential functional is evaluated via numerical integration.
#'
#' @export
mc_expectation <- function(
    x, 
    x0, 
    t, 
    n_samples = 1e3, 
    delta = 0.001,
    seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Integration grid
  u <- seq(delta, 1-delta, by = delta)
  n_u <- length(u)
  
  vals <- numeric(n_samples)
  
  for (i in seq_len(n_samples)) {
    
    # Brownian bridge on the SAME grid as u
    B <- rbridge(end = 1, frequency = n_u) %>% as.numeric()
    
    # Path inside nu
    Y <- (1 - u) * f_trans(x0) +
      u * f_trans(x) +
      sqrt(t) * B
    
    # Integral of nu(Y) = 3/4 * cot^2(Y)
    int_val <- sintegral(u, 1 / (tan(Y)^2))$int
    
    # Monte Carlo exponential
    vals[i] <- exp(-3 * t / 8 * int_val)
  }
  
  mean(vals)
}

#' Exact transition density of the Wright-Fisher diffusion
#'
#' Computes the exact transition density of the Wright-Fisher diffusion
#' @param x Target allele frequency.
#' @param x0 Initial allele frequency.
#' @param t Scaled time since the initial observation.
#' @param n_samples Number of Monte Carlo samples.
#' @param delta Step size for numerical integration.
#' @param seed Optional integer. If provided, sets the random seed
#'   for the Monte Carlo simulation.
#'
#' @return Numeric value of the transition density at \code{x}.
#'
#' @export
exact_transition_density <- function(
    x, 
    x0, 
    t, 
    n_samples = 1e3, 
    delta = 0.001,
    seed = NULL) {
  
  prefactor <- (x0 * (1 - x0))^(1/4) / (x * (1 - x))^(3/4)
  kernel <- qt(f_trans(x), f_trans(x0), t)
  
  mc <- mc_expectation(
    x = x,
    x0 = x0,
    t = t,
    n_samples = n_samples,
    delta = delta,
    seed = NULL
  )
  
  prefactor * kernel * exp(-t / 4) * mc
}

#' Evaluate the exact transition density on a grid
#'
#' Computes the exact transition density on a regular grid
#' of allele frequencies, suitable for visualization.
#'
#' @param x0 Initial allele frequency.
#' @param t Scaled time.
#' @param n_samples Number of Monte Carlo samples.
#' @param delta Integration step size.
#' @param n_grid Number of grid points.
#' @param eps Boundary trimming parameter.
#' @param seed Optional integer. If provided, sets the random seed
#'   for the Monte Carlo simulation.
#'
#' @return A data frame with columns \code{x}, \code{density},
#'   \code{generation}, and \code{gen_scaled}.
#'
#' @export
exact_transition_density_grid <- function(
    x0, 
    t,
    n_samples = 500,
    delta = 0.001,
    n_grid = 500,
    eps = 1e-6,
    seed = NULL
) {
  
  x_grid <- seq(eps, 1-eps, length.out = n_grid)
  
  density_vals <- sapply(x_grid, function(xi) {
    exact_transition_density(
      x = xi,
      x0 = x0,
      t = t,
      n_samples = n_samples,
      delta = delta,
      seed = seed
    )
  })
  
  data.frame(
    x = x_grid,
    density = density_vals,
    generation = t
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Small-time approximation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Unnormalized Small-Time Approximation Density
#'
#' Computes the unnormalized small-time approximation of the
#' Wright–Fisher diffusion transition density.
#'
#' @param x Vector of evaluation points in (0, 1).
#' @param x0 Initial allele frequency in (0, 1).
#' @param N Effective population size.
#' @param t Generation.
#' @param eps Small positive constant used to avoid boundary issues.
#'
#' @return A numeric vector containing the unnormalized density values
#' evaluated at \code{x}.
#'
st_approximation_unnorm <- function(x, x0, N, t, eps = 1e-6) {
  
  f_x  <- f_trans(x)
  f_x0 <- f_trans(x0)
  
  x_safe <- pmin(pmax(x, eps), 1 - eps)
  denom  <- (x_safe * (1 - x_safe))^(3/4)
  
  ratio <- (cot(f_x) - cot(f_x0)) / (f_x - f_x0)
  
  # Replace 0/0 with limit value
  singular <- abs(f_x - f_x0) < 1e-12
  ratio[singular] <- 1 + cot(f_x0)^2
  
  result <- (x0 * (1 - x0))^(1/4) /
    denom *
    exp(-N/(2*t) * (f_x - f_x0)^2) *
    (1 - t/(2*N) * (-1/4 + (3/4) * ratio))
  
  result
}

#' Normalized Small-Time Approximation Density
#'
#' Constructs a normalized small-time approximation of the
#' Wright–Fisher diffusion transition density via numerical integration.
#'
#' @param x0 Initial allele frequency in (0, 1).
#' @param N Effective population size.
#' @param t Time parameter.
#' @param n_grid Number of grid points used for numerical integration.
#' @param eps Small positive constant used to avoid boundary issues.
#' @param x_grid Optional user-supplied grid.
#'
#' @return A list containing:
#' \itemize{
#'   \item{normalizing_constant}{Normalization constant of the density}
#'   \item{density}{Normalized density values evaluated on the grid}
#'   \item{grid}{Grid points where the density is evaluated}
#' }
#'
st_approximation <- function(x0,
                             N,
                             t,
                             n_grid = 2000,
                             eps = 1e-3,
                             x_grid = NULL)
{
  
  if (is.null(x_grid)) {
    x_grid <- seq(eps, 1 - eps, length.out = n_grid)
  }
  
  unnorm_vals <- st_approximation_unnorm(
    x = x_grid,
    x0 = x0,
    N = N,
    t = t,
    eps = eps
  )
  
  valid_idx   <- is.finite(unnorm_vals)
  x_grid      <- x_grid[valid_idx]
  unnorm_vals <- unnorm_vals[valid_idx]
  
  normalizing_constant <- Bolstad2::sintegral(x_grid, unnorm_vals)$int
  
  list(
    normalizing_constant = normalizing_constant,
    density = unnorm_vals / normalizing_constant,
    grid = x_grid
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Small-time Gaussian approximation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Small-Time Gaussian Approximation
#'
#' Computes the small-time Gaussian approximation of the
#' Wright–Fisher diffusion transition density.
#'
#' @param x Evaluation points in [0,1].
#' @param x0 Initial allele frequency.
#' @param N Effective population size.
#' @param t Generation.
#'
#' @return Numeric vector of approximate density values.
#'
st_gaussian_approximation <- function(x, x0, N, t){
  
  result <- dnorm(
    x,
    mean = x0,
    sd   = sqrt((t / N) * x0 * (1 - x0))
  )
  
  return(result)
}

#' Small-Time Gaussian Approximation on a Grid
#'
#' Evaluates the small-time Gaussian approximation on a regular grid.
#'
#' @param x0 Initial allele frequency.
#' @param N Effective population size.
#' @param t Generation.
#' @param grid_length Number of grid points.
#'
#' @return A data frame with columns \code{x} and \code{density}.
#'
st_gaussian_approximation_grid <- function(x0, N, t, grid_length = 1e4){
  
  x <- seq(0, 1, length.out = grid_length)
  
  result <- st_gaussian_approximation(x, x0, N, t)
  
  data.frame(
    x = x,
    density = result
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Beta approximation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

beta_approximation <- function(x, x0, N, t, eps = 1e-6){
  x0 <- pmin(pmax(x0, eps), 1 - eps) # avoids spike at x0 = 0 or x0 = 1
  e_xt <- x0
  var_xt <- x0*(1-x0)*(1-(1-(1/N))^t)
  
  alpha_t <- ((e_xt * (1-e_xt))/var_xt - 1)*e_xt
  beta_t <- ((e_xt * (1-e_xt))/var_xt - 1)*(1-e_xt)
  
  return(dbeta(x, shape1 = alpha_t, shape2 = beta_t))
}

beta_approximation_grid <- function(x0, N, t, grid_length = 1e4){
  x <- seq(0, 1, length.out = grid_length)
  result <- beta_approximation(x, x0, N, t)
  
  return(data.frame(x = x, density = result))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Normal approximation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

normal_wf <- function(x, x0, N, t){
  result <- dnorm(
    x,
    mean = x0,
    #sd = sqrt((1-(1-(1/N))^t) * x0 * (1 - x0))
    sd = sqrt((1-exp(-t/N)) * x0 * (1 - x0))
  )
  
  return(result)
}

normal_wf_grid <- function(x0, N, t){
  x <- seq(0, 1, length.out = N)
  result <- normal_wf(x, x0, N, t)
  
  return(data.frame(x = x, density = result))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adaptive estimator functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Beta Kernel Function
#'
#' Computes the Beta kernel used for the adaptive beta kernel density
#' estimation on the unit interval.
#'
#' @param x Vector of evaluation points in [0, 1].
#' @param t Kernel center in [0, 1].
#' @param b Bandwidth parameter.
#'
#' @return A numeric vector containing the Beta kernel evaluated at \code{x}.

beta_kernel <- function(x, t, b){
  dbeta(x,
        shape1 = t/b + 1,
        shape2 = (1 - t)/b + 1)
}

#' Optimal Beta Parameter Estimator for the Adaptive Beta Kernel Estimator
#'
#' @param data Vector of observations in the interval [0, 1].
#' @param delta Step size of the grid used to discretize [0, 1].
#' @param n_kn Number of subdivisions for the bandwidth (gamma) grid.
#' @param C1 Penalization constant (tunable, default = 1).
#' @param gamma_01 Initial value for constructing the gamma grid.
#'
#' @return A list containing:
#' \itemize{
#'   \item{x}{Grid points where the density is evaluated}
#'   \item{y}{Estimated density using the selected bandwidth}
#'   \item{k_opt}{Index of the selected bandwidth}
#'   \item{b_opt}{Selected bandwidth value}
#'   \item{gamma_grid}{Grid of candidate bandwidths}
#'   \item{full_matrix}{Matrix of density estimates for all bandwidths}
#' }
optimal_beta_estimator <- function(
    data,
    delta,
    n_kn = 100,
    C1 = 1,
    gamma_01 = 0.1,
    grid = NULL
) {
  n <- length(data)
  
  # Evaluation grid
  if(is.null(grid)){
    grid_x <- seq(0, 1, by = delta)  
  } else{
    grid_x <- grid
  }
  m_points <- length(grid_x)
  
  # 1. Construction of the gamma grid (bandwidths b)
  k_values <- 1:n_kn
  gamma_grid <- gamma_01 + (k_values / n_kn) * (2 - gamma_01)
  gamma_grid <- gamma_grid[gamma_grid <= 1]  # b must be < 1
  K <- length(gamma_grid)
  
  # 2. Pre-computation of the penalty term psi_gamma
  # L2 convergence rate: n^(-2s / (2s + 1))
  psi_gamma <- n^(-(2 * gamma_grid) / (2 * gamma_grid + 1))
  
  # 3. Construction of the estimation matrix (m_points x K)
  est_matrix <- matrix(0, nrow = m_points, ncol = K)
  
  b_grid <- n^(-1 / (2 * gamma_grid + 1))
  for (j in 1:K) {
    b <- b_grid[j]
    est_matrix[, j] <- sapply(grid_x, function(x) {
      mean(dbeta(
        data,
        shape1 = (x / b) + 1,
        shape2 = ((1 - x) / b) + 1
      ))
    })
  }
  
  # 4. Global Lepski method (Admissible set A)
  admissible <- rep(FALSE, K)
  
  for (k in 1:K) {
    consistent <- TRUE
    if (k < K) {
      for (l in (k + 1):K) {
        # Squared L2 difference
        # Integral approximation: sum(diff^2) * delta
        diff_l2 <- sum((est_matrix[, k] - est_matrix[, l])^2) * delta
        
        if (diff_l2 > 2 * C1 * psi_gamma[k]) {
          consistent <- FALSE
          break
        }
      }
    }
    admissible[k] <- consistent
  }
  
  # Select the optimal k (first admissible = smoothest)
  k_opt <- max(which(admissible))
  if (length(k_opt) == 0) k_opt <- 1
  
  return(list(
    x = grid_x,
    y = est_matrix[, k_opt],
    k_opt = k_opt,
    b_opt = b_grid[k_opt],
    gamma_grid = gamma_grid,
    full_matrix = est_matrix
  ))
}

#' Beta Kernel Density Estimator
#'
#' Computes a Beta kernel density estimator evaluated at a grid of
#' points in [0, 1].
#'
#' @param x Vector of observations in [0, 1].
#' @param n Number of grid points.
#' @param b Bandwidth parameter.
#'
#' @return A numeric vector containing the estimated density evaluated
#' at the grid points.

adaptive_beta_kernel_estimator <- function(
    x,
    delta,
    n_kn = 100,
    C1 = 1,
    gamma_01 = 0.1,
    grid = NULL
) {
  # Step 1: select optimal bandwidth using Lepski
  adapt <- optimal_beta_estimator(
    data = x,
    delta = delta,
    n_kn = n_kn,
    C1 = C1,
    gamma_01 = gamma_01,
    grid = grid
  )
  
  b_opt <- adapt$b_opt
  
  # Step 2: evaluate beta kernel estimator using b_opt
  if(is.null(grid)){
    grid <- seq(0, 1, by = delta)  
  }
  
  # compute density vector
  y <- sapply(grid, function(tt) mean(beta_kernel(x, tt, b_opt)))
  
  # return as list
  list(
    x = grid,
    y = y,
    b_opt = b_opt
  )
}
