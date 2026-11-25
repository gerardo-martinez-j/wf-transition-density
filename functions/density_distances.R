#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading relevant packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(Bolstad2)  # sintegral

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Declaring file paths 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

main_folder <- "/home/gerardo/Dropbox/FW-simulations"


path_density_models <- file.path(
  main_folder, "functions/density_models.R"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading functions to calculate transition density models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(path_density_models)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Hellinger distance between two densities
#'
#' Computes the Hellinger distance between two probability densities
#' evaluated on a common grid. The integral is approximated numerically
#' using \code{sintegral}.
#'
#' Negative density values (due to numerical error) are truncated to zero
#' before computing the distance.
#'
#' @param x Numeric vector. Grid of evaluation points.
#' @param f Numeric vector. First density evaluated on \code{x}.
#' @param g Numeric vector. Second density evaluated on \code{x}.
#'
#' @return A single numeric value corresponding to the Hellinger distance
#' between \code{f} and \code{g}.
#'
#' @details
#' The Hellinger distance is defined as
#' \deqn{
#' H(f,g) = \sqrt{\frac{1}{2} \int (\sqrt{f(x)} - \sqrt{g(x)})^2 dx}.
#' }
#'
#' The integral is approximated numerically over the supplied grid.
#'
#' @examples
#' x <- seq(0, 1, length.out = 1000)
#' f <- dbeta(x, 2, 5)
#' g <- dbeta(x, 3, 4)
#' hellinger_distance(x, f, g)
#'
#' @export
hellinger_distance <- function(x, f, g) {
  
  # numerical safety
  f[f < 0] <- 0
  g[g < 0] <- 0
  
  integrand <- (sqrt(f) - sqrt(g))^2
  
  integral <- sintegral(x, integrand)$int
  
  sqrt(0.5 * integral)
}

#' Compute L2 and Hellinger distances for one trajectory
#'
#' Computes L2 and Hellinger distances between a nonparametric density
#' estimator (adaptive beta kernel estimator) and several Wright–Fisher
#' transition density approximations.
#'
#' @param sample_vec Numeric vector. Sample of allele frequencies from
#' one trajectory at generation \code{t_gen}.
#' @param x0 Numeric. Initial allele frequency.
#' @param t_gen Numeric. Time in generations.
#' @param N Numeric. Effective population size.
#' @param delta Numeric. Bandwidth parameter for the adaptive beta
#' kernel estimator.
#' @param x_grid Numeric vector. Grid on which densities are evaluated.
#'
#' @return A tibble with one row containing:
#' \itemize{
#'   \item \code{l2_asymptotic} — L2 distance to the asymptotic expansion.
#'   \item \code{l2_beta} — L2 distance to the beta approximation.
#'   \item \code{l2_normal} — L2 distance to the normal approximation.
#'   \item \code{l2_gaussian} — L2 distance to the Gaussian approximation.
#'   \item \code{hell_asymptotic} — Hellinger distance to the asymptotic expansion.
#'   \item \code{hell_beta} — Hellinger distance to the beta approximation.
#'   \item \code{hell_normal} — Hellinger distance to the normal approximation.
#'   \item \code{hell_gaussian} — Hellinger distance to the Gaussian approximation.
#' }
#'
#' @details
#' The adaptive density estimate is computed using
#' \code{adaptive_beta_kernel_estimator}. Model densities are obtained
#' from:
#' \itemize{
#'   \item \code{asymptotical_expansion()}
#'   \item \code{beta_approximation()}
#'   \item \code{normal_wf()}
#'   \item \code{gaussian_approximation()}
#' }
#'
#' L2 and Hellinger distances are computed using numerical
#' integration.
#'
#' @examples
#' \dontrun{
#' compute_distances(
#'   sample_vec = rbeta(1000, 10, 5),
#'   x0 = 0.5,
#'   t_gen = 10,
#'   N = 100,
#'   delta = 0.01,
#'   x_grid = seq(0, 1, length.out = 500)
#' )
#' }
#'
#' @importFrom tibble tibble
#' @export
compute_distances <- function(sample_vec, 
                              x0, 
                              t_gen, 
                              N, 
                              delta,
                              x_grid) {
  
  # Adaptive estimator
  adaptive_est <- adaptive_beta_kernel_estimator(
    x = sample_vec,
    delta = delta,
    grid = x_grid
  )$y
  
  # Compute model densities
  model_asymp  <- st_approximation(x0 = x0, N = N, t = t_gen, x_grid = x_grid)$density
  model_beta   <- beta_approximation(x = x_grid, x0 = x0, N = N, t = t_gen)
  model_normal <- normal_wf(x = x_grid, x0 = x0, N = N, t = t_gen)
  model_gauss  <- st_gaussian_approximation(x = x_grid, x0 = x0, N = N, t = t_gen)
  
  # Continuous L2 distances via numerical integration
  l2_asymp  <- sqrt(sintegral(x_grid, (adaptive_est - model_asymp)^2)$int)
  l2_beta   <- sqrt(sintegral(x_grid, (adaptive_est - model_beta)^2)$int)
  l2_normal <- sqrt(sintegral(x_grid, (adaptive_est - model_normal)^2)$int)
  l2_gauss  <- sqrt(sintegral(x_grid, (adaptive_est - model_gauss)^2)$int)
  
  # Hellinger distances
  hell_asymp  <- hellinger_distance(x_grid, adaptive_est, model_asymp)
  hell_beta   <- hellinger_distance(x_grid, adaptive_est, model_beta)
  hell_normal <- hellinger_distance(x_grid, adaptive_est, model_normal)
  hell_gauss  <- hellinger_distance(x_grid, adaptive_est, model_gauss)
  
  tibble(
    l2_asymptotic = l2_asymp,
    l2_beta = l2_beta,
    l2_normal = l2_normal,
    l2_gaussian = l2_gauss,
    hell_asymptotic = hell_asymp,
    hell_beta = hell_beta,
    hell_normal = hell_normal,
    hell_gaussian = hell_gauss
  )
}