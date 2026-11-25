library(tidyverse)

#' Simulate a Wright–Fisher allele frequency trajectory
#'
#' Simulates a single Wright–Fisher trajectory for a biallelic locus,
#' starting from an initial allele frequency and evolving through
#' discrete generations until fixation, loss, or a specified final generation.
#'
#' If \code{final_generation} is \code{NULL} and \code{generations} is provided,
#' the maximum value in \code{generations} is used as the final generation.
#'
#' @param N Integer. Effective population size.
#' @param x0 Numeric in \eqn{[0, 1]}. Initial allele frequency.
#' @param final_generation Integer. Total number of generations to simulate.
#'   Ignored if \code{generations} is provided.
#' @param generations Optional integer vector. Specific generations to return
#'   (starting at generation 0).
#'
#' @return
#' If \code{generations} is \code{NULL}, a numeric vector of length
#' \code{final_generation + 1} containing the full trajectory.
#' Otherwise, a numeric vector containing the allele frequencies at the
#' requested generations.
#'
#' @examples
#' wf_trajectory(N = 100, x0 = 0.2, final_generation = 50)
#' wf_trajectory(N = 100, x0 = 0.2, generations = c(0, 10, 20, 50))
#'
#' @export
wf_trajectory <- function(N,
                          x0,
                          final_generation = NULL,
                          generations = NULL) {

  if (!is.numeric(x0) || x0 < 0 || x0 > 1) {
    stop("x0 must be a numeric value between 0 and 1")
  }

  if (is.null(final_generation) && is.null(generations)) {
    stop("Please provide either final_generation or generations")
  }

  if (!is.null(generations)) {
    final_generation <- max(generations)
  }

  trajectory <- numeric(final_generation + 1)
  trajectory[1] <- x0

  for (i in seq_len(final_generation)) {
    trajectory[i + 1] <- rbinom(1, size = N, prob = trajectory[i]) / N

    if (trajectory[i + 1] == 0) {
      trajectory[(i + 1):(final_generation + 1)] <- 0
      break
    }

    if (trajectory[i + 1] == 1) {
      trajectory[(i + 1):(final_generation + 1)] <- 1
      break
    }
  }

  if (is.null(generations)) {
    trajectory
  } else {
    trajectory[generations + 1]
  }
}

#' Simulate multiple Wright–Fisher trajectories
#'
#' Runs \code{wf_trajectory()} multiple times to generate replicated
#' Wright–Fisher allele frequency trajectories.
#'
#' If \code{final_generation} is \code{NULL} and \code{generations} is provided,
#' the maximum value in \code{generations} is used as the final generation.
#'
#' @param n_sim Integer. Number of trajectories to simulate.
#' @param seed Optional integer. Random seed for reproducibility.
#' @param N Integer. Effective population size.
#' @param x0 Numeric in \eqn{[0, 1]}. Initial allele frequency.
#' @param final_generation Integer. Total number of generations to simulate.
#'   Ignored if \code{generations} is provided.
#' @param generations Optional integer vector. Specific generations to return
#'   (starting at generation 0).
#'
#' @return
#' A tibble with \code{n_sim} rows and one column per generation.
#' Columns are named \code{t_0}, \code{t_1}, ..., or according to
#' the requested generations.
#'
#' @examples
#' wf_trajectories(n_sim = 100, N = 100, x0 = 0.2, final_generation = 50)
#' wf_trajectories(n_sim = 100, N = 100, x0 = 0.2,
#'                 generations = c(0, 10, 20, 50))
#'
#' @export
wf_trajectories <- function(n_sim = 1e3,
                            seed = NULL,
                            N,
                            x0,
                            final_generation = NULL,
                            generations = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.null(generations)) {
    final_generation <- max(generations)
  }

  sims <- replicate(
    n_sim,
    wf_trajectory(
      N = N,
      x0 = x0,
      final_generation = final_generation,
      generations = generations
    ),
    simplify = TRUE
  )

  sims <- t(sims)

  if (is.null(generations)) {
    as_tibble(sims, .name_repair = ~ paste0("t_", 0:final_generation))
  } else {
    as_tibble(sims, .name_repair = ~ paste0("t_", generations))
  }
}

#' Plot Wright–Fisher allele frequency trajectories
#'
#' Simulates and plots multiple independent Wright–Fisher allele frequency
#' trajectories over a specified number of generations.
#'
#' Each line represents one simulated trajectory. The plot displays
#' individual trajectories.
#'
#' @param n_sim Integer. Number of trajectories to simulate. Defaults to 10.
#' @param seed Optional integer. Random seed for reproducibility.
#' @param N Integer. Effective population size for a haploid population.
#' @param x0 Numeric in \eqn{[0, 1]}. Initial allele frequency.
#' @param final_generation Integer. Total number of generations to simulate.
#' @param color Character. Color used for the trajectories.
#'   Defaults to \code{"#607B7D"}.
#'
#' @return A ggplot object displaying the simulated allele frequency
#'   trajectories over generations.
#'
#' @examples
#' plot_trajectories_wf(
#'   n_sim = 10,
#'   N = 100,
#'   x0 = 0.2,
#'   final_generation = 50
#' )
#'
#' @export
plot_trajectories_wf <- function(
    n_sim = 10,
    seed = NULL,
    N,
    x0,
    final_generation,
    color = "#607B7D"){

  df_trajectories <- wf_trajectories(
    n_sim = n_sim,
    seed = seed,
    N = N,
    x0 = x0,
    final_generation = final_generation) %>%
    rowid_to_column("replicate") %>%
    pivot_longer(
      cols = starts_with("t_"),
      names_to = "generation",
      values_to = "frequency"
    ) %>%
    mutate(generation = as.integer(sub("t_", "", generation)))

  if(n_sim < 20){
    ggplot(df_trajectories,
           aes(x = generation,
               y = frequency,
               group = replicate)) +
      geom_line(color = color) +
      theme_classic() +
      labs(x = "Generation", y = "Allele frequency") +
      theme(legend.position = "none") +
      coord_cartesian(ylim = c(0,1))
  } else {
    ggplot(df_trajectories,
           aes(x = generation,
               y = frequency,
               group = replicate)) +
      geom_line(color = color,
                alpha = 0.7) +
      theme_classic() +
      labs(x = "Generation", y = "Allele frequency") +
      theme(legend.position = "none") +
      coord_cartesian(ylim = c(0,1))
  }
}