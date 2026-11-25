wf_diffusion <- function(x0,
                         final_time,
                         precision = 1e4,
                         times = NULL) {
  if (missing(x0)) {
    stop("A value for x0 must be provided.")
  }
  
  if (missing(final_time)) {
    stop("A value for final_time must be provided.")
  }
  
  delta <- final_time / precision  # time step
  
  # case 1: full trajectory
  if (is.null(times)) {
    trajectory <- numeric(precision + 1)
    trajectory[1] <- x0
    
    for (i in 2:(precision + 1)) {
      first_order_term <- sqrt(trajectory[i - 1] * (1 - trajectory[i - 1])) * 
        (sqrt(delta) * rnorm(1))
      trajectory[i] <- trajectory[i - 1] + first_order_term
      trajectory[i] <- max(0, min(1, trajectory[i]))  # keep in [0,1]
      
      if (trajectory[i] %in% c(0, 1)) {
        trajectory[(i + 1):(precision + 1)] <- trajectory[i]
        break
      }
    }
    
    names(trajectory) <- seq(0, final_time, length.out = precision + 1)
    return(trajectory)
  }
  
  # case 2: specific times
  if (!all(times <= final_time)) {
    stop("All times provided must be less than or equal to final_time.")
  }
  
  # indices of requested times
  times <- times[times > 0 & times < final_time]
  index_times <- round(times / delta)
  unique_index_times <- unique(index_times[index_times <= precision])
  
  result <- numeric(length(unique_index_times) + 2)
  result[1] <- x0
  result_index <- 2
  
  trajectory <- x0
  for (i in 2:(precision + 1)) {
    first_order_term <- sqrt(trajectory * (1 - trajectory)) * (sqrt(delta) * rnorm(1))
    trajectory <- trajectory + first_order_term
    trajectory <- max(0, min(1, trajectory))
    
    if (trajectory %in% c(0, 1)) {
      result[result_index:length(result)] <- trajectory
      break
    }
    
    if (i %in% unique_index_times) {
      result[result_index] <- trajectory
      result_index <- result_index + 1
    }
  }
  
  result[length(result)] <- trajectory
  names(result) <- c(0, times[match(unique_index_times, index_times)], final_time)
  return(result)
}

wf_diffusions <- function(n_sim = 1e3,
                          seed = NULL,
                          x0,
                          final_time,
                          precision = 1e4,
                          times = NULL) {
  if(!is.null(seed)){
    set.seed(seed)  
  }
  
  # run all simulations and store in a list
  sims <- replicate(
    n_sim,
    wf_diffusion(x0 = x0,
                 final_time = final_time,
                 precision = precision,
                 times = times),
    simplify = FALSE
  ) %>%
    t() %>%
    as_tibble()
  
  if(!is.null(generations)){
    names(results) <- paste0("t_", generations)
  } else {
    names(results) <- paste0("t_", 0:final_time)
  }
  
  return(results)
}
