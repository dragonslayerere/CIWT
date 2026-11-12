#' Convert cfda Data Frame to FIWchi2 3D Array
#'
#' @param data A data frame in \code{cfda} format with columns \code{id}, \code{time}, and \code{state}.
#' @param time_grid A numeric vector representing the common evaluation time points.
#' @param states Optional character vector of all possible states. If NULL, extracted automatically from data.
#'
#' @return A 3-dimensional array of shape \code{[n_subj, q_states, n_times]} containing binary indicators.
#' @export
cfda_to_array <- function(data, time_grid, states = NULL) {
  required_cols <- c("id", "time", "state")
  if (!all(required_cols %in% names(data))) {
    stop("cfda data frame must contain 'id', 'time', and 'state' columns.")
  }

  unique_ids <- unique(data$id)
  n_subj <- length(unique_ids)

  if (is.null(states)) {
    states <- sort(unique(as.character(data$state)))
  }
  q_states <- length(states)
  n_times <- length(time_grid)

  # Initialize 3D array [n_subj, q_states, n_times] with NAs to detect missing subjects later
  Y_array <- array(0, dim = c(n_subj, q_states, n_times),
                   dimnames = list(unique_ids, states, time_grid))

  # Process each subject's trajectory
  for (i in seq_along(unique_ids)) {
    sub_id <- unique_ids[i]
    sub_data <- data[data$id == sub_id, ]

    # SAFEGUARD: If subject exists in ID list but has no rows, skip to prevent crash
    if (nrow(sub_data) == 0) {
      Y_array[i, , ] <- NA
      next
    }

    sub_data <- sub_data[order(sub_data$time), ]

    # SAFEGUARD: all.inside = TRUE and rightmost.closed = TRUE prevent boundary NAs
    idx <- findInterval(time_grid, sub_data$time,
                        rightmost.closed = TRUE, all.inside = FALSE)
    idx[idx == 0] <- 1 # Handle pre-first-transition default state

    active_states <- sub_data$state[idx]

    for (j in seq_along(states)) {
      st <- states[j]
      Y_array[i, j, ] <- ifelse(active_states == st, 1, 0)
    }
  }

  return(Y_array)
}
