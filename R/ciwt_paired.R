#' Run Paired Functional Interval-Wise Chi-Squared Test with Covariates
#'
#' @param Y1 A 3D array of shape \code{[n_subj, q_states, n_times]} for condition 1.
#' @param Y2 A 3D array of shape \code{[n_subj, q_states, n_times]} for condition 2.
#' @param X Optional target covariate vector of length \code{n_subj} or matrix of shape \code{[n_subj, n_times]}.
#' @param Z Optional nuisance covariate matrix or 3D array to partial out.
#' @param time_grid A numeric vector of common time points.
#' @param B Number of bootstrap replications (default: 1000).
#' @param scale_covariates Logical, whether to standardize continuous nuisance covariates.
#'
#' @return An object of class \code{fiw_result}.
#' @export
ciwt_paired <- function(Y1, Y2, X = NULL, Z = NULL, time_grid = NULL, B = 1000, scale_covariates = TRUE) {

  if (!identical(dim(Y1), dim(Y2))) {
    stop("Y1 and Y2 must have identical dimensions.")
  }

  n_subj <- dim(Y1)[1]
  q_states <- dim(Y1)[2]
  n_times <- dim(Y1)[3]

  # If X is provided, format it
  if (!is.null(X)) {
    if (is.vector(X)) X_mat <- replicate(n_times, X) else X_mat <- X
  }

  R_in <- array(1, dim = c(n_subj, q_states, n_times))
  E_in <- array(0, dim = c(n_subj, q_states, n_times))

  for (t in seq_len(n_times)) {
    # Delta Y is our response
    Delta_Y_t <- Y1[,, t] - Y2[,, t]

    # Extract and filter nuisance matrix Z
    z_t <- NULL
    if (!is.null(Z)) {
      if (is.matrix(Z)) z_t <- Z
      else if (is.array(Z) && length(dim(Z)) == 3) z_t <- Z[,, t]
      z_t <- ciwt_filter(z_t, scale_covariates)
    }

    # If X is provided at time t
    x_t <- if (!is.null(X)) X_mat[, t] else NULL

    for (j in seq_len(q_states)) {
      delta_y_j <- Delta_Y_t[, j]

      # 1. Base case: No target covariate X
      if (is.null(x_t)) {
        if (is.null(z_t)) {
          E_in[, j, t] <- delta_y_j
        } else {
          df_dy <- data.frame(dy = delta_y_j, z_t)
          form_dy <- as.formula(paste("dy ~", paste(colnames(z_t), collapse = " + ")))
          E_in[, j, t] <- stats::residuals(stats::lm(form_dy, data = df_dy))
        }
        R_in[, j, t] <- 1
      }

      # 2. Case with target covariate X
      else {
        if (is.null(z_t)) {
          E_in[, j, t] <- delta_y_j - mean(delta_y_j)
          R_in[, j, t] <- x_t - mean(x_t)
        } else {
          df_dy <- data.frame(dy = delta_y_j, z_t)
          form_dy <- as.formula(paste("dy ~", paste(colnames(z_t), collapse = " + ")))
          E_in[, j, t] <- stats::residuals(stats::lm(form_dy, data = df_dy))

          df_x <- data.frame(x = x_t, z_t)
          form_x <- as.formula(paste("x ~", paste(colnames(z_t), collapse = " + ")))
          R_in[, j, t] <- stats::residuals(stats::lm(form_x, data = df_x))
        }
      }
    }
  }

  # Calculate time step durations (dt)
  if (is.null(time_grid)) {
    time_grid <- seq_len(n_times)
  }

  dt <- diff(time_grid)
  dt <- c(dt, min(dt))

  # Pass dt to C++ backend
  results_df <- ciwt_cpp_backend(R_in = R_in, E_in = E_in, B = as.integer(B), dt = dt)

  # Assign proper time_grid instead of generic seq_len
  results_df$Time <- time_grid
  class(results_df) <- c("fiw_result", "data.frame")
  attr(results_df, "framework") <- "PAIRED"

  return(results_df)
}
