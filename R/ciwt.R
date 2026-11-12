#' Run Functional Interval-Wise Chi-Squared Test
#'
#' @param Y A 3D array of shape \code{[n_subj, q_states, n_times]}.
#' @param X Target covariate vector of length \code{n_subj} or matrix of shape \code{[n_subj, n_times]}.
#' @param Z Optional nuisance covariate matrix of shape \code{[n_subj, d]} or 3D array of shape \code{[n_subj, d, n_times]}.
#' @param framework Character string: \code{"TCATA"} or \code{"TDS"}.
#' @param B Number of bootstrap replications (default: 1000).
#' @param use_firth Logical, whether to use Firth's penalization (default: TRUE).
#' @param scale_covariates Logical, whether to standardize continuous nuisance covariates to zero mean and unit variance (default: TRUE).
#'
#' @return An object of class \code{fiw_result}.
#' @export
ciwt <- function(Y, X, Z = NULL, time_grid = NULL, framework = c("TCATA", "TDS"),
                         B = 1000, use_firth = TRUE, scale_covariates = TRUE) {
  framework <- match.arg(framework)

  clean_data <- remove_missing_subjects(Y = Y, X = X, Z = Z)
  Y <- clean_data$Y
  X <- clean_data$X
  Z <- clean_data$Z
  n_subj <- clean_data$n_subj

  n_subj <- dim(Y)[1]
  q_states <- dim(Y)[2]
  n_times <- dim(Y)[3]

  if (is.vector(X)) {
    X_mat <- replicate(n_times, X)
  } else {
    X_mat <- X
  }

  R_in <- array(0, dim = c(n_subj, q_states, n_times))
  E_in <- array(0, dim = c(n_subj, q_states, n_times))

  ctrl <- logistf::logistf.control(maxit = 500, maxstep = 5)

  for (t in seq_len(n_times)) {
    Y_t <- Y[,, t]
    x_t <- X_mat[, t]

    # Extract nuisance matrix Z for time t
    z_t <- NULL
    if (!is.null(Z)) {
      if (is.matrix(Z)) z_t <- Z
      else if (is.array(Z) && length(dim(Z)) == 3) z_t <- Z[,, t]
      z_t <- ciwt_filter(z_t, scale_covariates)
    }

    for (j in seq_len(q_states)) {
      y_j <- Y_t[, j]

      # Skip if outcome is invariant at time t (all 0 or all 1)
      if (length(unique(y_j)) == 1) {
        pi_val <- mean(y_j)
        E_in[, j, t] <- 0
        R_in[, j, t] <- 0
        next
      }

      if (is.null(z_t)) {
        pi_0 <- mean(y_j)
        pi_0 <- pmax(pmin(pi_0, 1 - 1e-15), 1e-15)
        pi_vec <- rep(pi_0, n_subj)

        E_in[, j, t] <- y_j - pi_vec
        R_in[, j, t] <- x_t - mean(x_t)
      } else {
        df_model <- data.frame(y = y_j, z_t)
        formula_null <- as.formula(paste("y ~", paste(colnames(z_t), collapse = " + ")))

        if (use_firth) {
          fit <- logistf::logistf(formula_null, data = df_model, pl = FALSE, control = ctrl)
          design_mat <- stats::model.matrix(formula_null, data = df_model)
          eta <- design_mat %*% fit$coefficients
          pi_vec <- 1 / (1 + exp(-eta))
        } else {
          fit <- stats::glm(formula_null, data = df_model, family = stats::binomial())
          pi_vec <- stats::fitted(fit)
        }
        pi_vec <- pmax(pmin(pi_vec, 1 - 1e-15), 1e-15)

        E_in[, j, t] <- y_j - pi_vec

        # WLS Orthogonalization
        w_j <- pi_vec * (1 - pi_vec) + 1e-6
        df_wls <- data.frame(x = x_t, z_t)
        formula_wls <- as.formula(paste("x ~", paste(colnames(z_t), collapse = " + ")))
        wls_fit <- stats::lm(formula_wls, data = df_wls, weights = w_j)
        R_in[, j, t] <- stats::residuals(wls_fit)
      }
    }
  }

  # Calculate time step durations (dt)
  if (is.null(time_grid)) {
    time_grid <- seq_len(n_times)
  }

  # State at time t is held until t+1. For the final observation, we assign a nominal weight
  # equal to the smallest observed gap to avoid dropping it entirely.
  dt <- diff(time_grid)
  dt <- c(dt, min(dt))

  # C++ Backend Execution - Now passing 'dt'
  results_df <- ciwt_cpp_backend(R_in = R_in, E_in = E_in, B = B, dt = dt)
  results_df$Time <- time_grid

  class(results_df) <- c("fiw_result", "data.frame")
  attr(results_df, "framework") <- framework

  return(results_df)
}
