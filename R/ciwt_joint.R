#' Run Joint Functional Interval-Wise Chi-Squared Test
#'
#' @param Y A 3D array of shape \code{[n_subj, q_states, n_times]}.
#' @param X A 3D array of shape \code{[n_subj, K_covariates, n_times]} or matrix \code{[n_subj, K_covariates]}.
#' @param Z Optional nuisance covariate matrix.
#' @param framework Character string: \code{"TCATA"} or \code{"TDS"}.
#' @param B Number of bootstrap replications (default: 1000).
#' @param use_firth Logical, whether to use Firth's penalization (default: TRUE).
#' @param scale_covariates Logical, scale nuisance covariates (default: TRUE).
#'
#' @return An object of class \code{fiw_result}.
#' @export
ciwt_joint <- function(Y, X, Z = NULL, framework = c("TCATA", "TDS"),
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

  # Format X into a 3D array to support K covariates
  if (length(dim(X)) == 2) {
    K <- dim(X)[2]
    X_arr <- array(0, dim = c(n_subj, K, n_times))
    for (t in seq_len(n_times)) X_arr[,,t] <- X
  } else {
    K <- dim(X)[2]
    X_arr <- X
  }

  # Initialize horizontally expanded arrays for the C++ backend
  # Dimension: n_subjects x (q_states * K_covariates) x n_times
  R_in <- array(0, dim = c(n_subj, q_states * K, n_times))
  E_in <- array(0, dim = c(n_subj, q_states * K, n_times))

  ctrl <- logistf::logistf.control(maxit = 500, maxstep = 5)

  for (t in seq_len(n_times)) {
    Y_t <- Y[,, t]
    X_t <- X_arr[,, t]

    z_t <- NULL
    if (!is.null(Z)) {
      if (is.matrix(Z)) z_t <- Z
      else if (is.array(Z) && length(dim(Z)) == 3) z_t <- Z[,, t]
      z_t <- ciwt_filter(z_t, scale_covariates)
    }

    for (j in seq_len(q_states)) {
      y_j <- Y_t[, j]

      # 1. Invariant state handling
      if (length(unique(y_j)) == 1) {
        for (k in seq_len(K)) {
          col_idx <- (k - 1) * q_states + j
          E_in[, col_idx, t] <- 0
          R_in[, col_idx, t] <- 0
        }
        next
      }

      if (is.null(z_t)) {
        pi_0 <- mean(y_j)
        pi_0 <- pmax(pmin(pi_0, 1 - 1e-15), 1e-15)
        e_j <- y_j - rep(pi_0, n_subj)

        for (k in seq_len(K)) {
          col_idx <- (k - 1) * q_states + j
          x_k <- X_t[, k]
          E_in[, col_idx, t] <- e_j
          R_in[, col_idx, t] <- x_k - mean(x_k)
        }
      } else {
        # 2. Fit Null Model (Strictly under H0 using Z)
        df_model <- data.frame(y = y_j, z_t)
        formula_null <- as.formula(paste("y ~", paste(colnames(z_t), collapse = " + ")))

        if (use_firth) {
          fit <- logistf::logistf(formula_null, data = df_model, pl = FALSE, control = ctrl)
          eta <- stats::model.matrix(formula_null, data = df_model) %*% fit$coefficients
          pi_vec <- 1 / (1 + exp(-eta))
        } else {
          fit <- stats::glm(formula_null, data = df_model, family = stats::binomial())
          pi_vec <- stats::fitted(fit)
        }
        pi_vec <- pmax(pmin(pi_vec, 1 - 1e-15), 1e-15)

        # 3. Extract Score Residuals
        e_j <- y_j - pi_vec
        w_j <- pi_vec * (1 - pi_vec) + 1e-6

        # 4. Iterative WLS Orthogonalization for each Target Covariate k
        for (k in seq_len(K)) {
          col_idx <- (k - 1) * q_states + j
          x_k <- X_t[, k]

          df_wls <- data.frame(x = x_k, z_t)
          formula_wls <- as.formula(paste("x ~", paste(colnames(z_t), collapse = " + ")))
          wls_fit <- stats::lm(formula_wls, data = df_wls, weights = w_j)

          # Populate the widened arrays
          E_in[, col_idx, t] <- e_j
          R_in[, col_idx, t] <- stats::residuals(wls_fit)
        }
      }
    }
  }

  # Execute C++ Backend
  results_df <- ciwt_cpp_backend(R_in = R_in, E_in = E_in, B = as.integer(B))
  results_df$Time <- seq_len(n_times)
  class(results_df) <- c("fiw_result", "data.frame")
  attr(results_df, "framework") <- framework

  return(results_df)
}
