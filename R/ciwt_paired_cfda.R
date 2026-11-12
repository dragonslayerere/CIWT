#' Run Functional Interval-Wise Chi-Squared Test Directly from cfda Data
#'
#' @param data1 A data frame in \code{cfda} format (`id`, `time`, `state`).
#' @param data2 A data frame in \code{cfda} format (`id`, `time`, `state`).
#' @param X A numeric vector or matrix of target covariates across subjects.
#' @param Z Optional nuisance covariate matrix or array.
#' @param time_grid A numeric vector of common time points to evaluate the test.
#' @param states Character vector of state names. If NULL, inferred automatically.
#' @param framework Character string: \code{"TCATA"} or \code{"TDS"}. (Kept for cfda_to_array compatibility if needed).
#' @param B Number of bootstrap replications (default: 1000).
#' @param use_firth Logical, ignored in paired test as OLS is used for differences.
#'
#' @return An object of class \code{fiw_result}.
#' @export
ciwt_paired_cfda <- function(data1, data2, X = NULL, Z = NULL, time_grid, states = NULL,
                                          framework = c("TCATA", "TDS"), B = 1000, use_firth = TRUE, scale_covariates = TRUE) {
  framework <- match.arg(framework)

  # 1. Convert cfda format to FIWchi2 3D array
  Y_array_1 <- cfda_to_array(data = data1, time_grid = time_grid, states = states)
  Y_array_2 <- cfda_to_array(data = data2, time_grid = time_grid, states = states)

  # 2. Execute main pipeline (Now passing time_grid)
  res <- ciwt_paired(
    Y1 = Y_array_1,
    Y2 = Y_array_2,
    X = X,
    Z = Z,
    time_grid = time_grid,
    B = B,
    scale_covariates = scale_covariates
  )
  return(res)
}
