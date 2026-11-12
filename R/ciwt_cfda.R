#' Run Functional Interval-Wise Chi-Squared Test Directly from cfda Data
#'
#' @param data A data frame in \code{cfda} format (`id`, `time`, `state`).
#' @param X A numeric vector or matrix of target covariates across subjects.
#' @param Z Optional nuisance covariate matrix or array.
#' @param time_grid A numeric vector of common time points to evaluate the test.
#' @param states Character vector of state names. If NULL, inferred automatically.
#' @param framework Character string: \code{"TCATA"} or \code{"TDS"}.
#' @param B Number of bootstrap replications (default: 1000).
#' @param use_firth Logical, whether to use Firth's penalization (default: TRUE).
#'
#' @return An object of class \code{fiw_result}.
#' @export
ciwt_cfda <- function(data, X, Z = NULL, time_grid, states = NULL,
                                   framework = c("TCATA", "TDS"), B = 1000, use_firth = TRUE, scale_covariates = TRUE) {
  framework <- match.arg(framework)

  Y_array <- cfda_to_array(data = data, time_grid = time_grid, states = states)

  # Pass time_grid here!
  res <- ciwt(Y = Y_array, X = X, Z = Z, time_grid = time_grid, framework = framework,
                      B = B, use_firth = use_firth, scale_covariates = scale_covariates)

  return(res)
}
