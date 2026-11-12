#' Filter, Scale, and Resolve Multicollinearity for Nuisance Covariates
#'
#' @param z_t Matrix or data frame of covariates at time t.
#' @param scale_covariates Logical, whether to scale non-constant continuous variables.
#' @return A matrix of valid, linearly independent covariates, or NULL if empty.
#' @keywords internal
ciwt_filter <- function(z_t, scale_covariates = TRUE) {
  if (is.null(z_t)) return(NULL)

  z_t <- as.data.frame(z_t)
  valid_cols <- c()

  # 1. Univariate variance check and scaling
  for (col_idx in seq_len(ncol(z_t))) {
    col_vals <- z_t[[col_idx]]
    # Check if numeric AND has strictly positive variance
    if (is.numeric(col_vals) && stats::sd(col_vals, na.rm = TRUE) > 1e-8) {
      if (scale_covariates) {
        z_t[[col_idx]] <- as.vector(scale(col_vals))
      }
      valid_cols <- c(valid_cols, names(z_t)[col_idx])
    }
  }

  # Subset to valid columns
  z_t <- z_t[, valid_cols, drop = FALSE]

  # 2. Strict Multicollinearity Check via QR Decomposition
  if (ncol(z_t) > 0) {
    # Include an intercept placeholder to check the true design matrix rank
    design_matrix <- cbind(1, as.matrix(z_t))
    qr_decomp <- qr(design_matrix)

    # If the matrix is rank-deficient, silently drop the collinear duplicates
    if (qr_decomp$rank < ncol(design_matrix)) {
      keep_idx <- qr_decomp$pivot[seq_len(qr_decomp$rank)]
      keep_idx <- keep_idx[keep_idx > 1] - 1 # Remove intercept from the index
      z_t <- z_t[, keep_idx, drop = FALSE]
    }
  }

  # Return NULL if no valid columns remain
  if (ncol(z_t) == 0) return(NULL)

  return(as.matrix(z_t))
}
