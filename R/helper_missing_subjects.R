#' Safely Remove Subjects with Missing Data (Listwise Deletion)
#' @keywords internal
remove_missing_subjects <- function(Y, X = NULL, Z = NULL, Y2 = NULL) {
  n_subj <- dim(Y)[1]
  na_mask <- rep(FALSE, n_subj)

  # 1. Check Y arrays
  na_mask <- na_mask | apply(is.na(Y), 1, any)
  if (!is.null(Y2)) na_mask <- na_mask | apply(is.na(Y2), 1, any)

  # 2. Check X
  if (!is.null(X)) {
    if (is.vector(X)) {
      na_mask <- na_mask | is.na(X)
    } else {
      na_mask <- na_mask | apply(is.na(X), 1, any)
    }
  }

  # 3. Check Z
  if (!is.null(Z)) {
    if (is.matrix(Z)) na_mask <- na_mask | apply(is.na(Z), 1, any)
    else if (is.array(Z) && length(dim(Z)) == 3) na_mask <- na_mask | apply(is.na(Z), 1, any)
  }

  # Apply filtering if NAs exist
  if (any(na_mask)) {
    n_dropped <- sum(na_mask)
    warning(sprintf("Dropped %d subject(s) due to missing (NA) values to maintain matrix alignment.", n_dropped))
    valid_idx <- !na_mask

    Y <- Y[valid_idx, , , drop = FALSE]
    if (!is.null(Y2)) Y2 <- Y2[valid_idx, , , drop = FALSE]
    if (!is.null(X)) {
      if (is.vector(X)) X <- X[valid_idx] else X <- X[valid_idx, , drop = FALSE]
    }
    if (!is.null(Z)) {
      if (is.matrix(Z)) Z <- Z[valid_idx, , drop = FALSE] else Z <- Z[valid_idx, , , drop = FALSE]
    }
  }

  return(list(Y = Y, Y2 = Y2, X = X, Z = Z, n_subj = sum(!na_mask)))
}
