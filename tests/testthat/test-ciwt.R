test_that("run_fiw_chi2_cpp produces correct dimensions and valid p-values", {
  n_subj <- 10
  q_states <- 2
  n_times <- 30
  B <- 200 # Smaller permutation count for fast testing

  # Generate random test cubes
  R_array <- array(rnorm(n_subj * q_states * n_times), dim = c(n_subj, q_states, n_times))
  E_array <- array(rnorm(n_subj * q_states * n_times), dim = c(n_subj, q_states, n_times))

  # Execute C++ wrapper
  res <- run_fiw_chi2_cpp(R_in = R_array, E_in = E_array, B = B)

  # 1. Structural checks
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), n_times)

  # 2. Required columns check
  expected_cols <- c("P_Value_Pointwise", "P_Value_IWT")
  for (col in expected_cols) {
    expect_true(col %in% names(res))
  }

  # 3. Statistical boundary checks (p-values must be within [0, 1])
  expect_true(all(res$P_Value_Pointwise >= 0 & res$P_Value_Pointwise <= 1))
  expect_true(all(res$P_Value_IWT >= 0 & res$P_Value_IWT <= 1))
})
test_that("plot.fiw_result enforces class validation and extracts intervals", {
  # Mock invalid object lacking 'fiw_result' class
  fake_df <- data.frame(
    Time = 1:10,
    P_Value_Pointwise = runif(10),
    P_Value_IWT = runif(10)
  )

  # Call the method directly to test the class check guard
  expect_error(plot.fiw_result(fake_df), "Input must be an object of class 'fiw_result'")

  # Test successful interval extraction structure with valid class
  class(fake_df) <- c("fiw_result", "data.frame")
  fake_df$P_Value_IWT[3:5] <- 0.01 # Inject significant region

  out <- plot(fake_df, alpha = 0.05)

  expect_type(out, "list")
  expect_true("plot" %in% names(out))
  expect_true("significant_intervals" %in% names(out))
  expect_true(is.data.frame(out$significant_intervals))
})
