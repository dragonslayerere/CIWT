#' Plot FIW-Chi2 Results
#'
#' @param x An object of class \code{fiw_result}
#' @param alpha Numeric threshold for significance (default: 0.05)
#' @param show_pointwise Logical, whether to display the unadjusted pointwise p-values (default: FALSE)
#' @param ... Additional arguments passed to the plot method
#'
#' @return A list containing the \code{ggplot} object and a data frame of significant intervals.
#' @export
plot.ciwt_result <- function(x, alpha = 0.05, show_pointwise = FALSE, ...) {

  if (!inherits(x, "fiw_result")) {
    stop("Input must be an object of class 'fiw_result'")
  }

  # 1. Extract Contiguous Significant Intervals
  filtered_sig <- x |>
    dplyr::arrange(Time) |>
    dplyr::mutate(
      is_sig = P_Value_IWT <= alpha,
      block = cumsum(is_sig & !dplyr::lag(is_sig, default = FALSE))
    ) |>
    dplyr::filter(is_sig)

  if (nrow(filtered_sig) > 0) {
    sig_intervals <- filtered_sig |>
      dplyr::group_by(block) |>
      dplyr::summarise(
        Start = min(Time),
        End = max(Time),
        .groups = "drop"
      ) |>
      dplyr::select(-block)
  } else {
    sig_intervals <- tibble::tibble(Start = numeric(0), End = numeric(0))
  }

  # 2. Compute boundaries for right-continuous (càdlàg) intervals
  time_diffs <- diff(x$Time)
  last_step <- if (length(time_diffs) > 0) min(time_diffs) else 1

  shade_df <- x |>
    dplyr::arrange(Time) |>
    dplyr::mutate(
      xmin = Time,
      xmax = dplyr::lead(Time, default = max(Time) + last_step)
    ) |>
    dplyr::filter(P_Value_IWT <= alpha)

  # 3. Format Data for Plotting (Dynamic based on show_pointwise flag)
  plot_df <- x |> dplyr::select(Time, IWT = P_Value_IWT)

  if (show_pointwise) {
    plot_df$Pointwise <- x$P_Value_Pointwise
  }

  plot_df <- plot_df |>
    tidyr::pivot_longer(
      cols = -Time,
      names_to = "Test_Type",
      values_to = "P_Value"
    ) |>
    dplyr::mutate(Test_Type = factor(Test_Type, levels = c("IWT", "Pointwise")))

  # 4. Generate the ggplot
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = shade_df,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = 0,
        ymax = 1
      ),
      fill = "#E74C3C",
      alpha = 0.15,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = alpha,
      color = "#E74C3C",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    # Replaced geom_line with geom_step(direction = "hv")
    ggplot2::geom_step(
      data = plot_df,
      ggplot2::aes(x = Time, y = P_Value, color = Test_Type, linetype = Test_Type),
      linewidth = 1,
      direction = "hv"
    ) +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = Time, y = P_Value, color = Test_Type, shape = Test_Type),
      size = 1.5
    ) +
    ggplot2::scale_color_manual(values = c("Pointwise" = "#3498DB", "IWT" = "#2C3E50")) +
    ggplot2::scale_linetype_manual(values = c("Pointwise" = "dashed", "IWT" = "solid")) +
    ggplot2::scale_shape_manual(values = c("Pointwise" = 16, "IWT" = 16)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, by = 0.1)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "Functional Interval-Wise \U03C7\U00B2 Test",
      subtitle = paste("Shaded regions denote strict IWT significance (\U03B1 \U2264", alpha, ")"),
      y = "P-Value",
      color = "Test Type",
      linetype = "Test Type",
      shape = "Test Type"
    )

  invisible(list(
    plot = p,
    significant_intervals = sig_intervals
  ))
}
