
plot_synth_comparison <- function(res,
                                  x_limits = c(1950, 1975),
                                  y_limits = NULL,
                                  y_breaks = waiver(),
                                  y_label = "Real GDP per capita",
                                  y_formatter = scales::dollar_format(),
                                  title = "Spain and Synthetic Spain") {
  xlim_dates <- NULL
  if (!is.null(x_limits)) {
    xlim_dates <- as.Date(paste0(x_limits, "-01-01"))
  }

  p <- suppressWarnings(
    ggplot2::ggplot(res,
      type = "comparison",
      ylab = y_label,
      xlab = "Year",
      main = title,
      labels = c("Spain", "Synthetic Spain")
    )
  )

  p +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_x_date(
      expand = ggplot2::expansion(mult = c(0, 0.01)),
      date_breaks = "5 years",
      minor_breaks = NULL,
      date_labels = "%Y"
    ) +
    ggplot2::coord_cartesian(xlim = xlim_dates) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      minor_breaks = waiver(),
      labels = y_formatter,
      expand = expansion(mult = c(0, 0.02))
    )
}

plot_synth_gaps <- function(res = NULL,
                            gaps_df = NULL,
                            treatment_identifier = "Spain",
                            band_df = NULL,
                            x_limits = c(1950, 1975),
                            y_limits = NULL,
                            y_breaks = waiver(),
                            y_label = "Differential GDP per capita",
                            y_formatter = scales::dollar_format(),
                            treatment_year = NULL,
                            title = "Spain minus Synthetic Spain") {
  xlim_dates <- NULL
  if (!is.null(x_limits)) {
    xlim_dates <- as.Date(paste0(x_limits, "-01-01"))
  }

  if (!is.null(gaps_df) && nrow(gaps_df) > 0) {
    treated_df <- gaps_df |>
      dplyr::filter(Country == treatment_identifier) |>
      dplyr::mutate(year = as.Date(paste0(year, "-01-01")))

    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.5)

    if (!is.null(band_df)) {
      band_df <- band_df |>
        dplyr::mutate(year = as.Date(paste0(year, "-01-01")))

      lower_col <- dplyr::case_when(
        "q05" %in% names(band_df) ~ "q05",
        "lower" %in% names(band_df) ~ "lower",
        TRUE ~ NA_character_
      )
      upper_col <- dplyr::case_when(
        "q95" %in% names(band_df) ~ "q95",
        "upper" %in% names(band_df) ~ "upper",
        TRUE ~ NA_character_
      )
      median_col <- dplyr::case_when(
        "q50" %in% names(band_df) ~ "q50",
        "median" %in% names(band_df) ~ "median",
        "treated_gap" %in% names(band_df) ~ "treated_gap",
        TRUE ~ NA_character_
      )

      if (!is.na(lower_col) && !is.na(upper_col)) {
        band_df$lower <- band_df[[lower_col]]
        band_df$upper <- band_df[[upper_col]]
        p <- p +
          ggplot2::geom_ribbon(
            data = band_df,
            ggplot2::aes(x = year, ymin = lower, ymax = upper),
            fill = "grey70",
            alpha = 0.25,
            inherit.aes = FALSE
          )
      }
      if (!is.na(median_col)) {
        band_df$median <- band_df[[median_col]]
        p <- p +
          ggplot2::geom_line(
            data = band_df,
            ggplot2::aes(x = year, y = median),
            linetype = "dashed",
            colour = "black",
            linewidth = 0.8,
            inherit.aes = FALSE
          )
      }
    }

    p <- p +
      ggplot2::geom_line(
        data = treated_df,
        ggplot2::aes(x = year, y = gap),
        colour = "black",
        linewidth = 1.1,
        na.rm = TRUE
      ) +
      ggplot2::labs(
        y = y_label,
        x = "Year",
        title = title
      )

    if (!is.null(treatment_year)) {
      treat_date <- as.Date(paste0(treatment_year, "-01-01"))
      p <- p + ggplot2::geom_vline(xintercept = treat_date, linetype = "dashed", colour = "darkred", linewidth = 0.6)
    }

    return(
      p +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.ticks = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::scale_x_date(
          expand = ggplot2::expansion(mult = c(0, 0.01)),
          date_breaks = "5 years",
          minor_breaks = NULL,
          date_labels = "%Y"
        ) +
        ggplot2::coord_cartesian(xlim = xlim_dates) +
        ggplot2::scale_y_continuous(
          limits = y_limits,
          breaks = y_breaks,
          minor_breaks = waiver(),
          labels = y_formatter,
          expand = expansion(mult = c(0, 0.02))
        )
    )
  }

  p <- suppressWarnings(
    ggplot2::ggplot(res,
      type = "gaps",
      ylab = y_label,
      xlab = "Year",
      main = title,
      labels = c("Spain", "Synthetic Spain")
    )
  ) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_x_date(
      expand = ggplot2::expansion(mult = c(0, 0.01)),
      date_breaks = "5 years",
      minor_breaks = NULL,
      date_labels = "%Y"
    ) +
    ggplot2::coord_cartesian(xlim = xlim_dates) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      minor_breaks = waiver(),
      labels = y_formatter,
      expand = expansion(mult = c(0, 0.02))
    )

  if (!is.null(treatment_year)) {
    treat_date <- as.Date(paste0(treatment_year, "-01-01"))
    p <- p + ggplot2::geom_vline(xintercept = treat_date, linetype = "dashed", colour = "darkred", linewidth = 0.6)
  }

  p
}

plot_placebo <- function(placebo_gaps_df,
                         treatment_identifier = "Spain",
                         x_limits = c(1950, 1975),
                         y_limits = NULL,
                         y_breaks = waiver(),
                         y_label = "Differential GDP per capita",
                         y_formatter = scales::dollar_format(),
                         treatment_year = NULL,
                         title = "Placebo Test") {
  xlim_dates <- NULL
  if (!is.null(x_limits)) {
    xlim_dates <- as.Date(paste0(x_limits, "-01-01"))
  }

  plot_df <- placebo_gaps_df |>
    dplyr::filter(Country != "Average") |>
    dplyr::mutate(
      year = as.Date(paste0(year, "-01-01")),
      is_treated = Country == treatment_identifier
    )

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.5) +
    ggplot2::geom_line(
      data = dplyr::filter(plot_df, !is_treated),
      ggplot2::aes(x = year, y = gap, group = Country),
      colour = "grey70",
      linewidth = 0.5,
      alpha = 0.9
    ) +
    ggplot2::geom_line(
      data = dplyr::filter(plot_df, is_treated),
      ggplot2::aes(x = year, y = gap, group = Country),
      colour = "red",
      linewidth = 1.1
    ) +
    ggplot2::labs(
      y = y_label,
      x = "Year",
      title = title
    )

  if (!is.null(treatment_year)) {
    treat_date <- as.Date(paste0(treatment_year, "-01-01"))
    p <- p + ggplot2::geom_vline(xintercept = treat_date, linetype = "dashed", colour = "darkred", linewidth = 0.6)
  }

  p +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::scale_x_date(
      expand = ggplot2::expansion(mult = c(0, 0.01)),
      date_breaks = "5 years",
      minor_breaks = NULL,
      date_labels = "%Y"
    ) +
    ggplot2::coord_cartesian(xlim = xlim_dates) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      minor_breaks = waiver(),
      labels = y_formatter,
      expand = expansion(mult = c(0, 0.02))
    )
}

plot_post_pre_ratio <- function(ratio_df,
                                value_col = "rmspe_ratio",
                                metric_label = "RMSPE",
                                treatment_identifier = "Spain",
                                title = NULL) {
  if (is.null(title)) title <- sprintf("post-pre %s ratio", metric_label)

  ratio_df <- ratio_df |>
    dplyr::filter(is.finite(.data[[value_col]])) |>
    dplyr::mutate(is_treated = Country == treatment_identifier)

  if (nrow(ratio_df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = title)
    )
  }

  x_max <- max(ratio_df[[value_col]], na.rm = TRUE)
  x_max <- if (is.finite(x_max) && x_max > 0) x_max * 1.05 else 1

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = ratio_df,
      ggplot2::aes(x = .data[[value_col]], y = reorder(Country, .data[[value_col]])),
      fill = "black"
    ) +
    ggplot2::geom_col(
      data = dplyr::filter(ratio_df, is_treated),
      ggplot2::aes(x = .data[[value_col]], y = reorder(Country, .data[[value_col]])),
      fill = "red"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Country",
      title = title
    ) +
    ggplot2::coord_cartesian(xlim = c(0, x_max), expand = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 6),
      minor_breaks = NULL,
      expand = c(0, 0)
    )
}

plot_drop_one_gap_paths <- function(gap_paths_df,
                                    x_limits = NULL,
                                    y_limits = NULL,
                                    y_breaks = scales::pretty_breaks(n = 6),
                                    y_label = "Spain - Synthetic Spain",
                                    y_formatter = scales::comma_format(),
                                    reference_year = NULL,
                                    reference_label = NULL,
                                    title = "Gap Paths: Baseline and Leave-One-Donor-Out") {
  plot_df <- gap_paths_df |>
    dplyr::filter(is.finite(gap), !is.na(scenario))

  if (nrow(plot_df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = title)
    )
  }

  scenario_levels <- unique(as.character(plot_df$scenario))
  if ("Baseline" %in% scenario_levels) {
    scenario_levels <- c("Baseline", setdiff(scenario_levels, "Baseline"))
  }
  plot_df <- plot_df |>
    dplyr::mutate(scenario = factor(as.character(scenario), levels = scenario_levels))

  if (is.null(x_limits)) {
    x_limits <- range(plot_df$year, na.rm = TRUE)
  }
  if (is.null(y_limits)) {
    y_range <- range(plot_df$gap, na.rm = TRUE)
    if (!all(is.finite(y_range))) {
      y_range <- c(-1, 1)
    }
    y_pad <- diff(y_range) * 0.08
    if (!is.finite(y_pad) || y_pad <= 0) {
      y_pad <- max(abs(y_range), na.rm = TRUE) * 0.08
    }
    if (!is.finite(y_pad) || y_pad <= 0) {
      y_pad <- 1
    }
    y_limits <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  }

  palette_vals <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
  if (length(scenario_levels) > length(palette_vals)) {
    palette_vals <- c(palette_vals, grDevices::hcl.colors(length(scenario_levels) - length(palette_vals), "Dark 3"))
  }
  palette_vals <- palette_vals[seq_along(scenario_levels)]
  names(palette_vals) <- scenario_levels

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.5) +
    ggplot2::geom_line(
      data = dplyr::filter(plot_df, scenario != "Baseline"),
      ggplot2::aes(x = year, y = gap, colour = scenario, group = scenario),
      linewidth = 0.95,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = dplyr::filter(plot_df, scenario == "Baseline"),
      ggplot2::aes(x = year, y = gap, colour = scenario, group = scenario),
      linewidth = 1.15,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      x = "Year",
      y = y_label,
      colour = NULL,
      title = title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_colour_manual(values = palette_vals, breaks = scenario_levels) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      breaks = scales::pretty_breaks(n = 6),
      minor_breaks = NULL,
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      minor_breaks = waiver(),
      labels = y_formatter,
      expand = ggplot2::expansion(mult = c(0, 0.02))
    )

  if (!is.null(reference_year) && is.finite(reference_year)) {
    label_text <- if (!is.null(reference_label)) as.character(reference_label) else as.character(reference_year)
    label_y <- y_limits[2] - 0.05 * diff(y_limits)
    p <- p +
      ggplot2::geom_vline(xintercept = reference_year, linetype = "dashed", colour = "black", linewidth = 0.6) +
      ggplot2::annotate("text", x = reference_year + 0.1, y = label_y, label = label_text, hjust = 0, vjust = 1, size = 3.5)
  }

  p
}

