export_journal_outputs <- function(session_dir,
                                   plots_dir,
                                   tables_dir,
                                   spec_outputs,
                                   specs,
                                   treatment_identifier = "Spain") {
  figures_main_dir <- file.path(plots_dir, "main")
  figures_appendix_dir <- file.path(plots_dir, "appendix")
  tables_main_dir <- file.path(tables_dir, "main")
  tables_appendix_dir <- file.path(tables_dir, "appendix")
  for (dir_path in c(figures_main_dir, figures_appendix_dir, tables_main_dir, tables_appendix_dir)) {
    ensure_dir(dir_path)
  }

  safe_ggsave_local <- function(path, plot_obj, width = 20, height = 11, units = "cm", dpi = 600) {
    tryCatch(
      ggplot2::ggsave(path, plot_obj, width = width, height = height, units = units, dpi = dpi),
      error = function(e) message(sprintf("ggsave failed for %s: %s", path, e$message))
    )
  }

  rank_pvalue_local <- function(treated_value, placebo_values, direction = "greater") {
    placebo_values <- placebo_values[is.finite(placebo_values)]
    treated_value <- treated_value[is.finite(treated_value)]
    if (length(treated_value) == 0 || length(placebo_values) == 0) {
      return(NA_real_)
    }
    treated_value <- treated_value[1]
    if (identical(direction, "greater")) {
      (1 + sum(placebo_values >= treated_value)) / (length(placebo_values) + 1)
    } else {
      (1 + sum(placebo_values <= treated_value)) / (length(placebo_values) + 1)
    }
  }

  spec_lookup <- specs
  names(spec_lookup) <- vapply(specs, `[[`, character(1), "name")

  outcome_meta <- list(
    gdpcap = list(dep_var = "gdpcap", y_label = "Real GDP per capita", y_formatter = scales::dollar_format()),
    rknacapita = list(dep_var = "rknacapita", y_label = "Real capital stock per capita", y_formatter = scales::dollar_format()),
    hc = list(dep_var = "hc", y_label = "Human Capital", y_formatter = scales::number_format(accuracy = 0.01))
  )

  get_pool_result <- function(outcome_id, spec_name, pool_label = "all") {
    out <- spec_outputs[[outcome_id]]
    if (is.null(out)) return(NULL)
    out <- out[[spec_name]]
    if (is.null(out)) return(NULL)
    out[[pool_label]]
  }

  calc_neff <- function(donor_weights) {
    if (is.null(donor_weights) || nrow(donor_weights) == 0) return(NA_real_)
    sumsq <- sum(donor_weights$Weight^2, na.rm = TRUE)
    if (!is.finite(sumsq) || sumsq <= 0) return(NA_real_)
    1 / sumsq
  }

  comparison_series_from_result <- function(mscmt_obj, dep_var) {
    if (is.null(mscmt_obj) || is.null(mscmt_obj$combined) || is.null(mscmt_obj$combined[[dep_var]])) {
      return(NULL)
    }
    comp <- as.data.frame(mscmt_obj$combined[[dep_var]])
    if (ncol(comp) < 3) return(NULL)
    years <- suppressWarnings(as.numeric(rownames(comp)))
    if (length(years) != nrow(comp) || any(is.na(years))) {
      years <- seq_len(nrow(comp))
    }
    tibble::tibble(
      year = years,
      actual = as.numeric(comp[[1]]),
      synthetic = as.numeric(comp[[2]]),
      gap = as.numeric(comp[[3]])
    )
  }

  build_post_avg_gaps <- function(pool_obj, spec_name) {
    spec <- spec_lookup[[spec_name]]
    gaps_long <- pool_obj$placebo_outputs$gaps_long
    if (is.null(gaps_long)) return(NULL)
    gaps_long |>
      dplyr::filter(Country != "Average") |>
      dplyr::filter(year >= spec$post_window[1], year <= spec$post_window[2]) |>
      dplyr::group_by(Country) |>
      dplyr::summarise(avg_gap_post = mean(gap, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(is_treated = Country == treatment_identifier)
  }

  build_gap_band <- function(pool_obj, spec_name, period = c("pre", "post")) {
    period <- match.arg(period)
    spec <- spec_lookup[[spec_name]]
    placebo_only <- pool_obj$placebo_outputs$placebo_only
    treated_df <- pool_obj$gaps_main
    if (is.null(placebo_only) || is.null(treated_df)) {
      return(list(band = NULL, treated = NULL))
    }
    window <- if (period == "pre") spec$pre_window else spec$post_window
    band_df <- placebo_only |>
      dplyr::filter(year >= window[1], year <= window[2]) |>
      dplyr::group_by(year) |>
      dplyr::summarise(
        q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
        q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
        q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    treated_df <- treated_df |>
      dplyr::filter(Country == treatment_identifier, year >= window[1], year <= window[2])
    list(band = band_df, treated = treated_df)
  }

  plot_post_gap_histogram <- function(post_avg_gaps, title) {
    ggplot2::ggplot(post_avg_gaps, ggplot2::aes(x = avg_gap_post, fill = is_treated)) +
      ggplot2::geom_histogram(alpha = 0.7, bins = 30, colour = "white") +
      ggplot2::geom_vline(
        data = dplyr::filter(post_avg_gaps, is_treated),
        ggplot2::aes(xintercept = avg_gap_post),
        colour = "red",
        linewidth = 1
      ) +
      ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey60"), guide = "none") +
      ggplot2::labs(x = "Average gap (post period)", y = "Count", title = title) +
      ggplot2::theme_minimal()
  }

  plot_post_gap_ecdf <- function(post_avg_gaps, title) {
    ggplot2::ggplot(post_avg_gaps, ggplot2::aes(x = avg_gap_post, colour = is_treated)) +
      ggplot2::stat_ecdf(linewidth = 0.9) +
      ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
      ggplot2::labs(x = "Average gap (post period)", y = "ECDF", title = title) +
      ggplot2::theme_minimal()
  }

  plot_post_gap_rank <- function(post_avg_gaps, p_gap, title) {
    rank_df <- post_avg_gaps |>
      dplyr::arrange(avg_gap_post) |>
      dplyr::mutate(rank = dplyr::row_number())
    ggplot2::ggplot(rank_df, ggplot2::aes(x = rank, y = avg_gap_post, colour = is_treated)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
      ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey30"), guide = "none") +
      ggplot2::labs(
        x = "Rank (ascending)",
        y = "Average gap (post)",
        title = title,
        subtitle = if (!is.na(p_gap)) sprintf("Placebo p-value (greater): %.3f", p_gap) else NULL
      ) +
      ggplot2::theme_minimal()
  }

  plot_gap_band <- function(band_df, treated_df, title) {
    if (is.null(band_df) || is.null(treated_df) || nrow(band_df) == 0 || nrow(treated_df) == 0) {
      return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = title))
    }
    ggplot2::ggplot() +
      ggplot2::geom_ribbon(
        data = band_df,
        ggplot2::aes(x = year, ymin = q05, ymax = q95),
        fill = "grey70",
        alpha = 0.3
      ) +
      ggplot2::geom_line(data = band_df, ggplot2::aes(x = year, y = q50), linetype = "dashed", colour = "black") +
      ggplot2::geom_line(data = treated_df, ggplot2::aes(x = year, y = gap), colour = "red", linewidth = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
      ggplot2::labs(title = title, x = "Year", y = "Gap") +
      ggplot2::theme_minimal()
  }

  fmt_int <- function(x) {
    if (!length(x) || is.na(x) || !is.finite(x)) return(NA_character_)
    format(round(x, 0), big.mark = ",", trim = TRUE, scientific = FALSE)
  }

  fmt_num <- function(x, digits = 2) {
    if (!length(x) || is.na(x) || !is.finite(x)) return(NA_character_)
    format(round(x, digits), nsmall = digits, trim = TRUE, scientific = FALSE)
  }

  fmt_pct <- function(x, digits = 2) {
    if (!length(x) || is.na(x) || !is.finite(x)) return(NA_character_)
    sprintf("%0.*f%%", digits, x)
  }

  fmt_p <- function(x, digits = 3) {
    if (!length(x) || is.na(x) || !is.finite(x)) return(NA_character_)
    sprintf("%0.*f", digits, x)
  }

  fmt_interval <- function(lower, upper, digits = 0) {
    if (any(!is.finite(c(lower, upper)))) return(NA_character_)
    sprintf("[%s; %s]", fmt_num(lower, digits), fmt_num(upper, digits))
  }

  fmt_top_donor <- function(donor, weight, digits = 3) {
    if (!length(donor) || is.na(donor) || donor == "") return(NA_character_)
    if (!length(weight) || is.na(weight) || !is.finite(weight)) return(donor)
    sprintf("%s (%s)", donor, fmt_num(weight, digits))
  }

  write_journal_table <- function(df, path) {
    ensure_dir(dirname(path))
    readr::write_csv(df, path)
  }
  summarize_pool <- function(outcome_id, spec_name, pool_label = "all") {
    pool_obj <- get_pool_result(outcome_id, spec_name, pool_label)
    if (is.null(pool_obj)) return(NULL)

    spec <- spec_lookup[[spec_name]]
    dep_var <- outcome_meta[[outcome_id]]$dep_var
    comparison_df <- comparison_series_from_result(pool_obj$mscmt, dep_var)
    fit_df <- pool_obj$placebo_outputs$fit_metrics_df
    post_avg_gaps <- build_post_avg_gaps(pool_obj, spec_name)

    treated_avg_gap <- NA_real_
    p_gap <- NA_real_
    placebos_retained <- NA_integer_
    if (!is.null(post_avg_gaps) && nrow(post_avg_gaps) > 0) {
      treated_avg_gap <- post_avg_gaps$avg_gap_post[post_avg_gaps$Country == treatment_identifier][1]
      placebo_gaps <- post_avg_gaps$avg_gap_post[post_avg_gaps$Country != treatment_identifier]
      p_gap <- rank_pvalue_local(treated_avg_gap, placebo_gaps, direction = "greater")
      placebos_retained <- length(placebo_gaps)
    }

    p_rmspe <- NA_real_
    rmspe_ratio <- NA_real_
    pre_mspe <- NA_real_
    if (!is.null(fit_df) && nrow(fit_df) > 0) {
      treated_fit <- dplyr::filter(fit_df, Country == treatment_identifier)
      rmspe_ratio <- treated_fit$rmspe_ratio[1]
      pre_mspe <- treated_fit$pre_mspe[1]
      placebo_ratios <- fit_df$rmspe_ratio[fit_df$Country != treatment_identifier & fit_df$Country != "Average"]
      p_rmspe <- rank_pvalue_local(rmspe_ratio, placebo_ratios, direction = "greater")
      if (is.na(placebos_retained)) {
        placebos_retained <- length(placebo_ratios)
      }
    }

    avg_gap_change_pct <- NA_real_
    terminal_gap <- NA_real_
    if (!is.null(comparison_df) && nrow(comparison_df) > 0) {
      post_mask <- comparison_df$year >= spec$post_window[1] & comparison_df$year <= spec$post_window[2]
      denom <- comparison_df$synthetic[post_mask]
      numer <- comparison_df$gap[post_mask]
      valid <- is.finite(numer) & is.finite(denom) & abs(denom) > .Machine$double.eps
      if (any(valid)) {
        avg_gap_change_pct <- mean((numer[valid] / denom[valid]) * 100, na.rm = TRUE)
      }
      terminal_gap <- comparison_df$gap[comparison_df$year == spec$post_window[2]][1]
    }
    if ((!is.finite(terminal_gap) || is.na(terminal_gap)) && !is.null(pool_obj$gaps_main)) {
      terminal_gap <- pool_obj$gaps_main$gap[
        pool_obj$gaps_main$Country == treatment_identifier & pool_obj$gaps_main$year == spec$post_window[2]
      ][1]
    }

    terminal_lower <- NA_real_
    terminal_upper <- NA_real_
    if (!is.null(pool_obj$conformal) && !is.null(pool_obj$conformal$band_df)) {
      term_row <- dplyr::filter(pool_obj$conformal$band_df, year == spec$post_window[2])
      if (nrow(term_row) > 0) {
        terminal_lower <- term_row$lower[1]
        terminal_upper <- term_row$upper[1]
      }
    }

    donor_weights <- pool_obj$donor_weights
    positive_donors <- if (!is.null(donor_weights)) nrow(donor_weights) else 0L
    top_donor <- if (positive_donors > 0) donor_weights$Country[1] else NA_character_
    top_weight <- if (positive_donors > 0) donor_weights$Weight[1] else NA_real_
    n_eff <- calc_neff(donor_weights)

    included_donors <- NA_integer_
    excluded_countries <- character(0)
    if (!is.null(pool_obj$panel_info)) {
      if (!is.null(pool_obj$panel_info$pool_status)) {
        included_donors <- sum(pool_obj$panel_info$pool_status$status == "included", na.rm = TRUE)
      }
      if (!is.null(pool_obj$panel_info$excluded) && nrow(pool_obj$panel_info$excluded) > 0) {
        excluded_countries <- sort(unique(pool_obj$panel_info$excluded$country))
      }
    }

    list(
      pool_obj = pool_obj,
      comparison = comparison_df,
      post_avg_gaps = post_avg_gaps,
      pre_band = build_gap_band(pool_obj, spec_name, "pre"),
      post_band = build_gap_band(pool_obj, spec_name, "post"),
      summary = list(
        outcome = outcome_id,
        spec = spec_name,
        pool = pool_label,
        treatment_year = spec$post_window[1],
        pre_window = spec$pre_window,
        post_window = spec$post_window,
        eligible_placebos = included_donors,
        retained_placebos = placebos_retained,
        positive_donors = positive_donors,
        n_eff = n_eff,
        top_donor = top_donor,
        top_weight = top_weight,
        pre_mspe = pre_mspe,
        avg_post_gap = treated_avg_gap,
        avg_gap_change_pct = avg_gap_change_pct,
        gap_p_value = p_gap,
        rmspe_ratio = rmspe_ratio,
        rmspe_p_value = p_rmspe,
        terminal_year = spec$post_window[2],
        terminal_gap = terminal_gap,
        terminal_lower = terminal_lower,
        terminal_upper = terminal_upper,
        included_donors = included_donors,
        excluded_countries = excluded_countries
      )
    )
  }

  baseline_gdp <- summarize_pool("gdpcap", "baseline", "all")
  # Drop-one dinamico: una entrada por donante positivo del baseline (en orden
  # de peso) cuya corrida drop_<donante> exista con placebos. Asi las tablas
  # 6/A1/A2 siguen al pool que seleccione el stability gate en lugar de una
  # lista fija de donantes.
  drop_one_summaries <- list()
  if (!is.null(baseline_gdp) && !is.null(baseline_gdp$pool_obj$donor_weights)) {
    for (donor in baseline_gdp$pool_obj$donor_weights$Country) {
      drop_lbl <- sprintf("drop_%s", donor)
      drop_pool_obj <- get_pool_result("gdpcap", "baseline", drop_lbl)
      if (is.null(drop_pool_obj) || is.null(drop_pool_obj$placebo)) next
      drop_summary <- summarize_pool("gdpcap", "baseline", drop_lbl)
      if (!is.null(drop_summary)) drop_one_summaries[[donor]] <- drop_summary
    }
  }
  treat_1970 <- summarize_pool("gdpcap", "treat_1970", "all")
  treat_1980 <- summarize_pool("gdpcap", "treat_1980", "all")
  capital_summary <- summarize_pool("rknacapita", "baseline", "all")
  human_summary <- summarize_pool("hc", "baseline", "all")

  holm_gap_p <- NA_real_
  if (!is.null(baseline_gdp)) {
    pvals <- c(baseline_gdp$summary$gap_p_value, baseline_gdp$summary$rmspe_p_value)
    holm_gap_p <- stats::p.adjust(pvals, method = "holm")[1]
  }

  if (!is.null(baseline_gdp)) {
    p1 <- plot_synth_comparison(
      res = baseline_gdp$pool_obj$mscmt,
      x_limits = c(1950, 1975),
      y_label = "Real GDP per capita",
      y_formatter = outcome_meta$gdpcap$y_formatter,
      title = "Real GDP per capita in Spain and synthetic Spain"
    )
    safe_ggsave_local(file.path(figures_main_dir, "Figure_1_gdp_per_capita_spain_and_synthetic_spain.png"), p1)
    safe_ggsave_local(file.path(figures_appendix_dir, "Figure_A1_gdp_per_capita_spain_and_synthetic_spain.png"), p1)

    p2 <- plot_synth_gaps(
      gaps_df = baseline_gdp$pool_obj$gaps_main,
      treatment_identifier = treatment_identifier,
      x_limits = c(1950, 1975),
      y_label = "Real GDP per capita gap",
      y_formatter = outcome_meta$gdpcap$y_formatter,
      treatment_year = baseline_gdp$summary$treatment_year,
      title = "Real GDP per capita gap: Spain over synthetic Spain"
    )
    safe_ggsave_local(file.path(figures_main_dir, "Figure_2_gdp_per_capita_gap.png"), p2)

    p3 <- plot_placebo(
      placebo_gaps_df = baseline_gdp$pool_obj$placebo_outputs$gaps_long,
      treatment_identifier = treatment_identifier,
      x_limits = c(1950, 1975),
      y_label = "Real GDP per capita gap",
      y_formatter = outcome_meta$gdpcap$y_formatter,
      treatment_year = baseline_gdp$summary$treatment_year,
      title = "Placebo Test"
    )
    safe_ggsave_local(file.path(figures_main_dir, "Figure_3_placebo_test.png"), p3)

    p4 <- plot_post_pre_ratio(
      ratio_df = baseline_gdp$pool_obj$placebo_outputs$fit_metrics_df |>
        dplyr::select(Country, rmspe_ratio),
      value_col = "rmspe_ratio",
      metric_label = "RMSPE",
      treatment_identifier = treatment_identifier
    )
    safe_ggsave_local(file.path(figures_main_dir, "Figure_4_post_pre_rmspe_ratio.png"), p4)

    safe_ggsave_local(
      file.path(figures_main_dir, "Figure_5_post_treatment_gap_ecdf.png"),
      plot_post_gap_ecdf(baseline_gdp$post_avg_gaps, "ECDF of post-treatment gaps (placebos vs treated)")
    )
    safe_ggsave_local(
      file.path(figures_main_dir, "Figure_6_post_treatment_gap_histogram.png"),
      plot_post_gap_histogram(baseline_gdp$post_avg_gaps, "Distribution of post-treatment gaps (placebos vs treated)")
    )
    safe_ggsave_local(
      file.path(figures_main_dir, "Figure_7_ranked_post_treatment_gaps.png"),
      plot_post_gap_rank(baseline_gdp$post_avg_gaps, baseline_gdp$summary$gap_p_value, "Ranked post-treatment gaps")
    )
    safe_ggsave_local(
      file.path(figures_main_dir, "Figure_8_post_treatment_gap_bands.png"),
      plot_gap_band(baseline_gdp$post_band$band, baseline_gdp$post_band$treated, "Post-treatment gaps with placebo bands")
    )
    safe_ggsave_local(
      file.path(figures_main_dir, "Figure_9_pre_treatment_gap_bands.png"),
      plot_gap_band(baseline_gdp$pre_band$band, baseline_gdp$pre_band$treated, "Pre-treatment gaps with placebo bands")
    )
  }

  if (!is.null(baseline_gdp) && length(drop_one_summaries) > 0) {
    gap_path_df <- dplyr::bind_rows(
      baseline_gdp$pool_obj$gaps_main |>
        dplyr::filter(Country == treatment_identifier) |>
        dplyr::transmute(year = year, gap = gap, scenario = "Baseline"),
      dplyr::bind_rows(lapply(names(drop_one_summaries), function(donor) {
        drop_one_summaries[[donor]]$pool_obj$gaps_main |>
          dplyr::filter(Country == treatment_identifier) |>
          dplyr::transmute(year = year, gap = gap, scenario = paste("Drop", donor))
      }))
    )
    pA2 <- plot_drop_one_gap_paths(
      gap_paths_df = gap_path_df,
      x_limits = range(gap_path_df$year, na.rm = TRUE),
      y_label = "Spain - Synthetic Spain",
      y_formatter = outcome_meta$gdpcap$y_formatter,
      reference_year = 1959.5,
      reference_label = 1959,
      title = "Real GDP per capita gap paths under the baseline and leave-one-donor-out specifications"
    )
    safe_ggsave_local(file.path(figures_appendix_dir, "Figure_A2_gdp_per_capita_leave_one_donor_out_gap_paths.png"), pA2)
  }

  if (!is.null(capital_summary)) {
    pB1 <- plot_synth_comparison(
      res = capital_summary$pool_obj$mscmt,
      x_limits = c(1950, 1975),
      y_label = "Real capital stock per capita",
      y_formatter = outcome_meta$rknacapita$y_formatter,
      title = "Real capital stock per capita in Spain and synthetic Spain"
    )
    safe_ggsave_local(file.path(figures_appendix_dir, "Figure_B1_real_capital_stock_per_capita_spain_and_synthetic_spain.png"), pB1)

    pB2 <- plot_post_pre_ratio(
      ratio_df = capital_summary$pool_obj$placebo_outputs$fit_metrics_df |>
        dplyr::select(Country, rmspe_ratio),
      value_col = "rmspe_ratio",
      metric_label = "RMSPE",
      treatment_identifier = treatment_identifier,
      title = "Post-pre RMSPE ratio: real capital stock per capita"
    )
    safe_ggsave_local(file.path(figures_appendix_dir, "Figure_B2_real_capital_stock_per_capita_rmspe_ratio_ranking.png"), pB2)
  }

  if (!is.null(human_summary)) {
    pC1 <- plot_synth_comparison(
      res = human_summary$pool_obj$mscmt,
      x_limits = c(1950, 1975),
      y_label = "Human Capital",
      y_formatter = outcome_meta$hc$y_formatter,
      title = "Human capital in Spain and synthetic Spain"
    )
    safe_ggsave_local(file.path(figures_appendix_dir, "Figure_C1_human_capital_spain_and_synthetic_spain.png"), pC1)

    pC2 <- plot_post_pre_ratio(
      ratio_df = human_summary$pool_obj$placebo_outputs$fit_metrics_df |>
        dplyr::select(Country, rmspe_ratio),
      value_col = "rmspe_ratio",
      metric_label = "RMSPE",
      treatment_identifier = treatment_identifier,
      title = "Post-pre RMSPE ratio: human capital index"
    )
    safe_ggsave_local(file.path(figures_appendix_dir, "Figure_C2_human_capital_index_rmspe_ratio_ranking.png"), pC2)
  }
  if (!is.null(baseline_gdp)) {
    table4 <- baseline_gdp$pool_obj$donor_weights |>
      dplyr::mutate(Weight = round(Weight, 4))
    write_journal_table(table4, file.path(tables_main_dir, "Table_4_baseline_donor_weights.csv"))

    table4a <- baseline_gdp$pool_obj$predictor_table
    write_journal_table(table4a, file.path(tables_main_dir, "Table_4A_predictor_balance.csv"))

    table5 <- tibble::tibble(
      Metric = c(
        "Preferred specification",
        "Treatment year",
        "Pre-treatment window",
        "Post-treatment window",
        "Eligible placebo units after coverage filters",
        "Placebo units retained after RMSPE filter",
        "Positive donor weights",
        "Effective number of donors",
        "Largest donor weight",
        "Pre-treatment fit statistic (MSPE)",
        "Average post-treatment gap",
        "Average gap change",
        "One-sided rank p-value (average post gap)",
        "Holm-adjusted p-value",
        "Post/pre-RMSPE ratio",
        "One-sided rank p-value (post/pre RMSPE ratio)"
      ),
      Value = c(
        "baseline / all (selected specification)",
        as.character(baseline_gdp$summary$treatment_year),
        sprintf("%s-%s (%s annual observations)", baseline_gdp$summary$pre_window[1], baseline_gdp$summary$pre_window[2], baseline_gdp$summary$pre_window[2] - baseline_gdp$summary$pre_window[1] + 1),
        sprintf("%s-%s (%s annual observations)", baseline_gdp$summary$post_window[1], baseline_gdp$summary$post_window[2], baseline_gdp$summary$post_window[2] - baseline_gdp$summary$post_window[1] + 1),
        fmt_int(baseline_gdp$summary$eligible_placebos),
        fmt_int(baseline_gdp$summary$retained_placebos),
        fmt_int(baseline_gdp$summary$positive_donors),
        fmt_num(baseline_gdp$summary$n_eff, 2),
        fmt_top_donor(baseline_gdp$summary$top_donor, baseline_gdp$summary$top_weight, 4),
        fmt_num(baseline_gdp$summary$pre_mspe, 2),
        paste(fmt_int(baseline_gdp$summary$avg_post_gap), "PPP dollars"),
        fmt_pct(baseline_gdp$summary$avg_gap_change_pct, 2),
        fmt_p(baseline_gdp$summary$gap_p_value, 4),
        fmt_p(holm_gap_p, 4),
        fmt_num(baseline_gdp$summary$rmspe_ratio, 2),
        fmt_p(baseline_gdp$summary$rmspe_p_value, 4)
      )
    )
    write_journal_table(table5, file.path(tables_main_dir, "Table_5_baseline_fit_effect_size_and_placebo_inference.csv"))
  }

  if (length(drop_one_summaries) > 0) {
    drop_list <- unname(drop_one_summaries)
    table6 <- tibble::tibble(
      `Excluded donor` = names(drop_one_summaries),
      `Included donors` = vapply(drop_list, function(x) fmt_int(x$summary$included_donors), character(1)),
      `Avg. post-gap (placebo p-value)` = vapply(drop_list, function(x) sprintf("%s (p=%s)", fmt_int(x$summary$avg_post_gap), fmt_p(x$summary$gap_p_value, 3)), character(1)),
      `1975 gap` = vapply(drop_list, function(x) fmt_int(x$summary$terminal_gap), character(1)),
      `1975 conformal band` = vapply(drop_list, function(x) fmt_interval(x$summary$terminal_lower, x$summary$terminal_upper, 0), character(1)),
      `Band above zero` = vapply(drop_list, function(x) ifelse(is.finite(x$summary$terminal_lower) && x$summary$terminal_lower > 0, "Yes", "No"), character(1))
    )
    write_journal_table(table6, file.path(tables_main_dir, "Table_6_leave_one_donor_out_robustness.csv"))
  }

  if (!is.null(treat_1970) && !is.null(treat_1980)) {
    table7 <- tibble::tibble(
      `Pseudo-treatment` = c("1970", "1980"),
      `Included donors` = c(fmt_int(treat_1970$summary$included_donors), fmt_int(treat_1980$summary$included_donors)),
      `Excluded countries` = c(
        if (length(treat_1970$summary$excluded_countries) == 0) "None" else paste(treat_1970$summary$excluded_countries, collapse = ", "),
        if (length(treat_1980$summary$excluded_countries) == 0) "None" else paste(treat_1980$summary$excluded_countries, collapse = ", ")
      ),
      `Pre-MSPE` = c(fmt_int(treat_1970$summary$pre_mspe), fmt_int(treat_1980$summary$pre_mspe)),
      `Avg. post-gap (placebo p-value)` = c(
        sprintf("%s (p=%s)", fmt_int(treat_1970$summary$avg_post_gap), fmt_p(treat_1970$summary$gap_p_value, 3)),
        sprintf("%s (p=%s)", fmt_int(treat_1980$summary$avg_post_gap), fmt_p(treat_1980$summary$gap_p_value, 3))
      ),
      `Terminal result (year: gap [band])` = c(
        sprintf("%s: %s %s", treat_1970$summary$terminal_year, fmt_int(treat_1970$summary$terminal_gap), fmt_interval(treat_1970$summary$terminal_lower, treat_1970$summary$terminal_upper, 0)),
        sprintf("%s: %s %s", treat_1980$summary$terminal_year, fmt_int(treat_1980$summary$terminal_gap), fmt_interval(treat_1980$summary$terminal_lower, treat_1980$summary$terminal_upper, 0))
      )
    )
    write_journal_table(table7, file.path(tables_main_dir, "Table_7_pseudo_treatment_dates.csv"))

    appendix_d1 <- tibble::tibble(
      `Pseudo-treatment` = c("1970", "1980"),
      `Terminal year` = c(treat_1970$summary$terminal_year, treat_1980$summary$terminal_year),
      `Terminal gap` = c(fmt_int(treat_1970$summary$terminal_gap), fmt_int(treat_1980$summary$terminal_gap)),
      `Avg. post-gap` = c(fmt_int(treat_1970$summary$avg_post_gap), fmt_int(treat_1980$summary$avg_post_gap)),
      `Gap p-value` = c(fmt_p(treat_1970$summary$gap_p_value, 3), fmt_p(treat_1980$summary$gap_p_value, 3)),
      `RMSPE ratio` = c(fmt_num(treat_1970$summary$rmspe_ratio, 2), fmt_num(treat_1980$summary$rmspe_ratio, 2)),
      `RMSPE p-value` = c(fmt_p(treat_1970$summary$rmspe_p_value, 3), fmt_p(treat_1980$summary$rmspe_p_value, 3)),
      `Top donor (weight)` = c(fmt_top_donor(treat_1970$summary$top_donor, treat_1970$summary$top_weight, 3), fmt_top_donor(treat_1980$summary$top_donor, treat_1980$summary$top_weight, 3)),
      `Positive donors` = c(fmt_int(treat_1970$summary$positive_donors), fmt_int(treat_1980$summary$positive_donors))
    )
    write_journal_table(appendix_d1, file.path(tables_appendix_dir, "Table_D1_pseudo_treatment_timing_summary.csv"))
  }

  if (!is.null(baseline_gdp) && length(drop_one_summaries) > 0) {
    all_specs <- c(list(Baseline = baseline_gdp), drop_one_summaries)
    appendix_a1 <- tibble::tibble(
      Specification = c("Baseline", paste("Drop", names(drop_one_summaries))),
      `Excluded donor` = c("-", names(drop_one_summaries)),
      `Terminal gap (1975)` = vapply(all_specs, function(x) fmt_int(x$summary$terminal_gap), character(1)),
      `Avg. post-gap` = vapply(all_specs, function(x) fmt_int(x$summary$avg_post_gap), character(1)),
      `Gap p-value` = vapply(all_specs, function(x) fmt_p(x$summary$gap_p_value, 3), character(1)),
      `RMSPE ratio` = vapply(all_specs, function(x) fmt_num(x$summary$rmspe_ratio, 2), character(1)),
      `RMSPE p-value` = vapply(all_specs, function(x) fmt_p(x$summary$rmspe_p_value, 3), character(1)),
      `Positive donors` = vapply(all_specs, function(x) fmt_int(x$summary$positive_donors), character(1))
    )
    write_journal_table(appendix_a1, file.path(tables_appendix_dir, "Table_A1_gdp_per_capita_robustness_summary.csv"))

    weight_specs <- c(
      list(Baseline = baseline_gdp$pool_obj$donor_weights),
      stats::setNames(
        lapply(drop_one_summaries, function(x) x$pool_obj$donor_weights),
        paste("Drop", names(drop_one_summaries))
      )
    )
    appendix_a2 <- dplyr::bind_rows(
      lapply(
        names(weight_specs),
        function(spec_label) {
          df <- weight_specs[[spec_label]]
          if (is.null(df) || nrow(df) == 0) return(NULL)
          dplyr::transmute(df, Country = Country, specification = spec_label, Weight = round(Weight, 3))
        }
      )
    ) |>
      tidyr::pivot_wider(names_from = specification, values_from = Weight, values_fill = 0) |>
      dplyr::arrange(Country)
    write_journal_table(appendix_a2, file.path(tables_appendix_dir, "Table_A2_gdp_per_capita_donor_weights.csv"))
  }

  # Tabla A3: robustez al set de predictores (Ferman-Pinto-Possebom).
  # Mismo pool que el baseline; solo cambia el set de predictores.
  predspec_labels <- c(
    "Outcome mean + covariates (baseline)" = "all",
    "All outcome lags" = "predspec_all_lags",
    "Odd outcome lags" = "predspec_odd_lags",
    "First/mid/last outcome lags" = "predspec_first_mid_last"
  )
  predspec_summaries <- list()
  for (i in seq_along(predspec_labels)) {
    ps <- summarize_pool("gdpcap", "baseline", predspec_labels[[i]])
    if (!is.null(ps)) predspec_summaries[[names(predspec_labels)[i]]] <- ps
  }
  if (length(predspec_summaries) > 1) {
    fmt_top3_weights <- function(dw) {
      if (is.null(dw) || nrow(dw) == 0) return("-")
      n <- min(3L, nrow(dw))
      paste(sprintf("%s (%.2f)", dw$Country[seq_len(n)], dw$Weight[seq_len(n)]), collapse = "; ")
    }
    appendix_a3 <- tibble::tibble(
      `Predictor set` = names(predspec_summaries),
      `Pre-RMSPE` = vapply(predspec_summaries, function(x) fmt_num(sqrt(x$summary$pre_mspe), 1), character(1)),
      `Avg. post-gap` = vapply(predspec_summaries, function(x) fmt_int(x$summary$avg_post_gap), character(1)),
      `Terminal gap (1975)` = vapply(predspec_summaries, function(x) fmt_int(x$summary$terminal_gap), character(1)),
      `Gap p-value` = vapply(predspec_summaries, function(x) fmt_p(x$summary$gap_p_value, 3), character(1)),
      `RMSPE ratio` = vapply(predspec_summaries, function(x) fmt_num(x$summary$rmspe_ratio, 2), character(1)),
      `RMSPE p-value` = vapply(predspec_summaries, function(x) fmt_p(x$summary$rmspe_p_value, 3), character(1)),
      `Top-3 donor weights` = vapply(predspec_summaries, function(x) fmt_top3_weights(x$pool_obj$donor_weights), character(1))
    )
    write_journal_table(appendix_a3, file.path(tables_appendix_dir, "Table_A3_gdp_per_capita_predictor_set_robustness.csv"))
  }

  # Tabla A4: pools de donantes restringidos (Europa / LatAm), mismos
  # predictores que el baseline; cambia solo la composicion del pool.
  restricted_labels <- c(
    "Full pool (baseline)" = "all",
    "Europe only" = "pool_europe",
    "Europe excl. Portugal" = "pool_europe_no_portugal",
    "Latin America only" = "pool_latam"
  )
  restricted_summaries <- list()
  for (i in seq_along(restricted_labels)) {
    rs <- summarize_pool("gdpcap", "baseline", restricted_labels[[i]])
    if (!is.null(rs)) restricted_summaries[[names(restricted_labels)[i]]] <- rs
  }
  if (length(restricted_summaries) > 1) {
    fmt_top3_pool <- function(dw) {
      if (is.null(dw) || nrow(dw) == 0) return("-")
      n <- min(3L, nrow(dw))
      paste(sprintf("%s (%.2f)", dw$Country[seq_len(n)], dw$Weight[seq_len(n)]), collapse = "; ")
    }
    appendix_a4 <- tibble::tibble(
      `Donor pool` = names(restricted_summaries),
      `Donors` = vapply(restricted_summaries, function(x) fmt_int(x$summary$included_donors), character(1)),
      `Pre-RMSPE` = vapply(restricted_summaries, function(x) fmt_num(sqrt(x$summary$pre_mspe), 1), character(1)),
      `Avg. post-gap` = vapply(restricted_summaries, function(x) fmt_int(x$summary$avg_post_gap), character(1)),
      `Terminal gap (1975)` = vapply(restricted_summaries, function(x) fmt_int(x$summary$terminal_gap), character(1)),
      `Gap p-value` = vapply(restricted_summaries, function(x) fmt_p(x$summary$gap_p_value, 3), character(1)),
      `RMSPE ratio` = vapply(restricted_summaries, function(x) fmt_num(x$summary$rmspe_ratio, 2), character(1)),
      `RMSPE p-value` = vapply(restricted_summaries, function(x) fmt_p(x$summary$rmspe_p_value, 3), character(1)),
      `Top-3 donor weights` = vapply(restricted_summaries, function(x) fmt_top3_pool(x$pool_obj$donor_weights), character(1))
    )
    write_journal_table(appendix_a4, file.path(tables_appendix_dir, "Table_A4_gdp_per_capita_restricted_donor_pools.csv"))
  }

  if (!is.null(capital_summary)) {
    appendix_b1 <- tibble::tibble(
      `Terminal gap (1975)` = fmt_int(capital_summary$summary$terminal_gap),
      `Avg. post-gap` = fmt_int(capital_summary$summary$avg_post_gap),
      `Gap p-value` = fmt_p(capital_summary$summary$gap_p_value, 3),
      `RMSPE ratio` = fmt_num(capital_summary$summary$rmspe_ratio, 2),
      `RMSPE p-value` = fmt_p(capital_summary$summary$rmspe_p_value, 3),
      `Top donor (weight)` = fmt_top_donor(capital_summary$summary$top_donor, capital_summary$summary$top_weight, 3),
      `Positive donors` = fmt_int(capital_summary$summary$positive_donors)
    )
    write_journal_table(appendix_b1, file.path(tables_appendix_dir, "Table_B1_real_capital_stock_per_capita_summary_statistics.csv"))
    write_journal_table(capital_summary$pool_obj$donor_weights |> dplyr::mutate(Weight = round(Weight, 3)), file.path(tables_appendix_dir, "Table_B2_real_capital_stock_per_capita_donor_weights.csv"))
  }

  if (!is.null(human_summary)) {
    appendix_c1 <- tibble::tibble(
      `Terminal gap (1975)` = fmt_num(human_summary$summary$terminal_gap, 3),
      `Avg. post-gap` = fmt_num(human_summary$summary$avg_post_gap, 3),
      `Gap p-value` = fmt_p(human_summary$summary$gap_p_value, 3),
      `RMSPE ratio` = fmt_num(human_summary$summary$rmspe_ratio, 2),
      `RMSPE p-value` = fmt_p(human_summary$summary$rmspe_p_value, 3),
      `Top donor (weight)` = fmt_top_donor(human_summary$summary$top_donor, human_summary$summary$top_weight, 3),
      `Positive donors` = fmt_int(human_summary$summary$positive_donors)
    )
    write_journal_table(appendix_c1, file.path(tables_appendix_dir, "Table_C1_human_capital_summary_statistics.csv"))
    write_journal_table(human_summary$pool_obj$donor_weights |> dplyr::mutate(Weight = round(Weight, 3)), file.path(tables_appendix_dir, "Table_C2_human_capital_donor_weights.csv"))
  }

  for (dir_path in c(plots_dir, tables_dir)) {
    keep <- c("main", "appendix")
    existing <- list.files(dir_path, full.names = TRUE, no.. = TRUE)
    for (item in existing) {
      if (!(basename(item) %in% keep)) {
        unlink(item, recursive = TRUE, force = TRUE)
      }
    }
  }

  invisible(list(figures_main = figures_main_dir, figures_appendix = figures_appendix_dir, tables_main = tables_main_dir, tables_appendix = tables_appendix_dir))
}


