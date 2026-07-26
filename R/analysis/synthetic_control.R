run_synthetic_control <- function(session_dir) {
  plots_dir <- file.path(session_dir, OUTPUT_PLOTS)
  tables_dir <- file.path(session_dir, OUTPUT_TABLES)
  ensure_dir(plots_dir)
  ensure_dir(tables_dir)

  write_working_outputs <- FALSE
  save_table_jpg_impl <- get("save_table_jpg", mode = "function")

  # Opciones del optimizador externo (DEoptim) inyectadas en cada mscmt().
  # Amplia el presupuesto de iteraciones para asegurar convergencia estable
  # con pools grandes (ver settings.R: deoptim_itermax / deoptim_steptol).
  build_deoptim_opar <- function() {
    opar <- list()
    itmax <- if (exists("deoptim_itermax", inherits = TRUE)) suppressWarnings(as.integer(deoptim_itermax)) else NA_integer_
    stept <- if (exists("deoptim_steptol", inherits = TRUE)) suppressWarnings(as.integer(deoptim_steptol)) else NA_integer_
    if (length(itmax) == 1 && is.finite(itmax) && itmax > 0) opar$itermax <- itmax
    if (length(stept) == 1 && is.finite(stept) && stept > 0) opar$steptol <- stept
    opar
  }
  deoptim_opar <- build_deoptim_opar()

  safe_ggsave <- function(path, plot_obj, width = 20, height = 11, units = "cm", dpi = 600) {
    if (!isTRUE(write_working_outputs)) return(invisible(NULL))
    tryCatch(
      ggplot2::ggsave(path, plot_obj, width = width, height = height, units = units, dpi = dpi),
      error = function(e) message(sprintf("ggsave failed for %s: %s", path, e$message))
    )
  }

  save_table_jpg <- function(...) {
    if (!isTRUE(write_working_outputs)) return(invisible(NULL))
    save_table_jpg_impl(...)
  }

  write_working_csv <- function(...) {
    if (!isTRUE(write_working_outputs)) return(invisible(NULL))
    readr::write_csv(...)
  }

  save_working_rds <- function(...) {
    if (!isTRUE(write_working_outputs)) return(invisible(NULL))
    saveRDS(...)
  }


  write_run_metadata <- function(session_dir, metadata) {
    val <- function(x) {
      if (length(x) == 0 || is.na(x) || is.null(x)) return("")
      as.character(x)
    }
    to_block <- function(name, items) {
      if (length(items) == 0) return(character(0))
      c(sprintf("%s:", name), sprintf(" - %s", items))
    }
    lines <- c(
      sprintf("timestamp: \"%s\"", val(metadata$timestamp)),
      sprintf("session_dir: \"%s\"", val(metadata$session_dir)),
      sprintf("r_version: \"%s\"", val(metadata$r_version)),
      sprintf("platform: \"%s\"", val(metadata$platform)),
      sprintf("renv_version: \"%s\"", val(metadata$renv_version)),
      sprintf("placebo_rmspe_ratio_cutoff: \"%s\"", val(metadata$placebo_rmspe_ratio_cutoff)),
      "env_flags:",
      sprintf(" OUTCOME_ONLY: \"%s\"", val(metadata$env_flags$OUTCOME_ONLY)),
      sprintf(" SPEC_ONLY: \"%s\"", val(metadata$env_flags$SPEC_ONLY)),
      sprintf(" SKIP_DROP_ONE: \"%s\"", val(metadata$env_flags$SKIP_DROP_ONE)),
      sprintf(" placebo_rmspe_ratio_cutoff_env: \"%s\"", val(metadata$env_flags$placebo_rmspe_ratio_cutoff))
    )
    if (!is.null(metadata$stability)) {
      lines <- c(
        lines,
        "stability_gate:",
        sprintf(" enabled: \"%s\"", val(metadata$stability$enabled)),
        sprintf(" specs: \"%s\"", val(paste(metadata$stability$specs, collapse = ","))),
        sprintf(" top_weight_max: \"%s\"", val(metadata$stability$top_weight_max)),
        sprintf(" neff_min: \"%s\"", val(metadata$stability$neff_min)),
        sprintf(" min_positive_donors: \"%s\"", val(metadata$stability$min_positive_donors)),
        sprintf(" drop_top1_tau_max_pct: \"%s\"", val(metadata$stability$drop_top1_tau_max_pct)),
        sprintf(" pre_mspe_max_increase_pct: \"%s\"", val(metadata$stability$pre_mspe_max_increase_pct)),
        sprintf(" top_k_candidates: \"%s\"", val(metadata$stability$top_k_candidates))
      )
    }
    if (!is.null(metadata$outcomes) && length(metadata$outcomes) > 0) {
      lines <- c(lines, to_block("outcomes", metadata$outcomes))
    } else {
      lines <- c(lines, sprintf("outcome: \"%s\"", val(metadata$outcome)))
    }
    lines <- c(
      lines,
      to_block("specs", metadata$specs)
    )
    if (!is.null(metadata$data_files) && length(metadata$data_files) > 0) {
      lines <- c(lines, "data_files:")
      df_lines <- unlist(lapply(metadata$data_files, function(df) {
        c(
          sprintf(" - key: \"%s\"", val(df$key)),
          sprintf(" path: \"%s\"", val(df$path)),
          sprintf(" expected_md5: \"%s\"", val(df$expected_md5))
        )
      }))
      lines <- c(lines, df_lines)
    }
    ensure_dir(session_dir)
    writeLines(lines, file.path(session_dir, "run_config.yml"))
  }

  extract_weights <- function(res, controls_identifier) {
    safe_get <- function(expr) tryCatch(expr, error = function(e) NULL)
    candidates <- list(
      safe_get(res$solution.w),
      safe_get(res$solution$w),
      safe_get(res$w),
      safe_get(res$weights),
      safe_get(res$W),
      safe_get(res$Results$W)
    )
    for (w in candidates) {
      if (is.null(w)) next
      w_mat <- as.matrix(w)
      if (nrow(w_mat) == length(controls_identifier)) {
        w_vec <- w_mat[, 1]
        nms <- rownames(w_mat)
      } else if (ncol(w_mat) == length(controls_identifier)) {
        w_vec <- w_mat[1, ]
        nms <- colnames(w_mat)
      } else if (length(w_mat) == length(controls_identifier)) {
        w_vec <- as.numeric(w_mat)
        nms <- names(w_mat)
      } else {
        next
      }
      if (is.null(nms) || any(is.na(nms)) || any(nms == "")) {
        nms <- controls_identifier
      }
      names(w_vec) <- nms
      return(w_vec)
    }
    NULL
  }

  extract_gap_matrix <- function(obj, dep_var = NULL) {
    candidates <- list(
      obj$gaps,
      obj$gap,
      obj$res_gaps,
      obj$Results$gaps,
      obj$Results$gap
    )
    for (g in candidates) {
      if (is.null(g)) next
      if (is.vector(g) && !is.list(g)) {
        return(matrix(as.numeric(g), ncol = 1, dimnames = list(names(g), NULL)))
      }
      if (is.matrix(g)) return(g)
      if (is.data.frame(g)) return(as.matrix(g))
    }
    # mscmt objects often store gaps as a list per variable; pick the dependent var
    if (inherits(obj, "mscmt") && is.list(obj$gaps)) {
      dep <- dep_var
      if (is.null(dep) && !is.null(obj$dependent)) dep <- obj$dependent
      if (!is.null(dep) && !is.null(obj$gaps[[dep]])) {
        g <- as.numeric(obj$gaps[[dep]])
        return(matrix(g, ncol = 1, dimnames = list(names(obj$gaps[[dep]]), NULL)))
      }
    }
    # MSCMT placebo objects may come back as a list of mscmt objects (one per unit)
    if (is.list(obj) && any(vapply(obj, inherits, logical(1), "mscmt"))) {
      valid_idx <- which(vapply(obj, inherits, logical(1), "mscmt"))
      if (length(valid_idx) == 0) return(NULL)
      dep <- dep_var
      if (is.null(dep)) {
        dep <- obj[[valid_idx[1]]]$dependent
      }
      if (is.null(dep)) return(NULL)
      # Collect gaps for the dependent variable from each mscmt element
      gap_list <- lapply(obj[valid_idx], function(el) {
        if (!inherits(el, "mscmt")) return(NULL)
        if (is.null(el$gaps)) return(NULL)
        g <- el$gaps[[dep]]
        if (is.null(g)) return(NULL)
        as.numeric(g)
      })
      # Remove NULLs and align lengths
      valid_idx <- which(!vapply(gap_list, is.null, logical(1)))
      if (length(valid_idx) == 0) return(NULL)
      gap_list <- gap_list[valid_idx]
      lens <- vapply(gap_list, length, integer(1))
      common_len <- min(lens)
      gap_mat <- do.call(cbind, lapply(gap_list, function(g) g[seq_len(common_len)]))
      colnames(gap_mat) <- names(obj)[valid_idx]
      return(gap_mat)
    }
    NULL
  }

  resolve_years <- function(gap_mat, fallback_years) {
    yrs <- fallback_years
    if (!is.null(rownames(gap_mat))) {
      suppressWarnings({
        rn_years <- as.numeric(rownames(gap_mat))
      })
      if (!any(is.na(rn_years))) {
        yrs <- rn_years
      }
    }
    if (length(yrs) != nrow(gap_mat)) {
      yrs <- seq_len(nrow(gap_mat))
    }
    yrs
  }

  gap_matrix_to_long <- function(gap_mat, years, unit_names = NULL) {
    if (is.null(colnames(gap_mat)) || any(colnames(gap_mat) == "")) {
      if (!is.null(unit_names) && length(unit_names) == ncol(gap_mat)) {
        colnames(gap_mat) <- unit_names
      } else {
        colnames(gap_mat) <- paste0("unit_", seq_len(ncol(gap_mat)))
      }
    }
    dplyr::as_tibble(gap_mat) |>
      dplyr::mutate(year = years) |>
      tidyr::pivot_longer(cols = -year, names_to = "Country", values_to = "gap")
  }

  summarize_gap_metrics <- function(gaps_long, pre_window, post_window, treatment_identifier, spec_name, pool_label, outcome) {
    if (is.null(gaps_long)) return(NULL)
    treated <- dplyr::filter(gaps_long, Country == treatment_identifier)
    if (nrow(treated) == 0) return(NULL)
    pre_mask <- treated$year >= pre_window[1] & treated$year <= pre_window[2]
    post_mask <- treated$year >= post_window[1] & treated$year <= post_window[2]
    tibble::tibble(
      outcome = outcome,
      spec = spec_name,
      pool = pool_label,
      Country = treatment_identifier,
      pre_mspe = mean(treated$gap[pre_mask]^2, na.rm = TRUE),
      post_mspe = mean(treated$gap[post_mask]^2, na.rm = TRUE),
      avg_gap_post = mean(treated$gap[post_mask], na.rm = TRUE),
      n_pre = sum(pre_mask & !is.na(treated$gap)),
      n_post = sum(post_mask & !is.na(treated$gap))
    )
  }

  compute_weight_stability <- function(donor_weights) {
    if (is.null(donor_weights) || nrow(donor_weights) == 0) {
      return(list(
        top_donor = NA_character_,
        top_weight = NA_real_,
        n_positive = 0L,
        n_eff = NA_real_
      ))
    }
    w <- donor_weights$Weight
    sumsq <- sum(w^2, na.rm = TRUE)
    list(
      top_donor = donor_weights$Country[1],
      top_weight = donor_weights$Weight[1],
      n_positive = nrow(donor_weights),
      n_eff = if (is.finite(sumsq) && sumsq > 0) 1 / sumsq else NA_real_
    )
  }

  compute_treated_gap_stats <- function(gaps_long_main, treatment_identifier, pre_window, post_window) {
    if (is.null(gaps_long_main)) {
      return(list(pre_mspe = NA_real_, tau_post = NA_real_))
    }
    treated <- dplyr::filter(gaps_long_main, Country == treatment_identifier)
    if (nrow(treated) == 0) {
      return(list(pre_mspe = NA_real_, tau_post = NA_real_))
    }
    pre_mask <- treated$year >= pre_window[1] & treated$year <= pre_window[2]
    post_mask <- treated$year >= post_window[1] & treated$year <= post_window[2]
    pre_mspe <- mean(treated$gap[pre_mask]^2, na.rm = TRUE)
    tau_post <- mean(treated$gap[post_mask], na.rm = TRUE)
    if (!is.finite(pre_mspe)) pre_mspe <- NA_real_
    if (!is.finite(tau_post)) tau_post <- NA_real_
    list(pre_mspe = pre_mspe, tau_post = tau_post)
  }

  evaluate_stability_candidate <- function(stats, thresholds, tau_change_pct = NA_real_, pre_mspe_ref = NA_real_, require_tau = TRUE, require_pre_cap = FALSE) {
    pass_top_weight <- is.finite(stats$top_weight) && stats$top_weight <= thresholds$top_weight_max
    pass_neff <- is.finite(stats$n_eff) && stats$n_eff >= thresholds$neff_min
    pass_n_positive <- is.finite(stats$n_positive) && stats$n_positive >= thresholds$min_positive_donors
    pass_tau <- TRUE
    if (require_tau) {
      pass_tau <- is.finite(tau_change_pct) && tau_change_pct <= thresholds$drop_top1_tau_max_pct
    }
    pass_pre_cap <- TRUE
    if (require_pre_cap && is.finite(pre_mspe_ref) && pre_mspe_ref > 0) {
      max_pre <- pre_mspe_ref * (1 + thresholds$pre_mspe_max_increase_pct / 100)
      pass_pre_cap <- is.finite(stats$pre_mspe) && stats$pre_mspe <= max_pre
    }
    list(
      pass_top_weight = pass_top_weight,
      pass_neff = pass_neff,
      pass_n_positive = pass_n_positive,
      pass_tau = pass_tau,
      pass_pre_cap = pass_pre_cap,
      eligible = isTRUE(pass_top_weight && pass_neff && pass_n_positive && pass_tau && pass_pre_cap)
    )
  }

  rank_pvalue <- function(treated_value, placebo_values, direction = "greater") {
    placebo_values <- placebo_values[is.finite(placebo_values)]
    if (is.null(treated_value) || length(treated_value) == 0) return(NA_real_)
    treated_value <- treated_value[1]
    if (!is.finite(treated_value)) return(NA_real_)
    if (length(placebo_values) == 0) {
      return(NA_real_)
    }
    if (direction == "greater") {
      (1 + sum(placebo_values >= treated_value)) / (length(placebo_values) + 1)
    } else {
      (1 + sum(placebo_values <= treated_value)) / (length(placebo_values) + 1)
    }
  }

  build_fit_metrics <- function(gaps_long, pre_period, post_period, include_average = FALSE) {
    if (is.null(gaps_long) || nrow(gaps_long) == 0) {
      return(NULL)
    }

    metrics <- gaps_long
    if (!isTRUE(include_average)) {
      metrics <- dplyr::filter(metrics, Country != "Average")
    }

    metrics <- metrics |>
      dplyr::mutate(
        period = dplyr::case_when(
          year >= pre_period[1] & year <= pre_period[2] ~ "pre",
          year >= post_period[1] & year <= post_period[2] ~ "post",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(period)) |>
      dplyr::group_by(Country, period) |>
      dplyr::summarise(mspe = mean(gap^2, na.rm = TRUE), .groups = "drop") |>
      tidyr::pivot_wider(names_from = period, values_from = mspe, names_glue = "{period}_mspe")

    if (!"pre_mspe" %in% names(metrics)) metrics$pre_mspe <- NA_real_
    if (!"post_mspe" %in% names(metrics)) metrics$post_mspe <- NA_real_

    metrics |>
      dplyr::mutate(
        pre_rmspe = ifelse(is.finite(pre_mspe) & pre_mspe >= 0, sqrt(pre_mspe), NA_real_),
        post_rmspe = ifelse(is.finite(post_mspe) & post_mspe >= 0, sqrt(post_mspe), NA_real_),
        mspe_ratio = dplyr::if_else(is.finite(pre_mspe) & pre_mspe > 0, post_mspe / pre_mspe, NA_real_),
        rmspe_ratio = dplyr::if_else(is.finite(pre_rmspe) & pre_rmspe > 0, post_rmspe / pre_rmspe, NA_real_)
      )
  }

  predictor_label_map <- c(
    gdpcap = "Real GDP per capita",
    hc = "Human capital index",
    csh_c = "Consumption share",
    csh_i = "Investment share",
    csh_g = "Government share",
    csh_x = "Exports share",
    csh_m = "Imports share",
    rknacapita = "Capital stock per capita",
    pl_gdpo = "Price level of output"
  )

  format_predictor_label <- function(predictor_name) {
    if (is.null(predictor_name) || length(predictor_name) == 0 || is.na(predictor_name)) {
      return(NA_character_)
    }
    base_name <- sub("\\.mean\\..*$", "", predictor_name)
    label <- unname(predictor_label_map[base_name])
    if (length(label) == 0 || is.na(label)) {
      label <- base_name
    }
    period_match <- regexec("^.+\\.mean\\.(\\d{4})\\.(\\d{4})$", predictor_name)
    period_parts <- regmatches(predictor_name, period_match)[[1]]
    if (length(period_parts) == 3) {
      return(sprintf("%s (%s-%s mean)", label, period_parts[2], period_parts[3]))
    }
    label
  }

  build_predictor_balance_table <- function(raw_predictor_table) {
    if (is.null(raw_predictor_table) || nrow(raw_predictor_table) == 0) {
      return(NULL)
    }

    predictor_table <- as.data.frame(raw_predictor_table)
    predictor_table <- tibble::rownames_to_column(predictor_table, var = "Predictor_raw")
    colnames(predictor_table)[2:4] <- c("Spain", "Synthetic Spain", "Sample mean")

    predictor_table |>
      dplyr::transmute(
        Predictor = vapply(Predictor_raw, format_predictor_label, character(1)),
        Spain = Spain,
        `Synthetic Spain` = `Synthetic Spain`,
        `Sample mean` = `Sample mean`
      ) |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))
  }

  write_fit_metric_outputs <- function(fit_metrics_df, tables_dir, treatment_identifier) {
    if (is.null(fit_metrics_df) || nrow(fit_metrics_df) == 0) {
      return(invisible(NULL))
    }

    placebo_fit_metrics_df <- fit_metrics_df |>
      dplyr::filter(Country != treatment_identifier, Country != "Average")

    write_working_csv(fit_metrics_df, file.path(tables_dir, "unit_fit_metrics.csv"))
    write_working_csv(
      fit_metrics_df |>
        dplyr::select(Country, pre_rmspe, post_rmspe, rmspe_ratio),
      file.path(tables_dir, "unit_rmspe.csv")
    )
    write_working_csv(
      placebo_fit_metrics_df |>
        dplyr::select(Country, pre_mspe, post_mspe, mspe_ratio),
      file.path(tables_dir, "placebo_mspe.csv")
    )
    write_working_csv(
      placebo_fit_metrics_df |>
        dplyr::select(Country, pre_rmspe, post_rmspe, rmspe_ratio),
      file.path(tables_dir, "placebo_rmspe.csv")
    )
  }

  build_placebo_outputs <- function(resplacebo, years, unit_names, treatment_identifier, pre_period, post_period, tables_dir, dep_var = NULL) {
    gap_mat <- extract_gap_matrix(resplacebo, dep_var = dep_var)
    if (is.null(gap_mat)) {
      return(list(gaps_long = NULL, band_df = NULL, fit_metrics_df = NULL, placebo_only = NULL))
    }
    years <- resolve_years(gap_mat, years)
    colnames(gap_mat) <- if (is.null(colnames(gap_mat))) unit_names else colnames(gap_mat)
    gaps_long <- gap_matrix_to_long(gap_mat, years, unit_names)
    write_working_csv(gaps_long, file.path(tables_dir, "placebo_gaps_long.csv"))

    placebo_only <- dplyr::filter(gaps_long, Country != treatment_identifier, Country != "Average")

    band_df <- placebo_only |>
      dplyr::group_by(year) |>
      dplyr::summarise(
        q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
        q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
        q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    write_working_csv(band_df, file.path(tables_dir, "placebo_gap_quantiles.csv"))

    fit_metrics_df <- build_fit_metrics(
      gaps_long = gaps_long,
      pre_period = pre_period,
      post_period = post_period,
      include_average = FALSE
    )
    write_fit_metric_outputs(fit_metrics_df, tables_dir, treatment_identifier)

    list(gaps_long = gaps_long, band_df = band_df, fit_metrics_df = fit_metrics_df, placebo_only = placebo_only)
  }

  conformal_bands_from_placebos <- function(treated_gaps, placebo_long, alpha = 0.1, treatment_identifier = "Spain") {
    if (is.null(treated_gaps) || is.null(placebo_long)) {
      return(list(band_df = NULL, pvals_df = NULL))
    }
    treated <- dplyr::filter(treated_gaps, Country == treatment_identifier)
    placebo_only <- placebo_long |>
      dplyr::filter(Country != treatment_identifier, Country != "Average")
    if (nrow(treated) == 0 || nrow(placebo_only) == 0) {
      return(list(band_df = NULL, pvals_df = NULL))
    }
    calc <- lapply(seq_len(nrow(treated)), function(i) {
      g_t <- treated$gap[i]
      y_t <- treated$year[i]
      pool <- placebo_only$gap[placebo_only$year == y_t]
      pool <- pool[!is.na(pool)]
      if (length(pool) == 0 || is.na(g_t)) {
        return(NULL)
      }
      q <- stats::quantile(abs(pool), probs = 1 - alpha, na.rm = TRUE, names = FALSE)
      p_val <- (1 + sum(abs(pool) >= abs(g_t)))/(length(pool) + 1)
      dplyr::tibble(
        year = y_t,
        treated_gap = g_t,
        conformal_q = q,
        lower = g_t - q,
        upper = g_t + q,
        placebo_n = length(pool),
        p_value = p_val
      )
    })
    out <- dplyr::bind_rows(calc)
    if (nrow(out) == 0) {
      return(list(band_df = NULL, pvals_df = NULL))
    }
    band_df <- out |>
      dplyr::select(year, lower, upper, median = treated_gap, conformal_q, placebo_n, p_value)
    pvals_df <- out |>
      dplyr::select(year, p_value, placebo_n)
    list(band_df = band_df, pvals_df = pvals_df)
  }

  prepare_panel_for_spec <- function(panel_full, pre_window, post_window, key_vars, outcome_var, spec_name, pool_label, tables_dir, treatment_identifier, donor_filter = NULL, write_outputs = TRUE) {
    df <- panel_full
    if (!is.null(donor_filter)) {
      df <- dplyr::filter(
        df,
        country %in% donor_filter | country %in% treatment_identifier | country == "Average"
      )
    }
    base_pool <- unique(df$country[df$country != treatment_identifier & df$country != "Average"])
    df <- dplyr::filter(df, year >= pre_window[1], year <= post_window[2])
    pre_len <- pre_window[2] - pre_window[1] + 1
    post_len <- post_window[2] - post_window[1] + 1

    coverage_pre <- df |>
      dplyr::filter(country != "Average", year >= pre_window[1], year <= pre_window[2]) |>
      dplyr::group_by(country) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(key_vars), \(x) sum(!is.na(x)), .names = "n_{.col}"),
        .groups = "drop"
      )

    coverage_post <- df |>
      dplyr::filter(country != "Average", year >= post_window[1], year <= post_window[2]) |>
      dplyr::group_by(country) |>
      dplyr::summarise(
        n_outcome = sum(!is.na(.data[[outcome_var]])),
        .groups = "drop"
      )

    if (write_outputs) {
      write_working_csv(coverage_pre, file.path(tables_dir, "pwt_coverage_pre.csv"))
      save_table_jpg(coverage_pre, file.path(tables_dir, "pwt_coverage_pre.jpg"), title = sprintf("PWT coverage pre (%s-%s)", pre_window[1], pre_window[2]))
      write_working_csv(coverage_post, file.path(tables_dir, "pwt_coverage_post_outcome.csv"))
      save_table_jpg(coverage_post, file.path(tables_dir, "pwt_coverage_post_outcome.jpg"), title = sprintf("PWT coverage post outcome (%s-%s)", post_window[1], post_window[2]))
    }

    excluded_pre <- coverage_pre |>
      tidyr::pivot_longer(cols = dplyr::starts_with("n_"), names_to = "var", values_to = "non_na") |>
      dplyr::mutate(var = gsub("^n_", "", var)) |>
      dplyr::filter(non_na < pre_len) |>
      dplyr::mutate(reason = paste0("Insufficient pre data for ", var, " (", non_na, "/", pre_len, ")")) |>
      dplyr::select(country, reason) |>
      dplyr::distinct()

    excluded_post <- coverage_post |>
      dplyr::filter(n_outcome < post_len) |>
      dplyr::mutate(
        var = outcome_var,
        reason = paste0("Insufficient post data for ", outcome_var, " (", n_outcome, "/", post_len, ")")
      ) |>
      dplyr::select(country, reason) |>
      dplyr::distinct()

    excluded <- dplyr::bind_rows(excluded_pre, excluded_post) |>
      dplyr::group_by(country) |>
      dplyr::summarise(reason = paste(unique(reason), collapse = "; "), .groups = "drop")

    if (nrow(excluded) > 0) {
      if (write_outputs) {
        write_working_csv(excluded, file.path(tables_dir, "excluded_countries.csv"))
        save_table_jpg(excluded, file.path(tables_dir, "excluded_countries.jpg"), title = "Excluded countries (coverage)")
        if (nrow(excluded_pre) > 0) {
          write_working_csv(excluded_pre, file.path(tables_dir, "excluded_countries_pre.csv"))
          save_table_jpg(excluded_pre, file.path(tables_dir, "excluded_countries_pre.jpg"), title = "Excluded (pre, predictors)")
        }
        if (nrow(excluded_post) > 0) {
          write_working_csv(excluded_post, file.path(tables_dir, "excluded_countries_post_outcome.csv"))
          save_table_jpg(excluded_post, file.path(tables_dir, "excluded_countries_post_outcome.jpg"), title = "Excluded (post, outcome)")
        }
      }
      df <- dplyr::filter(df, !country %in% excluded$country)
    }

    panel_no_avg <- df |>
      dplyr::filter(country != "Average")
    panel_donors_only <- panel_no_avg |>
      dplyr::filter(country != treatment_identifier)
    if (nrow(panel_donors_only) == 0) stop(sprintf("No donor countries remain after coverage filters for spec %s / pool %s.", spec_name, pool_label))

    # Pool estable de donantes: se fija una sola vez tras los filtros de cobertura
    remaining_pool <- sort(unique(panel_donors_only$country))

    region_avg <- max(panel_donors_only$region, na.rm = TRUE) + 1
    if (!is.finite(region_avg)) region_avg <- 1

    panel_included <- panel_no_avg |>
      dplyr::filter(country == treatment_identifier | country %in% remaining_pool)

    numeric_cols <- setdiff(names(dplyr::select(panel_donors_only, where(is.numeric))), c("region", "year", "pop"))
    avg_clean <- panel_donors_only |>
      dplyr::filter(country %in% remaining_pool) |>
      dplyr::group_by(year) |>
      dplyr::summarise(
        pop = sum(pop, na.rm = TRUE),
        dplyr::across(
          dplyr::all_of(numeric_cols),
          ~ {
            wts <- dplyr::cur_data_all()$pop
            valid <- !is.na(.x) & !is.na(wts)
            if (!any(valid)) return(NA_real_)
            stats::weighted.mean(.x[valid], w = wts[valid])
          }
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(country = "Average", isocode = "AVG", region = region_avg)

    panel_df <- dplyr::bind_rows(avg_clean, panel_included) |>
      dplyr::mutate(country = forcats::fct_relevel(country, "Average")) |>
      dplyr::arrange(country, year)

    panel_df$region <- as.integer(panel_df$region)
    panel_df <- as.data.frame(panel_df)
    panel_df$country <- as.character(panel_df$country)
    panel_df$region <- as.numeric(panel_df$region)

    pool_status <- tibble::tibble(
      country = sort(unique(c(base_pool, remaining_pool, excluded$country))),
      status = dplyr::if_else(country %in% remaining_pool, "included", "excluded"),
      phase = dplyr::case_when(
        country %in% excluded_pre$country ~ "pre_predictors",
        country %in% excluded_post$country ~ "post_outcome",
        TRUE ~ NA_character_
      )
    ) |>
      dplyr::left_join(
        dplyr::bind_rows(
          dplyr::mutate(excluded_pre, phase = "pre_predictors"),
          dplyr::mutate(excluded_post, phase = "post_outcome")
        ),
        by = c("country", "phase")
      )

    if (write_outputs) {
      final_pool <- tibble::tibble(country = remaining_pool)
      write_working_csv(pool_status, file.path(tables_dir, "pool_status.csv"))
      save_table_jpg(pool_status, file.path(tables_dir, "pool_status.jpg"), title = "Pool status (included/excluded)")
      write_working_csv(final_pool, file.path(tables_dir, "final_pool.csv"))
      save_table_jpg(final_pool, file.path(tables_dir, "final_pool.jpg"), title = "Final donor pool")
    }

    pool_status_rows[[length(pool_status_rows) + 1]] <<- dplyr::mutate(
      pool_status,
      outcome = outcome$id,
      spec = spec$name,
      pool = pool_label
    ) |>
      dplyr::relocate(outcome, spec, pool)

    pool_summary_rows[[length(pool_summary_rows) + 1]] <<- tibble::tibble(
      outcome = outcome$id,
      spec = spec$name,
      pool = pool_label,
      donors_included = sum(pool_status$status == "included", na.rm = TRUE),
      donors_excluded = sum(pool_status$status != "included", na.rm = TRUE),
      exclusion_reasons = paste(sort(unique(stats::na.omit(pool_status$reason))), collapse = "; ")
    )

    list(
      panel = panel_df,
      coverage_pre = coverage_pre,
      coverage_post = coverage_post,
      excluded = excluded,
      pool_status = pool_status,
      key_vars = key_vars
    )
  }

  specs <- list(
    list(name = "baseline", label = "Baseline (1950-1959 pre, 1960-1975 post)", pre_window = c(1950, 1959), post_window = c(1960, 1975)),
    list(name = "treat_1970", label = "Treatment 1970 (1960-1969 pre, 1970-1985 post)", pre_window = c(1960, 1969), post_window = c(1970, 1985)),
    list(name = "treat_1980", label = "Treatment 1980 (1970-1979 pre, 1980-1995 post)", pre_window = c(1970, 1979), post_window = c(1980, 1995))
  )

  spec_env <- Sys.getenv("SPEC_ONLY", "")
  if (nzchar(spec_env)) {
    wanted_specs <- trimws(strsplit(spec_env, ",")[[1]])
    specs <- Filter(function(s) s$name %in% wanted_specs, specs)
  }

  outcomes <- list(
    list(
      id = "gdpcap",
      dep_var = "gdpcap",
      y_label = "Real GDP per capita",
      y_formatter = scales::dollar_format(),
      specs = c("baseline", "treat_1970", "treat_1980"),
      run_drop_one = TRUE,
      use_stability_gate = TRUE
    ),
    list(
      id = "rknacapita",
      dep_var = "rknacapita",
      y_label = "Real capital stock per capita",
      y_formatter = scales::dollar_format(),
      specs = c("baseline"),
      run_drop_one = FALSE,
      use_stability_gate = FALSE
    ),
    list(
      id = "hc",
      dep_var = "hc",
      y_label = "Human capital index",
      y_formatter = scales::number_format(accuracy = 0.01),
      specs = c("baseline"),
      run_drop_one = FALSE,
      use_stability_gate = FALSE
    )
  )

  # Filtro de configuracion: outcomes activos segun settings.R (enabled_outcomes).
  # OUTCOME_ONLY (variable de entorno) se aplica despues como override ad hoc.
  if (exists("enabled_outcomes", inherits = TRUE) && length(enabled_outcomes) > 0) {
    outcomes <- Filter(function(o) o$id %in% as.character(enabled_outcomes), outcomes)
  }

  outcome_env <- Sys.getenv("OUTCOME_ONLY", "")
  if (nzchar(outcome_env)) {
    wanted_outcomes <- trimws(strsplit(outcome_env, ",")[[1]])
    outcomes <- Filter(function(o) o$id %in% wanted_outcomes, outcomes)
  }
  if (length(outcomes) == 0) {
    stop("No outcomes selected after applying enabled_outcomes / OUTCOME_ONLY.")
  }

  env_flags <- Sys.getenv(c("OUTCOME_ONLY", "SPEC_ONLY", "SKIP_DROP_ONE", "placebo_mspe_ratio_cutoff"))
  names(env_flags)[names(env_flags) == "placebo_mspe_ratio_cutoff"] <- "placebo_rmspe_ratio_cutoff"
  stability_gate <- list(
    enabled = if (exists("stability_gate_enabled", inherits = TRUE)) isTRUE(stability_gate_enabled) else FALSE,
    specs = if (exists("stability_gate_specs", inherits = TRUE) && length(stability_gate_specs) > 0) as.character(stability_gate_specs) else "baseline",
    top_weight_max = if (exists("stability_top_weight_max", inherits = TRUE)) as.numeric(stability_top_weight_max) else 0.40,
    neff_min = if (exists("stability_neff_min", inherits = TRUE)) as.numeric(stability_neff_min) else 3,
    min_positive_donors = if (exists("stability_min_positive_donors", inherits = TRUE)) as.integer(stability_min_positive_donors) else 4L,
    drop_top1_tau_max_pct = if (exists("stability_drop_top1_tau_max_pct", inherits = TRUE)) as.numeric(stability_drop_top1_tau_max_pct) else 25,
    pre_mspe_max_increase_pct = if (exists("stability_pre_mspe_max_increase_pct", inherits = TRUE)) as.numeric(stability_pre_mspe_max_increase_pct) else 30,
    top_k_candidates = if (exists("stability_top_k_candidates", inherits = TRUE)) as.integer(stability_top_k_candidates) else 3L
  )
  checksum_info <- tryCatch(read_checksums(), error = function(e) e)
  checksum_lookup <- function(key) {
    if (inherits(checksum_info, "error") || is.null(checksum_info[[key]])) return(NA_character_)
    checksum_info[[key]]
  }

  metadata <- list(
    timestamp = basename(session_dir),
    session_dir = normalizePath(session_dir, winslash = "/", mustWork = FALSE),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    renv_version = tryCatch(as.character(utils::packageVersion("renv")), error = function(e) NA_character_),
    placebo_rmspe_ratio_cutoff = placebo_rmspe_ratio_cutoff,
    stability = stability_gate,
    env_flags = as.list(env_flags),
    outcomes = vapply(outcomes, `[[`, character(1), "id"),
    specs = vapply(specs, `[[`, character(1), "name"),
    data_files = list(
      list(key = PWT11_CHECKSUM_KEY, path = PWT11_LOCAL, expected_md5 = checksum_lookup(PWT11_CHECKSUM_KEY))
    )
  )
  write_run_metadata(session_dir, metadata)

  max_post_year <- max(vapply(specs, function(s) s$post_window[2], numeric(1)))
  raw <- load_pwt_data()
  panel_full <- prepare_pwt11_panel(raw, max_year = max_post_year)
  panel_full$region <- as.numeric(panel_full$region)
  panel_full <- as.data.frame(panel_full)
  panel_full$country <- as.character(panel_full$country)
  panel_full$region <- as.numeric(panel_full$region)

  treatment_identifier <- "Spain"
  base_key_vars <- c("gdpcap", "hc", "csh_c", "csh_i", "csh_g", "csh_x", "csh_m", "rknacapita", "pl_gdpo")
  build_outcome_key_vars <- function(dep_var) {
    unique(c("gdpcap", setdiff(base_key_vars, dep_var)))
  }

  metrics_accum <- list()
  predictor_balance_tables <- list()
  predictor_balance_rows <- list()
  donor_weights_summary <- list()
  donor_weights_rows <- list()
  pool_status_rows <- list()
  pool_summary_rows <- list()
  placebo_filter_rows <- list()
  stability_selection_rows <- list()
  spec_outputs <- list()
  pval_accum <- list()

  for (outcome in outcomes) {
    key_vars <- build_outcome_key_vars(outcome$dep_var)
    outcome_specs <- Filter(function(s) s$name %in% outcome$specs, specs)
    if (length(outcome_specs) == 0) next

    outcome_root_tables <- file.path(tables_dir, outcome$id)
    outcome_root_plots <- file.path(plots_dir, outcome$id)
    ensure_dir(outcome_root_tables)
    ensure_dir(outcome_root_plots)

    spec_outputs[[outcome$id]] <- list()
    use_stability_gate <- isTRUE(outcome$use_stability_gate)

    for (spec in outcome_specs) {
      spec_tables_root <- file.path(outcome_root_tables, "specs", spec$name)
      spec_plots_root <- file.path(outcome_root_plots, "specs", spec$name)
      ensure_dir(spec_tables_root)
      ensure_dir(spec_plots_root)

      pool_results <- list()
      top_donor <- NULL
      dep_var <- outcome$dep_var

      run_spec_pool <- function(pool_label, donor_filter, run_placebos = TRUE,
                                times_pred_override = NULL, agg_fns_override = NULL) {
        pool_tables <- file.path(spec_tables_root, pool_label)
        pool_plots <- file.path(spec_plots_root, pool_label)
        ensure_dir(pool_tables)
        ensure_dir(pool_plots)

        panel_res <- prepare_panel_for_spec(
          panel_full = panel_full,
          pre_window = spec$pre_window,
          post_window = spec$post_window,
          key_vars = key_vars,
          outcome_var = dep_var,
          spec_name = spec$name,
          pool_label = pool_label,
          tables_dir = pool_tables,
          treatment_identifier = treatment_identifier,
          donor_filter = donor_filter
        )

        key_vars <- panel_res$key_vars
        panel_df <- panel_res$panel
        keep_vars <- unique(c("country", "isocode", "region", "year", dep_var, key_vars))
        panel_df <- panel_df[, intersect(names(panel_df), keep_vars), drop = FALSE]

        datprep <- MSCMT::listFromLong(
          panel_df,
          unit.variable = which(names(panel_df) == "region"),
          time.variable = which(names(panel_df) == "year"),
          unit.names.variable = which(names(panel_df) == "country")
        )

        controls_identifier <- setdiff(colnames(datprep[[1]]), c(treatment_identifier, "Average"))
        pre_period <- range(spec$pre_window)
        post_period <- range(spec$post_window)
        pred_vars <- key_vars

        times_dep <- matrix(pre_period, ncol = 1)
        colnames(times_dep) <- dep_var
        times_pred <- do.call(cbind, lapply(pred_vars, function(x) pre_period))
        colnames(times_pred) <- pred_vars
        agg_fns <- rep("mean", ncol(times_pred))
        # Bateria de predictores: permite sustituir el set de predictores
        # (p. ej. lags del outcome) manteniendo panel, pool y ventanas.
        if (!is.null(times_pred_override)) {
          times_pred <- times_pred_override
          agg_fns <- agg_fns_override
        }

        res <- suppressWarnings(
          MSCMT::mscmt(
            datprep,
            treatment.identifier = treatment_identifier,
            controls.identifier = controls_identifier,
            times.dep = times_dep,
            times.pred = times_pred,
            agg.fns = agg_fns,
            seed = 1,
            outer.optim = "DEoptim",
            outer.opar = deoptim_opar,
            verbose = FALSE
          )
        )

        predictor_table_raw <- as.data.frame(res$predictor.table)
        predictor_table <- build_predictor_balance_table(predictor_table_raw)
        if (!is.null(predictor_table)) {
          predictor_balance_tables[[paste(outcome$id, spec$name, pool_label, sep = "_")]] <<- predictor_table
          predictor_balance_rows[[length(predictor_balance_rows) + 1]] <<- dplyr::mutate(
            predictor_table,
            outcome = outcome$id,
            spec = spec$name,
            pool = pool_label
          ) |>
            dplyr::relocate(outcome, spec, pool)
          write_working_csv(predictor_table, file.path(pool_tables, "predictor_balance.csv"))
          save_table_jpg(
            predictor_table,
            file.path(pool_tables, "predictor_balance.jpg"),
            title = sprintf("Predictor balance table (%s / %s / %s)", outcome$id, spec$name, pool_label)
          )
        }

        w_vec <- extract_weights(res, controls_identifier)
        donor_weights <- NULL
        if (!is.null(w_vec)) {
          donor_weights <- tibble::tibble(
            Country = names(w_vec),
            Weight = as.numeric(w_vec)
          ) |>
            dplyr::filter(Weight > 0) |>
            dplyr::arrange(dplyr::desc(Weight))
          write_working_csv(donor_weights, file.path(pool_tables, "donor_weights.csv"))
          save_table_jpg(donor_weights, file.path(pool_tables, "donor_weights.jpg"), title = sprintf("Donor weights (%s / %s / %s)", outcome$id, spec$name, pool_label))
          donor_weights_summary[[paste(outcome$id, spec$name, pool_label, sep = "_")]] <<- donor_weights
          donor_weights_rows[[length(donor_weights_rows) + 1]] <<- dplyr::mutate(
            donor_weights,
            outcome = outcome$id,
            spec = spec$name,
            pool = pool_label
          ) |>
            dplyr::relocate(outcome, spec, pool)
        }

        if (pool_label == "all" && !is.null(donor_weights) && nrow(donor_weights) > 0) {
          top_donor <<- donor_weights$Country[1]
        }

        gap_mat_main <- extract_gap_matrix(res, dep_var = dep_var)
        if (!is.null(gap_mat_main) && ncol(gap_mat_main) == 1) {
          colnames(gap_mat_main) <- treatment_identifier
        }
        gaps_long_main <- NULL
        if (!is.null(gap_mat_main)) {
          years_guess <- suppressWarnings(as.numeric(rownames(datprep[[1]])))
          gaps_long_main <- gap_matrix_to_long(
            gap_mat_main,
            years = resolve_years(gap_mat_main, years_guess),
            unit_names = colnames(datprep[[1]])
          )
        }

        metrics_main <- summarize_gap_metrics(
          gaps_long = gaps_long_main,
          pre_window = pre_period,
          post_window = post_period,
          treatment_identifier = treatment_identifier,
          spec_name = spec$name,
          pool_label = pool_label,
          outcome = outcome$id
        )
        if (!is.null(metrics_main)) metrics_accum[[length(metrics_accum) + 1]] <<- metrics_main

        x_limits <- c(spec$pre_window[1], spec$post_window[2])
        if (!isTRUE(run_placebos)) {
          return(list(
            mscmt = res,
            placebo = NULL,
            predictor_table = predictor_table,
            donor_weights = donor_weights,
            gaps_main = gaps_long_main,
            placebo_outputs = list(gaps_long = NULL, band_df = NULL, fit_metrics_df = NULL, placebo_only = NULL, good_units = NULL),
            panel_info = panel_res,
            conformal = list(band_df = NULL, pvals_df = NULL)
          ))
        }

        gap_ylim <- NULL
        if (!is.null(gaps_long_main)) {
          g_range <- range(gaps_long_main$gap, na.rm = TRUE)
          g_pad <- diff(g_range) * 0.08
          if (!is.finite(g_pad)) g_pad <- 100
          gap_ylim <- c(g_range[1] - g_pad, g_range[2] + g_pad)
        }
        comp_ylim <- NULL
        if (!is.null(gap_mat_main) && !is.null(res$Y1plot)) {
          comp_range <- range(res$Y1plot, res$Y1plot + gap_mat_main, na.rm = TRUE)
          comp_pad <- diff(comp_range) * 0.05
          if (!is.finite(comp_pad)) comp_pad <- 500
          comp_ylim <- c(comp_range[1] - comp_pad, comp_range[2] + comp_pad)
        }

        p3 <- plot_synth_comparison(res, x_limits = x_limits, y_limits = comp_ylim, y_label = outcome$y_label, y_formatter = outcome$y_formatter, title = sprintf("%s: Spain and Synthetic Spain", outcome$y_label))
        p4 <- plot_synth_gaps(gaps_df = gaps_long_main, treatment_identifier = treatment_identifier, x_limits = x_limits, y_limits = gap_ylim, y_label = paste0(outcome$y_label, " gap"), y_formatter = outcome$y_formatter, treatment_year = post_period[1], title = sprintf("%s: Spain minus Synthetic Spain", outcome$y_label))

        safe_ggsave(file.path(pool_plots, "Figure_3_synth_comparison.png"), p3)
        safe_ggsave(file.path(pool_plots, "Figure_4_synth_gaps.png"), p4)

        cl <- tryCatch(parallel::makeCluster(max(1, min(4, parallel::detectCores()))), error = function(e) NULL)
        resplacebo <- suppressWarnings(
          MSCMT::mscmt(
            datprep,
            treatment.identifier = treatment_identifier,
            controls.identifier = controls_identifier,
            times.dep = times_dep,
            times.pred = times_pred,
            agg.fns = agg_fns,
            cl = cl,
            placebo = TRUE,
            seed = 1,
            outer.optim = "DEoptim",
            outer.opar = deoptim_opar
          )
        )
        if (!is.null(cl)) parallel::stopCluster(cl)

        time_index <- suppressWarnings(as.numeric(rownames(datprep[[1]])))
        if (any(is.na(time_index))) time_index <- sort(unique(panel_df$year))

        placebo_outputs <- build_placebo_outputs(
          resplacebo,
          years = time_index,
          unit_names = colnames(datprep[[1]]),
          treatment_identifier = treatment_identifier,
          pre_period = pre_period,
          post_period = post_period,
          tables_dir = pool_tables,
          dep_var = dep_var
        )

        # Filtro ex ante: descarta placebos con RMSPE pre muy superior al tratado (Abadie et al. 2010)
        if (!is.null(placebo_outputs$gaps_long) && !is.null(placebo_outputs$fit_metrics_df)) {
          pre_rmspe_df <- placebo_outputs$fit_metrics_df |>
            dplyr::select(Country, pre_rmspe)
          treated_rmspe <- pre_rmspe_df$pre_rmspe[pre_rmspe_df$Country == treatment_identifier]
          cutoff <- placebo_rmspe_ratio_cutoff
          use_filter <- length(treated_rmspe) > 0 && is.finite(treated_rmspe) && treated_rmspe > 0 &&
            is.finite(cutoff) && cutoff > 0
          pre_rmspe_df <- pre_rmspe_df |>
            dplyr::mutate(
              ratio_vs_treated = if (use_filter) pre_rmspe / treated_rmspe else NA_real_,
              keep = if (use_filter) (is.na(ratio_vs_treated) | ratio_vs_treated <= cutoff) else TRUE
            )
          good_units <- pre_rmspe_df$Country[pre_rmspe_df$keep | pre_rmspe_df$Country == treatment_identifier | pre_rmspe_df$Country == "Average"]
          good_units <- unique(c(good_units, treatment_identifier))
          placebos_total <- sum(pre_rmspe_df$Country != treatment_identifier & pre_rmspe_df$Country != "Average", na.rm = TRUE)
          placebos_kept <- sum(pre_rmspe_df$Country %in% good_units & pre_rmspe_df$Country != treatment_identifier & pre_rmspe_df$Country != "Average", na.rm = TRUE)
          placebo_excluded <- pre_rmspe_df |>
            dplyr::filter(!keep, Country != treatment_identifier, Country != "Average")
          if (nrow(placebo_excluded) > 0) {
            write_working_csv(placebo_excluded, file.path(pool_tables, "placebo_excluded_badfit.csv"))
            save_table_jpg(placebo_excluded, file.path(pool_tables, "placebo_excluded_badfit.jpg"), title = sprintf("Placebos excluded (pre RMSPE > %sx treated)", cutoff))
          }
          message(sprintf(
            "Placebos retenidos: %d de %d (cutoff pre-RMSPE x tratado = %s, tratado=%.3f)",
            placebos_kept,
            placebos_total,
            if (!use_filter) "disabled" else format(cutoff, digits = 4),
            ifelse(length(treated_rmspe) > 0, treated_rmspe[1], NA_real_)
          ))
          placebo_filter_rows[[length(placebo_filter_rows) + 1]] <<- tibble::tibble(
            outcome = outcome$id,
            spec = spec$name,
            pool = pool_label,
            filter_applied = isTRUE(use_filter),
            cutoff = if (!use_filter) NA_real_ else cutoff,
            treated_pre_rmspe = ifelse(length(treated_rmspe) > 0, treated_rmspe[1], NA_real_),
            placebos_total = placebos_total,
            placebos_kept = placebos_kept
          )
          gaps_long_filt <- placebo_outputs$gaps_long |>
            dplyr::filter(Country %in% good_units)
          placebo_only_filt <- gaps_long_filt |>
            dplyr::filter(Country != treatment_identifier, Country != "Average")
          band_df_filt <- NULL
          if (nrow(placebo_only_filt) > 0) {
            band_df_filt <- placebo_only_filt |>
              dplyr::group_by(year) |>
              dplyr::summarise(
                q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
                q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
                q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
                .groups = "drop"
              )
          }
          fit_metrics_df_filt <- placebo_outputs$fit_metrics_df |>
            dplyr::filter(Country %in% good_units)
          placebo_outputs$gaps_long <- gaps_long_filt
          placebo_outputs$band_df <- band_df_filt
          placebo_outputs$fit_metrics_df <- fit_metrics_df_filt
          placebo_outputs$placebo_only <- placebo_only_filt
          placebo_outputs$good_units <- good_units
          write_working_csv(gaps_long_filt, file.path(pool_tables, "placebo_gaps_long.csv"))
          if (!is.null(band_df_filt)) write_working_csv(band_df_filt, file.path(pool_tables, "placebo_gap_quantiles.csv"))
          write_fit_metric_outputs(fit_metrics_df_filt, pool_tables, treatment_identifier)
        }

        if (is.null(placebo_outputs$gaps_long) || is.null(placebo_outputs$fit_metrics_df)) {
          placebo_filter_rows[[length(placebo_filter_rows) + 1]] <<- tibble::tibble(
            outcome = outcome$id,
            spec = spec$name,
            pool = pool_label,
            filter_applied = FALSE,
            cutoff = NA_real_,
            treated_pre_rmspe = NA_real_,
            placebos_total = NA_real_,
            placebos_kept = NA_real_
          )
        }

        conformal <- conformal_bands_from_placebos(
          treated_gaps = gaps_long_main,
          placebo_long = placebo_outputs$gaps_long,
          alpha = 0.1,
          treatment_identifier = treatment_identifier
        )

        # Escala del eje Y del placebo: fija para gdpcap y dinamica para el resto.
        p5_ylim <- NULL
        p5_breaks <- waiver()
        p5_range <- range(placebo_outputs$gaps_long$gap, na.rm = TRUE)
        max_abs <- suppressWarnings(max(abs(p5_range), na.rm = TRUE))
        if (!is.finite(max_abs)) max_abs <- 0
        if (identical(outcome$id, "gdpcap")) {
          target <- max(7500, ceiling(max_abs / 2500) * 2500)
          if (!is.finite(target) || target <= 0) target <- 7500
          p5_ylim <- c(-target, target)
          p5_breaks <- seq(-target, target, by = 2500)
        } else if (max_abs > 0) {
          target <- max_abs * 1.08
          p5_ylim <- c(-target, target)
          p5_breaks <- scales::pretty_breaks(n = 6)
        }

        p5 <- plot_placebo(
          placebo_gaps_df = placebo_outputs$gaps_long,
          treatment_identifier = treatment_identifier,
          x_limits = x_limits,
          y_limits = p5_ylim,
          y_breaks = p5_breaks,
          y_label = paste0(outcome$y_label, " gap"),
          y_formatter = outcome$y_formatter,
          treatment_year = post_period[1],
          title = sprintf("Placebo Test: %s gap", outcome$y_label)
        )
        safe_ggsave(file.path(pool_plots, "Figure_5_placebo.png"), p5)

        mspe_ratio_df <- placebo_outputs$fit_metrics_df |>
          dplyr::select(Country, mspe_ratio)
        write_working_csv(mspe_ratio_df, file.path(pool_tables, "post_pre_mspe_ratio.csv"))
        save_table_jpg(mspe_ratio_df, file.path(pool_tables, "post_pre_mspe_ratio.jpg"), title = "Post/Pre MSPE ratio")

        rmspe_ratio_df <- placebo_outputs$fit_metrics_df |>
          dplyr::select(Country, rmspe_ratio)
        write_working_csv(rmspe_ratio_df, file.path(pool_tables, "post_pre_rmspe_ratio.csv"))
        save_table_jpg(rmspe_ratio_df, file.path(pool_tables, "post_pre_rmspe_ratio.jpg"), title = "Post/Pre RMSPE ratio")

        pval_df <- tryCatch({
          mspe_cutoff <- if (is.finite(placebo_rmspe_ratio_cutoff) && placebo_rmspe_ratio_cutoff > 0) placebo_rmspe_ratio_cutoff^2 else 5
          vals <- MSCMT::pvalue(resplacebo, exclude.ratio = mspe_cutoff, ratio.type = "mspe", alternative = "less")
          tibble::enframe(vals, name = "Country", value = "p_value_less")
        }, error = function(e) NULL)
        if (!is.null(pval_df)) {
          write_working_csv(pval_df, file.path(pool_tables, "placebo_pvalues.csv"))
          save_table_jpg(pval_df, file.path(pool_tables, "placebo_pvalues.jpg"), title = "Placebo p-values")
        }

        if (!is.null(conformal$band_df)) {
          write_working_csv(conformal$band_df, file.path(pool_tables, "conformal_bands.csv"))
          save_table_jpg(conformal$band_df, file.path(pool_tables, "conformal_bands.jpg"), title = "Conformal bands (placebos espaciales)")
        }

        # Para graficar, priorizamos la banda de placebos (q05/q50/q95);
        # las bandas conformal se guardan en CSV/JPG pero no reemplazan la mediana placebo en la figura.
        band_for_plot <- placebo_outputs$band_df
        if (is.null(band_for_plot)) band_for_plot <- conformal$band_df
        if (!is.null(band_for_plot)) {
          band_for_plot <- dplyr::filter(band_for_plot, year >= post_period[1])
        }
        if (!is.null(band_for_plot)) {
          p4_band <- plot_synth_gaps(
            gaps_df = gaps_long_main,
            treatment_identifier = treatment_identifier,
            band_df = band_for_plot,
            x_limits = x_limits,
            y_limits = NULL,
            y_label = paste0(outcome$y_label, " gap"),
            y_formatter = outcome$y_formatter,
            treatment_year = post_period[1],
            title = sprintf("%s: Spain minus Synthetic Spain", outcome$y_label)
          )
          safe_ggsave(file.path(pool_plots, "Figure_4_synth_gaps.png"), p4_band)
        }

        if (!is.null(placebo_outputs$band_df)) {
          save_table_jpg(placebo_outputs$band_df, file.path(pool_tables, "placebo_gap_quantiles.jpg"), title = "Placebo gap quantiles")
        }
        if (!is.null(placebo_outputs$fit_metrics_df)) {
          save_table_jpg(
            placebo_outputs$fit_metrics_df |>
              dplyr::filter(Country != treatment_identifier) |>
              dplyr::select(Country, pre_mspe, post_mspe, mspe_ratio),
            file.path(pool_tables, "placebo_mspe.jpg"),
            title = "Placebo MSPE pre/post"
          )
          save_table_jpg(
            placebo_outputs$fit_metrics_df |>
              dplyr::filter(Country != treatment_identifier) |>
              dplyr::select(Country, pre_rmspe, post_rmspe, rmspe_ratio),
            file.path(pool_tables, "placebo_rmspe.jpg"),
            title = "Placebo RMSPE pre/post"
          )
        }
        if (!is.null(placebo_outputs$gaps_long)) {
          save_table_jpg(placebo_outputs$gaps_long, file.path(pool_tables, "placebo_gaps_long.jpg"), title = "Placebo gaps (truncated)")
        }

        if (!is.null(placebo_outputs$gaps_long)) {
          post_avg_gaps <- placebo_outputs$gaps_long |>
            dplyr::filter(Country != "Average") |>
            dplyr::filter(year >= post_period[1], year <= post_period[2]) |>
            dplyr::group_by(Country) |>
            dplyr::summarise(avg_gap_post = mean(gap, na.rm = TRUE), .groups = "drop") |>
            dplyr::mutate(
              is_treated = Country == treatment_identifier,
              is_spain = Country == "Spain"
            )
          write_working_csv(post_avg_gaps, file.path(pool_tables, "placebo_post_avg_gaps.csv"))
          save_table_jpg(post_avg_gaps, file.path(pool_tables, "placebo_post_avg_gaps.jpg"), title = "Average post-treatment gaps (all units)")

          p_gap <- NA_real_
          if (any(post_avg_gaps$is_treated)) {
            treated_gap <- post_avg_gaps$avg_gap_post[post_avg_gaps$is_treated]
            placebo_gaps <- post_avg_gaps$avg_gap_post[!post_avg_gaps$is_treated]
            p_gap <- rank_pvalue(treated_gap, placebo_gaps, direction = "greater")
            p_gap_df <- tibble::tibble(metric = "avg_post_gap", p_value = p_gap, treated = treated_gap, n_placebos = length(placebo_gaps))
            write_working_csv(p_gap_df, file.path(pool_tables, "pvalue_post_gap.csv"))
            save_table_jpg(p_gap_df, file.path(pool_tables, "pvalue_post_gap.jpg"), title = "p-value (avg post gap, ranking placebo)")
            pval_accum[[length(pval_accum) + 1]] <<- tibble::tibble(
              outcome = outcome$id,
              spec = spec$name,
              pool = pool_label,
              metric = "avg_post_gap",
              p_value_raw = p_gap
            )
          }

          if (nrow(post_avg_gaps) > 0) {
            p_gap_hist <- ggplot2::ggplot(post_avg_gaps, ggplot2::aes(x = avg_gap_post, fill = is_spain)) +
              ggplot2::geom_histogram(alpha = 0.7, bins = 30, colour = "white") +
              ggplot2::geom_vline(
                data = dplyr::filter(post_avg_gaps, is_spain),
                ggplot2::aes(xintercept = avg_gap_post),
                colour = "red",
                linewidth = 1
              ) +
              ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
              ggplot2::labs(x = "Average gap (post period)", y = "Count", title = "Distribution of post-treatment gaps (placebos vs treated)") +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_post_gap_hist.png"), p_gap_hist)

            p_gap_ecdf <- ggplot2::ggplot(post_avg_gaps, ggplot2::aes(x = avg_gap_post, colour = is_spain)) +
              ggplot2::stat_ecdf(linewidth = 0.9) +
              ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
              ggplot2::labs(x = "Average gap (post period)", y = "ECDF", title = "ECDF of post-treatment gaps (placebos vs treated)") +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_post_gap_ecdf.png"), p_gap_ecdf)

            rank_df <- post_avg_gaps |>
              dplyr::arrange(avg_gap_post) |>
              dplyr::mutate(rank = dplyr::row_number())
            p_rank <- ggplot2::ggplot(rank_df, ggplot2::aes(x = rank, y = avg_gap_post, colour = is_spain)) +
              ggplot2::geom_point(size = 2) +
              ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
              ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey30")) +
              ggplot2::labs(x = "Rank (ascending)", y = "Average gap (post)", title = "Ranked post-treatment gaps", subtitle = if (!is.na(p_gap)) sprintf("Placebo p-value (greater): %.3f", p_gap) else NULL) +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_post_gap_rank.png"), p_rank)
          }
        }

        if (!is.null(placebo_outputs$fit_metrics_df)) {
          rmspe_df_plot <- placebo_outputs$fit_metrics_df |>
            dplyr::mutate(is_treated = Country == treatment_identifier) |>
            dplyr::filter(is.finite(rmspe_ratio))
          write_working_csv(rmspe_df_plot, file.path(pool_tables, "unit_rmspe_ratios.csv"))
          if (nrow(rmspe_df_plot) > 0 && any(rmspe_df_plot$is_treated)) {
            treated_rmspe_ratio <- rmspe_df_plot$rmspe_ratio[rmspe_df_plot$is_treated]
            placebo_rmspe_ratio <- rmspe_df_plot$rmspe_ratio[!rmspe_df_plot$is_treated]
            p_rmspe <- rank_pvalue(treated_rmspe_ratio, placebo_rmspe_ratio, direction = "greater")
            p_rmspe_df <- tibble::tibble(metric = "rmspe_ratio", p_value = p_rmspe, treated = treated_rmspe_ratio, n_placebos = length(placebo_rmspe_ratio))
            write_working_csv(p_rmspe_df, file.path(pool_tables, "pvalue_rmspe_ratio.csv"))
            save_table_jpg(p_rmspe_df, file.path(pool_tables, "pvalue_rmspe_ratio.jpg"), title = "p-value (post/pre RMSPE ratio)")
            pval_accum[[length(pval_accum) + 1]] <<- tibble::tibble(
              outcome = outcome$id,
              spec = spec$name,
              pool = pool_label,
              metric = "rmspe_ratio",
              p_value_raw = p_rmspe
            )
            p_rmspe_hist <- ggplot2::ggplot(rmspe_df_plot, ggplot2::aes(x = rmspe_ratio, fill = is_treated)) +
              ggplot2::geom_histogram(alpha = 0.7, bins = 30, colour = "white") +
              ggplot2::geom_vline(
                data = dplyr::filter(rmspe_df_plot, is_treated),
                ggplot2::aes(xintercept = rmspe_ratio),
                colour = "red",
                linewidth = 1
              ) +
              ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
              ggplot2::labs(x = "Post/Pre RMSPE ratio", y = "Count", title = "Distribution of post/pre RMSPE ratios") +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_rmspe_ratio_hist.png"), p_rmspe_hist)
            p_rmspe_ecdf <- ggplot2::ggplot(rmspe_df_plot, ggplot2::aes(x = rmspe_ratio, colour = is_treated)) +
              ggplot2::stat_ecdf(linewidth = 0.9) +
              ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
              ggplot2::labs(x = "Post/Pre RMSPE ratio", y = "ECDF", title = "ECDF of post/pre RMSPE ratios") +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_rmspe_ratio_ecdf.png"), p_rmspe_ecdf)
            rank_rmspe <- rmspe_df_plot |>
              dplyr::arrange(rmspe_ratio) |>
              dplyr::mutate(rank = dplyr::row_number())
            p_rmspe_rank <- ggplot2::ggplot(rank_rmspe, ggplot2::aes(x = rank, y = rmspe_ratio, colour = is_treated)) +
              ggplot2::geom_point(size = 2) +
              ggplot2::scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey30")) +
              ggplot2::labs(x = "Rank (ascending)", y = "Post/Pre RMSPE ratio", title = "Ranked post/pre RMSPE ratios", subtitle = if (!is.na(p_rmspe)) sprintf("Placebo p-value (greater): %.3f", p_rmspe) else NULL) +
              ggplot2::theme_minimal()
            safe_ggsave(file.path(pool_plots, "placebo_rmspe_ratio_rank.png"), p_rmspe_rank)
          }
        }

        # Resumen de RMSPE y tendencias pre
        if (!is.null(placebo_outputs$gaps_long)) {
          # Bandas placebo en el pre
          pre_band <- placebo_outputs$placebo_only |>
            dplyr::filter(year >= pre_period[1], year <= pre_period[2]) |>
            dplyr::group_by(year) |>
            dplyr::summarise(
              q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
              q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
              q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
              .groups = "drop"
            )
          if (nrow(pre_band) > 0 && !is.null(gaps_long_main)) {
            write_working_csv(pre_band, file.path(pool_tables, "pre_placebo_gap_quantiles.csv"))
            save_table_jpg(pre_band, file.path(pool_tables, "pre_placebo_gap_quantiles.jpg"), title = "Pre-treatment placebo bands (gaps)")
            treated_pre <- gaps_long_main |>
              dplyr::filter(year >= pre_period[1], year <= pre_period[2], Country == treatment_identifier)
            if (nrow(treated_pre) > 0) {
              pre_band$year_num <- as.integer(round(pre_band$year))
              treated_pre$year_num <- as.integer(round(treated_pre$year))
              p_pre <- ggplot2::ggplot() +
                ggplot2::geom_ribbon(
                  data = pre_band,
                  ggplot2::aes(x = year_num, ymin = q05, ymax = q95),
                  fill = "grey70",
                  alpha = 0.3
                ) +
                ggplot2::geom_line(
                  data = pre_band,
                  ggplot2::aes(x = year_num, y = q50),
                  linetype = "dashed",
                  colour = "black"
                ) +
                ggplot2::geom_line(
                  data = treated_pre,
                  ggplot2::aes(x = year_num, y = gap),
                  colour = "red",
                  linewidth = 1
                ) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
                ggplot2::labs(title = "Pre-treatment gaps with placebo bands", x = "Year", y = "Gap") +
                ggplot2::theme_minimal()
              safe_ggsave(file.path(pool_plots, "pre_gap_bands.png"), p_pre)
            }
          }

          # Bandas placebo en el post
          post_band <- placebo_outputs$placebo_only |>
            dplyr::filter(year >= post_period[1], year <= post_period[2]) |>
            dplyr::group_by(year) |>
            dplyr::summarise(
              q05 = stats::quantile(gap, 0.05, na.rm = TRUE),
              q50 = stats::quantile(gap, 0.50, na.rm = TRUE),
              q95 = stats::quantile(gap, 0.95, na.rm = TRUE),
              .groups = "drop"
            )
          if (nrow(post_band) > 0 && !is.null(gaps_long_main)) {
            write_working_csv(post_band, file.path(pool_tables, "post_placebo_gap_quantiles.csv"))
            save_table_jpg(post_band, file.path(pool_tables, "post_placebo_gap_quantiles.jpg"), title = "Post-treatment placebo bands (gaps)")
            treated_post <- gaps_long_main |>
              dplyr::filter(year >= post_period[1], year <= post_period[2], Country == treatment_identifier)
            if (nrow(treated_post) > 0) {
              post_band$year_num <- as.integer(round(post_band$year))
              treated_post$year_num <- as.integer(round(treated_post$year))
              p_post <- ggplot2::ggplot() +
                ggplot2::geom_ribbon(
                  data = post_band,
                  ggplot2::aes(x = year_num, ymin = q05, ymax = q95),
                  fill = "grey70",
                  alpha = 0.3
                ) +
                ggplot2::geom_line(
                  data = post_band,
                  ggplot2::aes(x = year_num, y = q50),
                  linetype = "dashed",
                  colour = "black"
                ) +
                ggplot2::geom_line(
                  data = treated_post,
                  ggplot2::aes(x = year_num, y = gap),
                  colour = "red",
                  linewidth = 1
                ) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
                ggplot2::labs(title = "Post-treatment gaps with placebo bands", x = "Year", y = "Gap") +
                ggplot2::theme_minimal()
              safe_ggsave(file.path(pool_plots, "post_gap_bands.png"), p_post)
            }
          }

          # Test de tendencia pre (pendiente)
          slope_fun <- function(df) {
            if (nrow(df) < 2) return(NA_real_)
            fit <- try(stats::lm(gap ~ year, data = df), silent = TRUE)
            if (inherits(fit, "try-error")) return(NA_real_)
            as.numeric(coef(fit)[["year"]])
          }
          slopes <- placebo_outputs$gaps_long |>
            dplyr::filter(year >= pre_period[1], year <= pre_period[2]) |>
            dplyr::group_by(Country) |>
            dplyr::summarise(
              n_pre = dplyr::n(),
              slope = slope_fun(dplyr::cur_data_all()),
              .groups = "drop"
            )
          if (nrow(slopes) > 0) {
            write_working_csv(slopes, file.path(pool_tables, "pre_gap_slopes.csv"))
            treated_slope <- slopes$slope[slopes$Country == treatment_identifier]
            placebo_slopes <- slopes$slope[slopes$Country != treatment_identifier & slopes$Country != "Average"]
            p_slope <- rank_pvalue(abs(treated_slope), abs(placebo_slopes), direction = "greater")
            slope_summary <- tibble::tibble(
              metric = "pre_gap_slope",
              treated_slope = treated_slope,
              p_value = p_slope,
              n_placebos = length(placebo_slopes)
            )
            write_working_csv(slope_summary, file.path(pool_tables, "pre_gap_slope_summary.csv"))
            save_table_jpg(slope_summary, file.path(pool_tables, "pre_gap_slope_summary.jpg"), title = "Pre-treatment gap slope (rank placebo p-value)")
          }
        }

        if (!is.null(placebo_outputs$fit_metrics_df)) {
          # Resumen RMSPE pre/post con percentiles placebo
          fit_metrics_df_plot <- placebo_outputs$fit_metrics_df |>
            dplyr::mutate(is_treated = Country == treatment_identifier) |>
            dplyr::filter(is.finite(pre_rmspe) | is.finite(post_rmspe))
          treated_pre_rmspe <- fit_metrics_df_plot$pre_rmspe[fit_metrics_df_plot$is_treated]
          treated_post_rmspe <- fit_metrics_df_plot$post_rmspe[fit_metrics_df_plot$is_treated]
          treated_ratio <- fit_metrics_df_plot$rmspe_ratio[fit_metrics_df_plot$is_treated]
          treated_pre_rmspe_val <- if (length(treated_pre_rmspe) > 0) treated_pre_rmspe[1] else NA_real_
          treated_post_rmspe_val <- if (length(treated_post_rmspe) > 0) treated_post_rmspe[1] else NA_real_
          treated_ratio_val <- if (length(treated_ratio) > 0) treated_ratio[1] else NA_real_
          placebo_ratios <- fit_metrics_df_plot$rmspe_ratio[!fit_metrics_df_plot$is_treated]
          placebo_pre_vals <- fit_metrics_df_plot$pre_rmspe[!fit_metrics_df_plot$is_treated]
          safe_quant <- function(x) {
            if (length(x) == 0 || all(is.na(x))) return(rep(NA_real_, 3))
            stats::quantile(x, probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
          }
          ratio_quant <- safe_quant(placebo_ratios)
          pre_quant <- safe_quant(placebo_pre_vals)
          rmspe_summary <- tibble::tibble(
            row = c("treated", "placebo_q05", "placebo_q50", "placebo_q95"),
            pre_rmspe = c(treated_pre_rmspe_val, pre_quant),
            post_rmspe = c(treated_post_rmspe_val, NA, NA, NA),
            rmspe_ratio = c(treated_ratio_val, ratio_quant)
          )
          p_rmspe_val <- rank_pvalue(treated_ratio_val, placebo_ratios, direction = "greater")
          p_rmspe_summary <- tibble::tibble(
            metric = "rmspe_ratio_rank",
            treated = treated_ratio_val,
            p_value = p_rmspe_val,
            n_placebos = length(placebo_ratios)
          )
          write_working_csv(rmspe_summary, file.path(pool_tables, "rmspe_summary.csv"))
          save_table_jpg(rmspe_summary, file.path(pool_tables, "rmspe_summary.jpg"), title = "RMSPE pre/post summary")
          write_working_csv(p_rmspe_summary, file.path(pool_tables, "rmspe_pvalue_summary.csv"))
          save_table_jpg(p_rmspe_summary, file.path(pool_tables, "rmspe_pvalue_summary.jpg"), title = "RMSPE ratio p-value (placebo rank)")
        }

        p6 <- plot_post_pre_ratio(
          ratio_df = rmspe_ratio_df,
          value_col = "rmspe_ratio",
          metric_label = "RMSPE",
          treatment_identifier = treatment_identifier
        )
        safe_ggsave(file.path(pool_plots, "Figure_6_rmspe_ratio.png"), p6)
        save_working_rds(resplacebo, file.path(pool_tables, "resplacebo_object.rds"))

        metrics_placebo <- summarize_gap_metrics(
          gaps_long = placebo_outputs$gaps_long,
          pre_window = pre_period,
          post_window = post_period,
          treatment_identifier = treatment_identifier,
          spec_name = spec$name,
          pool_label = pool_label,
          outcome = outcome$id
        )
        if (!is.null(metrics_placebo)) {
          metrics_placebo$metric_source <- "placebo_gaps"
          metrics_accum[[length(metrics_accum) + 1]] <<- metrics_placebo
        }

        list(
          mscmt = res,
          placebo = resplacebo,
          predictor_table = predictor_table,
          donor_weights = donor_weights,
          gaps_main = gaps_long_main,
          placebo_outputs = placebo_outputs,
          panel_info = panel_res,
          conformal = conformal
        )
      }

      pool_stats <- function(pool_obj) {
        w_stats <- compute_weight_stability(pool_obj$donor_weights)
        g_stats <- compute_treated_gap_stats(
          gaps_long_main = pool_obj$gaps_main,
          treatment_identifier = treatment_identifier,
          pre_window = spec$pre_window,
          post_window = spec$post_window
        )
        c(w_stats, g_stats)
      }

      add_stability_row <- function(candidate_pool, candidate_source, dropped_donor, stats, tau_change_pct, eval, selected, reason, pre_mspe_ref = NA_real_) {
        max_pre_allowed <- if (is.finite(pre_mspe_ref) && pre_mspe_ref > 0) {
          pre_mspe_ref * (1 + stability_gate$pre_mspe_max_increase_pct / 100)
        } else {
          NA_real_
        }
        stability_selection_rows[[length(stability_selection_rows) + 1]] <<- tibble::tibble(
          outcome = outcome$id,
          spec = spec$name,
          candidate_pool = candidate_pool,
          candidate_source = candidate_source,
          dropped_donor = dropped_donor,
          top_donor = stats$top_donor,
          top_weight = stats$top_weight,
          n_eff = stats$n_eff,
          n_positive_donors = stats$n_positive,
          treated_pre_mspe = stats$pre_mspe,
          treated_tau_post = stats$tau_post,
          tau_change_pct = tau_change_pct,
          pre_mspe_ref = pre_mspe_ref,
          pre_mspe_max_allowed = max_pre_allowed,
          pass_top_weight = eval$pass_top_weight,
          pass_neff = eval$pass_neff,
          pass_n_positive = eval$pass_n_positive,
          pass_tau = eval$pass_tau,
          pass_pre_cap = eval$pass_pre_cap,
          eligible = eval$eligible,
          selected = selected,
          decision_reason = reason
        )
      }

      apply_stability_gate <- use_stability_gate && isTRUE(stability_gate$enabled) && spec$name %in% stability_gate$specs
      if (apply_stability_gate) {
        pool_results[["all_unconstrained"]] <- run_spec_pool(pool_label = "all_unconstrained", donor_filter = NULL, run_placebos = FALSE)

        all_stats <- pool_stats(pool_results[["all_unconstrained"]])
        all_top_donor <- all_stats$top_donor
        tau_change_top1 <- NA_real_
        top1_drop_label <- NA_character_
        if (!is.na(all_top_donor) && nzchar(all_top_donor)) {
          top1_drop_label <- sprintf("drop_%s", all_top_donor)
          if (is.null(pool_results[[top1_drop_label]])) {
            drop_countries <- setdiff(unique(panel_full$country), all_top_donor)
            pool_results[[top1_drop_label]] <- run_spec_pool(pool_label = top1_drop_label, donor_filter = drop_countries, run_placebos = FALSE)
          }
          top1_stats <- pool_stats(pool_results[[top1_drop_label]])
          if (is.finite(all_stats$tau_post) && abs(all_stats$tau_post) > .Machine$double.eps && is.finite(top1_stats$tau_post)) {
            tau_change_top1 <- abs(top1_stats$tau_post - all_stats$tau_post) / abs(all_stats$tau_post) * 100
          }
        }

        all_eval <- evaluate_stability_candidate(
          stats = all_stats,
          thresholds = stability_gate,
          tau_change_pct = tau_change_top1,
          pre_mspe_ref = NA_real_,
          require_tau = TRUE,
          require_pre_cap = FALSE
        )
        add_stability_row(
          candidate_pool = "all_unconstrained",
          candidate_source = "all_unconstrained",
          dropped_donor = NA_character_,
          stats = all_stats,
          tau_change_pct = tau_change_top1,
          eval = all_eval,
          selected = FALSE,
          reason = if (all_eval$eligible) "all_unconstrained_pass" else "all_unconstrained_fail"
        )

        selected_candidate <- list(
          source = "all_unconstrained",
          dropped_donor = NA_character_,
          donor_filter = NULL,
          stats = all_stats,
          eval = all_eval,
          tau_change_pct = tau_change_top1,
          reason = "accepted_unconstrained"
        )

        if (!isTRUE(all_eval$eligible)) {
          donor_rank <- pool_results[["all_unconstrained"]]$donor_weights$Country
          donor_rank <- donor_rank[!is.na(donor_rank) & nzchar(donor_rank)]
          top_k <- suppressWarnings(as.integer(stability_gate$top_k_candidates))
          if (!is.finite(top_k) || top_k <= 0) top_k <- 3L
          candidate_donors <- unique(utils::head(donor_rank, n = min(length(donor_rank), top_k)))
          feasible_candidates <- list()

          for (cand_donor in candidate_donors) {
            cand_label <- sprintf("drop_%s", cand_donor)
            if (is.null(pool_results[[cand_label]])) {
              drop_countries <- setdiff(unique(panel_full$country), cand_donor)
              pool_results[[cand_label]] <- run_spec_pool(pool_label = cand_label, donor_filter = drop_countries, run_placebos = FALSE)
            }
            cand_stats <- pool_stats(pool_results[[cand_label]])
            cand_eval <- evaluate_stability_candidate(
              stats = cand_stats,
              thresholds = stability_gate,
              tau_change_pct = NA_real_,
              pre_mspe_ref = all_stats$pre_mspe,
              require_tau = FALSE,
              require_pre_cap = TRUE
            )
            add_stability_row(
              candidate_pool = cand_label,
              candidate_source = cand_label,
              dropped_donor = cand_donor,
              stats = cand_stats,
              tau_change_pct = NA_real_,
              eval = cand_eval,
              selected = FALSE,
              reason = if (cand_eval$eligible) "fallback_candidate_pass" else "fallback_candidate_fail",
              pre_mspe_ref = all_stats$pre_mspe
            )
            if (isTRUE(cand_eval$eligible)) {
              feasible_candidates[[length(feasible_candidates) + 1]] <- list(
                source = cand_label,
                dropped_donor = cand_donor,
                donor_filter = setdiff(unique(panel_full$country), cand_donor),
                stats = cand_stats,
                eval = cand_eval,
                tau_change_pct = NA_real_,
                reason = "fallback_drop_topk"
              )
            }
          }

          if (length(feasible_candidates) == 0) {
            stop(sprintf(
              "No stable baseline pool found for outcome=%s spec=%s under configured thresholds.",
              outcome$id,
              spec$name
            ))
          }

          pre_vals <- vapply(feasible_candidates, function(x) x$stats$pre_mspe, numeric(1))
          pre_vals[!is.finite(pre_vals)] <- Inf
          best_idx <- which.min(pre_vals)
          selected_candidate <- feasible_candidates[[best_idx]]
        }

        pool_results[["all"]] <- run_spec_pool(pool_label = "all", donor_filter = selected_candidate$donor_filter)
        all_final_stats <- pool_stats(pool_results[["all"]])
        all_final_eval <- evaluate_stability_candidate(
          stats = all_final_stats,
          thresholds = stability_gate,
          tau_change_pct = if (identical(selected_candidate$source, "all_unconstrained")) selected_candidate$tau_change_pct else NA_real_,
          pre_mspe_ref = if (identical(selected_candidate$source, "all_unconstrained")) NA_real_ else all_stats$pre_mspe,
          require_tau = identical(selected_candidate$source, "all_unconstrained"),
          require_pre_cap = !identical(selected_candidate$source, "all_unconstrained")
        )
        add_stability_row(
          candidate_pool = "all",
          candidate_source = selected_candidate$source,
          dropped_donor = selected_candidate$dropped_donor,
          stats = all_final_stats,
          tau_change_pct = if (identical(selected_candidate$source, "all_unconstrained")) selected_candidate$tau_change_pct else NA_real_,
          eval = all_final_eval,
          selected = TRUE,
          reason = selected_candidate$reason,
          pre_mspe_ref = if (identical(selected_candidate$source, "all_unconstrained")) NA_real_ else all_stats$pre_mspe
        )
        message(sprintf(
          "Stability gate selected pool for %s/%s: source=%s dropped=%s top_weight=%.3f n_eff=%.2f",
          outcome$id,
          spec$name,
          selected_candidate$source,
          ifelse(is.na(selected_candidate$dropped_donor), "none", selected_candidate$dropped_donor),
          ifelse(is.finite(all_final_stats$top_weight), all_final_stats$top_weight, NA_real_),
          ifelse(is.finite(all_final_stats$n_eff), all_final_stats$n_eff, NA_real_)
        ))
      } else {
        pool_results[["all"]] <- run_spec_pool(pool_label = "all", donor_filter = NULL)
      }

      # Bateria de predictores (Ferman-Pinto-Possebom): mismo pool que el
      # baseline seleccionado, variando solo el set de predictores. Solo para
      # gdpcap/baseline; alimenta la Tabla A3.
      pred_battery_on <- if (exists("predictor_battery_enabled", inherits = TRUE)) isTRUE(predictor_battery_enabled) else FALSE
      if (pred_battery_on && identical(outcome$id, "gdpcap") && identical(spec$name, "baseline")) {
        baseline_donor_filter <- if (apply_stability_gate) selected_candidate$donor_filter else NULL
        pre_years <- seq(spec$pre_window[1], spec$pre_window[2])
        mid_year <- spec$pre_window[1] + ceiling((spec$pre_window[2] - spec$pre_window[1]) / 2)
        battery_defs <- list(
          predspec_all_lags = pre_years,
          predspec_odd_lags = pre_years[seq(2, length(pre_years), by = 2)],
          predspec_first_mid_last = unique(c(spec$pre_window[1], mid_year, spec$pre_window[2]))
        )
        for (bat_label in names(battery_defs)) {
          if (!is.null(pool_results[[bat_label]])) next
          bat_years <- battery_defs[[bat_label]]
          bat_times_pred <- do.call(cbind, lapply(bat_years, function(y) c(y, y)))
          colnames(bat_times_pred) <- rep(dep_var, length(bat_years))
          pool_results[[bat_label]] <- run_spec_pool(
            pool_label = bat_label,
            donor_filter = baseline_donor_filter,
            times_pred_override = bat_times_pred,
            agg_fns_override = rep("mean", ncol(bat_times_pred))
          )
        }
      }

      # Pools de donantes restringidos (Europa / Europa sin Portugal / LatAm):
      # mismos predictores y ventanas que el baseline, cambia solo el pool.
      # Solo para gdpcap/baseline; alimenta la Tabla A4.
      restricted_on <- if (exists("restricted_pools_enabled", inherits = TRUE)) isTRUE(restricted_pools_enabled) else FALSE
      if (restricted_on && identical(outcome$id, "gdpcap") && identical(spec$name, "baseline")) {
        restricted_defs <- list(
          pool_europe = europe_donor_countries,
          pool_europe_no_portugal = setdiff(europe_donor_countries, "Portugal"),
          pool_latam = latam_donor_countries
        )
        for (rp_label in names(restricted_defs)) {
          if (!is.null(pool_results[[rp_label]])) next
          pool_results[[rp_label]] <- run_spec_pool(
            pool_label = rp_label,
            donor_filter = restricted_defs[[rp_label]]
          )
        }
        # Variante "pool completo sin gate": los 36 donantes (EEUU incluido),
        # sin la seleccion del stability gate. Documenta que el gate liga por
        # concentracion del donante principal, no por EEUU (peso ~0.06).
        if (is.null(pool_results[["pool_full_no_gate"]])) {
          pool_results[["pool_full_no_gate"]] <- run_spec_pool(
            pool_label = "pool_full_no_gate",
            donor_filter = NULL
          )
        }
      }

      donor_weights_all <- pool_results[["all"]]$donor_weights
      drop_env <- tolower(Sys.getenv("SKIP_DROP_ONE", "0"))
      skip_drop_one <- drop_env %in% c("1", "true", "yes")
      run_drop_one <- isTRUE(outcome$run_drop_one) && identical(spec$name, "baseline") && !skip_drop_one
      donor_pos <- NULL
      if (!is.null(donor_weights_all) && nrow(donor_weights_all) > 0) {
        donor_pos <- unique(donor_weights_all$Country[donor_weights_all$Weight > 0])
      }
      drop_one_exclusions <- character(0)
      if (exists("drop_one_excluded_donors", inherits = TRUE) && length(drop_one_excluded_donors) > 0) {
        drop_one_exclusions <- as.character(drop_one_excluded_donors)
      }
      donor_pos <- setdiff(donor_pos, drop_one_exclusions)
      if (run_drop_one && !is.null(donor_pos) && length(donor_pos) > 0) {
        # Leave-one-out relativo al pool seleccionado: si el stability gate
        # restringio el pool baseline, los drop-one heredan esa restriccion.
        # Las corridas de candidatos del gate (sin placebos) no se reutilizan:
        # se re-estiman con placebos para que la inferencia de la Tabla 6 exista.
        base_pool_countries <- unique(panel_full$country)
        if (apply_stability_gate && !is.null(selected_candidate$donor_filter)) {
          base_pool_countries <- selected_candidate$donor_filter
        }
        for (donor in donor_pos) {
          drop_label <- sprintf("drop_%s", donor)
          existing <- pool_results[[drop_label]]
          if (!is.null(existing) && !is.null(existing$placebo)) next
          drop_countries <- setdiff(base_pool_countries, donor)
          pool_results[[drop_label]] <- run_spec_pool(pool_label = drop_label, donor_filter = drop_countries)
        }
      }


      spec_outputs[[outcome$id]][[spec$name]] <- pool_results
    }
  }

  if (length(metrics_accum) > 0) {
    metrics_df <- dplyr::bind_rows(metrics_accum)
    metrics_path <- file.path(tables_dir, "specs", "gap_metrics_summary.csv")
    ensure_dir(dirname(metrics_path))
    write_working_csv(metrics_df, metrics_path)
    save_table_jpg(metrics_df, file.path(tables_dir, "specs", "gap_metrics_summary.jpg"), title = "Gap/MSPE summary by outcome/spec/pool")
  }

  # Resumen de p-values crudos y ajustados (familia principal: baseline/pool all)
  if (length(pval_accum) > 0) {
    pvals_df <- dplyr::bind_rows(pval_accum)
    pvals_df$p_value_adj <- NA_real_
    family_mask <- pvals_df$spec == "baseline" & pvals_df$pool == "all"
    if (any(family_mask, na.rm = TRUE)) {
      pvals_df$p_value_adj[family_mask] <- stats::p.adjust(pvals_df$p_value_raw[family_mask], method = "holm")
    }
    pvals_path <- file.path(tables_dir, "specs", "pvalues_summary.csv")
    ensure_dir(dirname(pvals_path))
    write_working_csv(pvals_df, pvals_path)
    save_table_jpg(pvals_df, file.path(tables_dir, "specs", "pvalues_summary.jpg"), title = "P-values (raw and Holm-adjusted for baseline/all)")
  }

  specs_dir <- file.path(tables_dir, "specs")
  ensure_dir(specs_dir)

  if (length(pool_status_rows) > 0) {
    pool_status_df <- dplyr::bind_rows(pool_status_rows)
    pool_status_path <- file.path(specs_dir, "pool_status_all.csv")
    write_working_csv(pool_status_df, pool_status_path)
    save_table_jpg(pool_status_df, file.path(specs_dir, "pool_status_all.jpg"), title = "Pool status (all outcomes/specs)")
  }
  if (length(pool_summary_rows) > 0) {
    pool_summary_df <- dplyr::bind_rows(pool_summary_rows)
    pool_summary_path <- file.path(specs_dir, "pool_summary.csv")
    write_working_csv(pool_summary_df, pool_summary_path)
    save_table_jpg(pool_summary_df, file.path(specs_dir, "pool_summary.jpg"), title = "Pool summary by outcome/spec/pool")
  }
  if (length(donor_weights_rows) > 0) {
    donor_weights_df <- dplyr::bind_rows(donor_weights_rows)
    donor_weights_path <- file.path(specs_dir, "donor_weights_summary.csv")
    write_working_csv(donor_weights_df, donor_weights_path)
    save_table_jpg(donor_weights_df, file.path(specs_dir, "donor_weights_summary.jpg"), title = "Donor weights (all specs/pools)")
  }
  if (length(predictor_balance_rows) > 0) {
    predictor_balance_df <- dplyr::bind_rows(predictor_balance_rows)
    predictor_balance_path <- file.path(specs_dir, "predictor_balance_summary.csv")
    write_working_csv(predictor_balance_df, predictor_balance_path)
    save_table_jpg(predictor_balance_df, file.path(specs_dir, "predictor_balance_summary.jpg"), title = "Predictor balance (all specs/pools)")
  }
  if (length(placebo_filter_rows) > 0) {
    placebo_filter_df <- dplyr::bind_rows(placebo_filter_rows)
    placebo_filter_path <- file.path(specs_dir, "placebo_filter_summary.csv")
    write_working_csv(placebo_filter_df, placebo_filter_path)
    save_table_jpg(placebo_filter_df, file.path(specs_dir, "placebo_filter_summary.jpg"), title = "Placebo filter summary (RMSPE cutoffs)")
  }
  stability_df <- NULL
  if (length(stability_selection_rows) > 0) {
    stability_df <- dplyr::bind_rows(stability_selection_rows)
    stability_path <- file.path(specs_dir, "stability_selection_summary.csv")
    write_working_csv(stability_df, stability_path)
    save_table_jpg(stability_df, file.path(specs_dir, "stability_selection_summary.jpg"), title = "Stability gate selection summary")
  }

  export_journal_outputs(
    session_dir = session_dir,
    plots_dir = plots_dir,
    tables_dir = tables_dir,
    spec_outputs = spec_outputs,
    specs = specs,
    treatment_identifier = treatment_identifier,
    stability_selection = stability_df
  )

  list(
    specs = spec_outputs,
    predictor_tables = predictor_balance_tables,
    predictor_balance_tables = predictor_balance_tables,
    donor_weights = donor_weights_summary
  )
}







