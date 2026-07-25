
prepare_pwt11_panel <- function(raw_df,
                                countries = synthetic_countries,
                                code_map = region_codes,
                                max_year = NULL) {
  # Normaliza nombres ANTES de filtrar: PWT usa nombres oficiales largos para
  # Bolivia/Venezuela, "Turkiye" (isocode TUR) y "United States" (isocode USA);
  # si se filtra primero, estos paises quedan silenciosamente fuera del pool.
  df <- raw_df |>
    dplyr::mutate(
      country = as.character(country),
      country = textclean::mgsub(
        country,
        c("Bolivia (Plurinational State of)", "Venezuela (Bolivarian Republic of)"),
        c("Bolivia", "Venezuela")
      ),
      country = dplyr::if_else(isocode == "TUR", "Turkey", country),
      country = dplyr::if_else(isocode == "USA", "United States of America", country)
    )

  if (!is.null(countries)) {
    df <- dplyr::filter(df, country %in% countries)
  }

  df <- df |>
    dplyr::mutate(
      gdpcap = rgdpna / pop,
      rknacapita = rnna / pop
    )

  if (!is.null(max_year)) {
    df <- dplyr::filter(df, year <= max_year)
  }

  vars_keep <- c(
    "country", "isocode", "year", "pop", "hc", "gdpcap",
    "csh_c", "csh_i", "csh_g", "csh_x", "csh_m",
    "rknacapita", "pl_gdpo"
  )

  df <- dplyr::select(df, dplyr::any_of(vars_keep))

  avg <- df |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      dplyr::across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(country = "Average", isocode = "AVG")

  df <- dplyr::bind_rows(avg, df) |>
    dplyr::mutate(country = forcats::fct_relevel(country, "Average")) |>
    dplyr::arrange(country, year)

  if (is.null(code_map)) {
    codes <- sort(unique(df$isocode))
    code_map <- seq_along(codes)
    names(code_map) <- codes
  }

  df$region <- as.integer(unname(code_map[df$isocode]))
  df
}

