check_r_version <- function(required = REQUIRED_R_VERSION) {
  current <- paste(getRversion()[1, 1:2], collapse = ".")
  if (current != required) {
    stop(sprintf(
      paste(
        "Version de R no soportada: %s (se requiere R %s.x, version fijada en renv.lock).",
        "R 4.6 cambio su estandar C por defecto a C23, lo que rompe la compilacion de",
        "varios paquetes CRAN necesarios aqui (cli, data.table, rlang, curl, ...).",
        "Instala R %s.x desde https://cran.r-project.org/bin/windows/base/old/ y vuelve a intentarlo."
      ),
      as.character(getRversion()), required, required
    ))
  }
  invisible(TRUE)
}

ensure_dir <- function(path) {
  if (is.null(path) || !nzchar(path)) return(invisible(path))
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

read_checksums <- function(path = CHECKSUM_FILE) {
  if (!file.exists(path)) {
    stop(sprintf("No existe el archivo de checksums: %s", path))
  }

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) == 0) {
    return(setNames(character(0), character(0)))
  }

  parsed <- strsplit(lines, "\\s+", perl = TRUE)
  keys <- vapply(parsed, function(x) x[[1]], character(1))
  vals <- vapply(parsed, function(x) x[[2]], character(1))
  stats::setNames(vals, keys)
}

check_file_md5 <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("No existe el archivo: %s", path))
  }
  unname(tools::md5sum(path))
}

validate_file_hash <- function(path, checksum_key, checksums_path = CHECKSUM_FILE) {
  expected <- read_checksums(checksums_path)[[checksum_key]]
  if (is.null(expected) || !nzchar(expected)) {
    stop(sprintf("No hay checksum registrado para la clave '%s'", checksum_key))
  }

  actual <- check_file_md5(path)
  if (!identical(tolower(actual), tolower(expected))) {
    stop(sprintf(
      "Checksum incorrecto para %s. Esperado: %s. Observado: %s",
      path,
      expected,
      actual
    ))
  }

  invisible(TRUE)
}

save_table_jpg <- function(df, path, title = NULL, max_rows = 40) {
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  if (nrow(df) == 0) {
    warning(sprintf("save_table_jpg skipped for %s: zero rows", path))
    return(invisible(NULL))
  }
  if (nrow(df) > max_rows) {
    df <- df[seq_len(max_rows), , drop = FALSE]
  }

  font_size <- 0.7
  tbl <- gridExtra::tableGrob(
    df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core = list(fg_params = list(cex = font_size)),
      colhead = list(fg_params = list(fontface = "bold", cex = font_size))
    )
  )

  if (!is.null(title)) {
    tbl <- gridExtra::arrangeGrob(
      tbl,
      top = grid::textGrob(title, gp = grid::gpar(fontsize = 12, fontface = "bold"))
    )
  }

  width_cm <- min(60, 3 + 2.8 * ncol(df))
  height_cm <- min(60, 2 + 0.6 * nrow(df))

  ggplot2::ggsave(
    filename = path,
    plot = tbl,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 300
  )
}

