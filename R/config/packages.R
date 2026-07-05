# Packages required across the project
required_packages <- c(
  "rio",
  "lubridate",
  "zoo",
  "dplyr",
  "pracma",
  "quantmod",
  "tidyr",
  "ggplot2",
  "scales",
  "forcats",
  "textclean",
  "readxl",
  "gridExtra",
  "Synth",
  "MSCMT",
  "reactable",
  "parallel",
  "tibble",
  "readr",
  "DEoptim"
)

load_required_packages <- function() {
  missing_pkgs <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
  if (length(missing_pkgs) > 0) {
    stop(sprintf(
      "Faltan paquetes: %s. Ejecuta renv::restore() para instalarlos con versiones fijadas.",
      paste(missing_pkgs, collapse = ", ")
    ))
  }
  lapply(required_packages, library, character.only = TRUE)
  invisible(TRUE)
}
