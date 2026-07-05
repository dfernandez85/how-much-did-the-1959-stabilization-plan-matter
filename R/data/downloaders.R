raw_data_targets <- function() {
  list(
    list(url = PWT11_URL, dest = PWT11_LOCAL, key = PWT11_CHECKSUM_KEY)
  )
}

download_raw_target <- function(url, dest) {
  ensure_dir(dirname(dest))
  message(sprintf("Downloading %s -> %s", url, dest))
  utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
}

raw_data_is_ready <- function(targets = raw_data_targets(), checksums_path = CHECKSUM_FILE) {
  if (!file.exists(checksums_path)) {
    return(FALSE)
  }

  all(vapply(targets, function(tgt) {
    if (!file.exists(tgt$dest)) {
      return(FALSE)
    }

    tryCatch({
      validate_file_hash(tgt$dest, tgt$key, checksums_path = checksums_path)
      TRUE
    }, error = function(e) {
      FALSE
    })
  }, logical(1)))
}

fetch_pinned_data <- function(force = FALSE,
                              targets = raw_data_targets(),
                              checksums_path = CHECKSUM_FILE) {
  if (!file.exists(checksums_path)) {
    stop(sprintf(
      "Pinned checksum file not found: %s. The replication package expects versioned checksums and will not regenerate them automatically.",
      checksums_path
    ))
  }

  if (!force && raw_data_is_ready(targets = targets, checksums_path = checksums_path)) {
    message("Pinned raw data already present with valid checksums.")
    return(invisible(targets))
  }

  invisible(lapply(targets, function(tgt) {
    download_raw_target(tgt$url, tgt$dest)
    validate_file_hash(tgt$dest, tgt$key, checksums_path = checksums_path)
  }))
  message("Pinned raw data downloaded and validated against versioned checksums.")

  invisible(targets)
}

ensure_raw_data <- function(targets = raw_data_targets(), checksums_path = CHECKSUM_FILE) {
  if (!raw_data_is_ready(targets = targets, checksums_path = checksums_path)) {
    message("Required raw data missing or invalid. Bootstrapping pinned sources.")
    fetch_pinned_data(force = TRUE, targets = targets, checksums_path = checksums_path)
  }

  invisible(TRUE)
}

load_pwt_data <- function(local_path = PWT11_LOCAL,
                          checksum_key = PWT11_CHECKSUM_KEY) {
  validate_file_hash(local_path, checksum_key)

  read_pwt11 <- function(path) {
    readxl::read_excel(path, sheet = "Data") |>
      dplyr::rename(isocode = countrycode)
  }

  tibble::as_tibble(read_pwt11(local_path))
}
