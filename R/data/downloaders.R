raw_data_targets <- function() {
  list(
    list(url = PWT11_URL, dest = PWT11_LOCAL, key = PWT11_CHECKSUM_KEY)
  )
}

# The pinned PWT file is served through a Dataverse redirect to an object store
# whose latency is occasionally well above R's 60-second default timeout. The
# retries and the raised timeout keep a transient stall from aborting a full
# replication run; the checksum validation downstream is unchanged.
download_raw_target <- function(url, dest, attempts = 3L, timeout_seconds = 600L) {
  ensure_dir(dirname(dest))

  old_timeout <- getOption("timeout")
  options(timeout = max(old_timeout, timeout_seconds))
  on.exit(options(timeout = old_timeout), add = TRUE)

  for (attempt in seq_len(attempts)) {
    message(sprintf("Downloading %s -> %s (attempt %d of %d)", url, dest, attempt, attempts))

    status <- tryCatch(
      utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE),
      error = function(e) {
        message(sprintf("Download attempt %d failed: %s", attempt, conditionMessage(e)))
        1L
      }
    )

    if (identical(as.integer(status), 0L) && file.exists(dest) && file.size(dest) > 0) {
      return(invisible(dest))
    }

    if (file.exists(dest)) unlink(dest)
    if (attempt < attempts) Sys.sleep(10 * attempt)
  }

  stop(sprintf(
    "Failed to download %s after %d attempts. The pinned source may be temporarily unavailable; retry later or place the file manually at %s.",
    url, attempts, dest
  ))
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
