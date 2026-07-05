source("R/config/constants.R")
source("R/config/packages.R")
source("R/config/settings.R")
source("R/utils/helpers.R")
source("R/data/downloaders.R")
source("R/data/processors.R")
source("R/visualization/plots.R")
source("R/analysis/journal_exports.R")
source("R/analysis/synthetic_control.R")

check_r_version()
load_required_packages()
ensure_raw_data()

ensure_dir(OUTPUT_SESSIONS)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
session_dir <- file.path(OUTPUT_SESSIONS, timestamp)
ensure_dir(session_dir)

session_info_path <- file.path(session_dir, "session_info.txt")
writeLines(capture.output(utils::sessionInfo()), session_info_path)

cat(sprintf("Starting session at: %s\n", session_dir))


synthetic_outputs <- run_synthetic_control(session_dir)

cat("Analysis completed. Outputs saved in:\n")
cat(sprintf("Plots: %s\n", file.path(session_dir, OUTPUT_PLOTS)))
cat(sprintf("Tables: %s\n", file.path(session_dir, OUTPUT_TABLES)))
