# Utility script to download the pinned PWT input and validate it against the
# versioned checksum registry.
# Usage: Rscript scripts/fetch_data.R

source("R/config/constants.R")
source("R/utils/helpers.R")
source("R/data/downloaders.R")

fetch_pinned_data(force = TRUE)
