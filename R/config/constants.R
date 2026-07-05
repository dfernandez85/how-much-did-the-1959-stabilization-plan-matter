# Pinned in renv.lock. R 4.6 switched its default C dialect to C23, which
# breaks compilation of several required CRAN packages (cli, data.table,
# rlang, curl, ...), so the major.minor version must match exactly.
REQUIRED_R_VERSION <- "4.5"

# PWT 11.0 hosted in dataverse.nl (public file id 554105)
PWT11_URL <- "https://dataverse.nl/api/access/datafile/554105"
PWT11_LOCAL <- file.path("data", "raw", "pwt110.xlsx")

PWT11_CHECKSUM_KEY <- "pwt110"
CHECKSUM_FILE <- file.path("data", "raw", "checksums.txt")

# Output folders
OUTPUT_BASE   <- "output"
OUTPUT_SESSIONS <- file.path(OUTPUT_BASE, "sessions")
OUTPUT_PLOTS  <- "plots"
OUTPUT_TABLES <- "tables"
