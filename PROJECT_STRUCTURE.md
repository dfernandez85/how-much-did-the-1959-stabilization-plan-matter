# Project Structure

```text
repository-root/
  main.R                         # Entry point for the replication package
  README.md                      # Repository overview and execution instructions
  PROJECT_STRUCTURE.md           # This document
  scripts/
    fetch_data.R                 # Forces a clean download of the pinned PWT input
  R/
    config/
      constants.R                # Output paths and checksum keys
      packages.R                 # Required packages and loader
      settings.R                 # Countries, windows, and pipeline parameters
    data/
      downloaders.R              # Pinned data bootstrap, checksum validation, and local loading
      processors.R               # Data transformations and estimation-ready panels
    analysis/
      synthetic_control.R        # Estimation pipeline
      journal_exports.R          # Export of final manuscript and appendix outputs
    utils/
      helpers.R                  # General-purpose helper functions
    visualization/
      plots.R                    # Synthetic-control plotting functions
  data/
    raw/
      checksums.txt              # Required md5 checksum registry for raw inputs
      pwt110.xlsx                # Penn World Table 11.0 (not versioned in git)
    processed/                   # Reserved for derived intermediates
  output/
    sessions/
      <timestamp>/
        plots/
          main/
          appendix/
        tables/
          main/
          appendix/
        run_config.yml
        session_info.txt
  renv/
  renv.lock
```

## Execution flow
1. `main.R` loads configuration and packages; if the pinned raw input is missing, it downloads it automatically and validates it against the versioned checksum registry.
2. `R/analysis/synthetic_control.R` estimates the specifications required for the manuscript and the appendix.
3. `R/analysis/journal_exports.R` exports only the final figures and tables and removes intermediate plot and table trees.
4. Each run stores `session_info.txt` and `run_config.yml` as session-level reproducibility metadata.
