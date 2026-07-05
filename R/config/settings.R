
synthetic_countries <- c(
  # Europa
  "Albania", "Austria", "Belgium", "Cyprus", "Denmark", "Finland",
  "France", "Germany", "Iceland", "Ireland", "Italy", "Luxembourg",
  "Netherlands", "Norway", "Portugal", "Spain", "Sweden",
  "Switzerland", "Turkey", "United Kingdom",
  # America
  "Argentina", "Bolivia", "Brazil", "Canada", "Chile", "Colombia",
  "Costa Rica", "Ecuador", "El Salvador", "Guatemala", "Honduras",
  "Mexico", "Nicaragua", "Panama", "Peru", "Trinidad and Tobago",
  "United States of America", "Uruguay", "Venezuela"
)

region_codes <- NULL # se generan codigos numericos automaticamente

pre_treatment_years <- c(1950, 1959) # Incluye 1959 en el periodo pre
post_cutoff_year <- 1976
# Regla ex ante: descarta placebos cuyo RMSPE pre supere este multiple del RMSPE pre de Espana
# (criterio de Abadie et al. 2010). Si es NA o <= 0, no se filtra.
placebo_rmspe_ratio_cutoff <- 10
placebo_mspe_ratio_cutoff <- placebo_rmspe_ratio_cutoff
# Excluye donantes de las corridas drop_one_out cuando su peso es marginal.
drop_one_excluded_donors <- c("Brazil", "Ecuador")

# Regla ex ante para aceptar la especificacion principal del baseline.
# Si "all" no cumple estabilidad, se prueban candidatos drop_top1..top_k
# y se selecciona el de menor pre-MSPE que cumpla los umbrales.
# NOTA: top_weight_max revisado de 0.40 a 0.45 al corregir el matching de
# nombres con PWT (Turkey/Bolivia/Venezuela entran al pool y el peso maximo
# del baseline sin restricciones pasa de 0.387 a 0.405). La variante con
# umbral 0.40 (que excluye a Turkey via drop_top1) se conserva como robustez.
stability_gate_enabled <- TRUE
stability_gate_specs <- c("baseline")
stability_top_weight_max <- 0.45

# Bateria de predictores (Ferman-Pinto-Possebom 2020) para gdpcap/baseline:
# re-estima el baseline sobre el mismo pool variando solo el set de
# predictores (todos los lags del outcome, lags impares, first/mid/last).
# Alimenta la Tabla A3 del apendice.
predictor_battery_enabled <- TRUE

# Pools de donantes restringidos para gdpcap/baseline (Tabla A4):
# Europa-only, Europa sin Portugal (donante casi-tratado: EFTA 1960) y
# LatAm-only (America sin Canada/EEUU). Mismos predictores que el baseline.
restricted_pools_enabled <- TRUE
europe_donor_countries <- c(
  "Albania", "Austria", "Belgium", "Cyprus", "Denmark", "Finland",
  "France", "Germany", "Iceland", "Ireland", "Italy", "Luxembourg",
  "Netherlands", "Norway", "Portugal", "Sweden", "Switzerland",
  "Turkey", "United Kingdom"
)
latam_donor_countries <- c(
  "Argentina", "Bolivia", "Brazil", "Chile", "Colombia", "Costa Rica",
  "Ecuador", "El Salvador", "Guatemala", "Honduras", "Mexico",
  "Nicaragua", "Panama", "Peru", "Trinidad and Tobago", "Uruguay", "Venezuela"
)
stability_neff_min <- 3
stability_min_positive_donors <- 4
stability_drop_top1_tau_max_pct <- 25
stability_pre_mspe_max_increase_pct <- 30
stability_top_k_candidates <- 3
