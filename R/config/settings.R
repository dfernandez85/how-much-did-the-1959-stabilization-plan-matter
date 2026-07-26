
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
# La solucion sin restricciones concentra 0.664 en Nicaragua y no pasa el
# limite de concentracion; el fallback selecciona el pool sin Estados Unidos
# (peso maximo 0.373 en Nicaragua). La sensibilidad de esta seleccion al
# umbral se exporta como Tabla A3 del apendice: el pool elegido no cambia
# en 0.40/0.45/0.50, y la variante con el gate desactivado es la fila
# "Full pool incl. US" de la Tabla 8.
stability_gate_enabled <- TRUE
stability_gate_specs <- c("baseline")
stability_top_weight_max <- 0.45
# Umbrales alternativos evaluados a posteriori sobre los mismos candidatos
# ya estimados (no re-estiman nada): alimentan la tabla de sensibilidad.
stability_top_weight_sensitivity <- c(0.40, 0.45, 0.50)

# Bateria de predictores (Ferman-Pinto-Possebom 2020) para gdpcap/baseline:
# re-estima el baseline sobre el mismo pool variando solo el set de
# predictores (todos los lags del outcome, lags impares, first/mid/last).
# Alimenta la Tabla 7 del manuscrito.
predictor_battery_enabled <- TRUE

# Pools de donantes restringidos para gdpcap/baseline (Tabla 8):
# Europa-only, Europa sin Portugal (donante casi-tratado: EFTA 1960) y
# LatAm-only (America sin Canada/EEUU). Mismos predictores que el baseline.
# Incluye ademas la variante "pool completo sin stability gate" (EEUU dentro),
# que documenta que el gate liga por concentracion (Nicaragua), no por EEUU.
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

# Outcomes activos del pipeline. gdpcap es el outcome principal (tablas y
# figuras del manuscrito); rknacapita alimenta el Apendice B (Figuras B1-B2,
# Tablas B1-B2) y hc el Apendice C (Figuras C1-C2, Tablas C1-C2). Verificado
# contra el manuscrito v4: AMBOS outcomes auxiliares se usan (Seccion de
# robustez y conclusiones: contraste "capital deepening si / break en capital
# humano no" como argumento de eficiencia asignativa). No desactivar sin
# revisar el texto; hc es ademas PREDICTOR del baseline en todo caso.
enabled_outcomes <- c("gdpcap", "rknacapita", "hc")

# Presupuesto del optimizador externo (DEoptim) para la busqueda de la matriz V.
# MSCMT fija por defecto control=list(reltol=1e-14, steptol=500) pero deja
# itermax en el default de DEoptim (200), que es lo que ata la busqueda. Al
# ampliar el pool (p. ej. incluyendo Estados Unidos, isocode USA) el paisaje de
# optimizacion se vuelve mas rugoso y 200 iteraciones no bastan: el ajuste pre
# se dispara y deja de ser reproducible entre semillas. Subimos itermax y
# steptol para recuperar convergencia estable; se inyectan via outer.opar en
# cada llamada a MSCMT::mscmt. NA/<=0 -> se usa el default de MSCMT.
deoptim_itermax <- 2000
deoptim_steptol <- 1000
