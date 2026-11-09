SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## 00_parameter_grids.R -- the five Latin-hypercube grids of kinetic parameters,
## 2,000 vectors each, written to ./data/tau_grid_v6_<n>.rds.
##
## NO RANDOM SEED WAS SET when the stored grids were drawn, so re-running this
## script draws DIFFERENT grids. The five stored grids are the record of the
## parameter sets every result rests on; they are kept in ./data and this script
## refuses to overwrite them unless OVERWRITE=1. SEED=<integer> makes a new draw
## reproducible. Only grids 1-3 were swept.
## RUN: OVERWRITE=1 SEED=1 Rscript pipeline/00_parameter_grids.R
## ============================================================================
suppressMessages(library(data.table))
.SEED = Sys.getenv("SEED", "")
if (any(file.exists(sprintf("./data/tau_grid_v6_%d.rds", 1:5))) && Sys.getenv("OVERWRITE") != "1")
  stop("./data/tau_grid_v6_*.rds exist and are the stored record of the swept parameter sets -- set OVERWRITE=1 to replace them")
if (nzchar(.SEED)) set.seed(as.integer(.SEED))

# =========================================================================
# OVERVIEW
# =========================================================================
# v4 denormalized: ODE now uses raw state variables (M0, M1, Tr, Ea)
# instead of normalized fractions (f_1 = M1/K_M, f_r = Tr/K_T, etc.).
# Carrying capacities: K_M = 100, K_T = 100, K_R = 1.
#
# Parameters whose equation terms changed from fraction to raw variable
# have their sampling ranges divided by the relevant carrying capacity
# to preserve equivalent dynamical regimes.
#
# Rescaled (/ K_M=100):  tau_t_a, alpha_0, alpha_1_e
# Rescaled (/ K_T=100):  sigma_S
# Rescaled (x K_M=100):  tau_eff  -- see note below, this one is CORRECT
# Unchanged:             everything else (see inline comments)
#
# ---------------------------------------------------------------------
# CORRECTION (this revision): R IS NOT O(1). The previous header claimed
# "R in [0, K_R=1] in both versions". It is not, and that error propagated
# through everything.
#
#   dR = beta_ros * M1 - mu_R * R   =>   R* = beta_ros * M1 / mu_R
#
# In the normalized model M1 was f_1 in [0,1], giving R* <~ 10. In v4 M1 is
# raw (0..100), so R* <~ 1000. beta_ros multiplies M1 and therefore needed
# the same / K_M rescaling as alpha_0 / alpha_1_e / tau_t_a -- but it is
# swept as an explicit grid rather than sampled here, so it fell outside
# the rescaling pass and was never corrected.
#
# Measured consequence in the v4 sweep: attr_R had median 20, max 959, and
# the epithelial fixed point
#       Eh* = K_E * kappa_E / (kappa_E + tau_r*R + tau_P_B*PLP/K_LP_P)
# collapsed to a median of 1.97 -- i.e. 2% of capacity. 95% of all
# attractors sat below Eh = 26. The model never healed anywhere in 20000
# parameter sets; "resolved" vs "chronic" meant 2% vs 0.5% epithelium.
#
# THE FIX IS NOT IN THIS FILE. beta_ros is swept in the sweep scripts:
#       BETA_ROS_VALS = seq(0.0015, 0.0095, by = 0.001)   # v6; see 01_sweep.R
# With that, tau_r*R becomes comparable to kappa_E and alpha_1_r*R becomes
# comparable to r_L_P -- both balances straddle their transition instead of
# sitting ~100x on one side. The ranges below were evidently chosen for that
# regime, so NONE of them need widening to compensate.
#
# WHY tau_eff IS x100 AND NOT /100 (it looks wrong, it is not):
#   fb = tau_eff * g * Ea * M1, and g was rescaled DOWN by K_M
#   (a0 = 0.0045 = 0.45/100) while M1 was scaled UP by K_M. Those cancel in
#   the product g*M1, so recovering fb_denorm = K_M * fb_norm requires
#   tau_eff = K_M * tau_eff_norm. Verified against lib/ODE_SYSTEM_v6.R.
# ---------------------------------------------------------------------
#
# SAMPLING SCALE (this revision): rate constants are now drawn LOG-uniformly.
# Linear LHS on a range spanning decades never visits the slow end -- on
# [0.001, 0.5] only 0.2% of draws fall below 0.002, so the slow-rate regime
# was unsampled. Several parameters span 100-500x (tau_D, tau_S at 500x;
# tau_01/tau_10/tau_1_d/tau_eff/tau_t_a/tau_t_d at 190x; gamma_B_* at 100x).
# Log-uniform is the standard non-informative choice for a scale parameter
# and gives even coverage per decade, which is what actually matters when
# the quantity of interest is a RATIO of rates. The rule applied here:
#     log   = every strictly-positive rate constant
#     linear = bounded dimensionless weights (w_E_over_M) and fixed knobs
# This changes the study population -- expect a different bistability yield.
#
# UPDATED: alpha_1_p removed. The model now uses the unified-Hill
# efferocytosis switch g(L) = a_floor + (a0-a_floor)*Kg^n_g/(Kg^n_g+L^n_g)
# in place of the old flat alpha_1_p clearance rate (see lib/ODE_SYSTEM_v6.R).
# a0, a_floor, Kg, n_g, and beta_ros are NOT sampled here via LHS -- they
# are swept explicitly as a grid in the sweep scripts, applied on top of
# every LHS-sampled parameter vector in this file. tau_r and mu_R
# (ROS tissue damage / decay) ARE still sampled here and are now actually
# used by the RHS (they previously were sampled but unused).
# =========================================================================

for (ind in 1:5){
  rm(list = setdiff(ls(), c("ind", ".SEED")))
  # 2000 per batch (x5 batches = 10000 tau rows). The v4 sweep used 20000 per
  # batch, which was affordable only with the 6-rung fan; at 101 Eh0 rungs x 4
  # PLP0 levels a tau row costs hours, so 20000 would never finish.
  n_keep = 2000

  # =========================================================================
  # Latin Hypercube Sampling helper
  # =========================================================================
  lhs_sample = function(n, k) {
    mat = matrix(0, nrow = n, ncol = k)
    for (j in 1:k) {
      mat[, j] = (sample(n) - runif(n)) / n
    }
    mat
  }

  # =========================================================================
  # Parameter definitions: name, min, max, beta shape1, shape2, log
  # All parameters sampled freely -- nothing back-solved.
  # log = TRUE  -> uniform in log10 between min and max (rate constants)
  # log = FALSE -> uniform on the linear scale (bounded weights, fixed knobs)
  # =========================================================================
  param_specs = list(
    # Epithelial dynamics — unchanged (Eh, R, PLP/K_LP_P appear identically)
    kappa_E     = list(min = 0.05,   max = 0.15,   a = 1, b = 1, log = TRUE),  # healing rate
    tau_r       = list(min = 0.05,   max = 0.50,   a = 1, b = 1, log = TRUE),  # ROS tissue damage rate
    tau_P_B     = list(min = 0.05,   max = 0.50,   a = 1, b = 1, log = TRUE),  # pathogen-driven primary necrosis
    mu_Ea       = list(min = 0.05,   max = 0.95,   a = 1, b = 1, log = TRUE),  # secondary necrosis rate (Ea -> lysis + DAMPs)

    # ROS decay. NOTE: R* = beta_ros*M1/mu_R, so mu_R sets the ROS scale
    # jointly with the swept beta_ros -- see the header correction.
    mu_R        = list(min = 0.05,   max = 0.50,   a = 1, b = 1, log = TRUE),  # ROS natural decay rate

    # Macrophage dynamics — tau_01, tau_10, tau_1_d unchanged (D*M0, S*M1, M1 same in both)
    tau_01      = list(min = 0.005,  max = 0.95,   a = 1, b = 1, log = TRUE),  # M0->M1 activation
    tau_10      = list(min = 0.005,  max = 0.95,   a = 1, b = 1, log = TRUE),  # Treg-driven M1 suppression
    tau_1_d     = list(min = 0.005,  max = 0.95,   a = 1, b = 1, log = TRUE),  # M1 natural turnover
    tau_eff     = list(min = 100*0.005, max = 100*0.95, a = 1, b = 1, log = TRUE),  # efferocytosis (x100: see header)

    # Treg dynamics
    tau_t_a     = list(min = 0.005/100, max = 0.95/100, a = 1, b = 1, log = TRUE),  # Treg activation rate (/100)
    tau_t_d     = list(min = 0.005,  max = 0.95,   a = 1, b = 1, log = TRUE),  # Treg deactivation rate

    # Clearance. alpha_0 and alpha_1_e are ~1000x weaker than alpha_1_r*R and
    # contribute <0.1% of pathogen clearance each; that is a property of the
    # model, not of these bounds, and is NOT fixed by changing them.
    alpha_0     = list(min = 0.005/100, max = 0.01/100, a = 1, b = 1, log = TRUE),  # M0 basal clearance of microbes
    alpha_1_e   = list(min = 0.01/100,  max = 0.05/100, a = 1, b = 1, log = TRUE),  # M1 direct clearance of microbes
    alpha_1_r   = list(min = 0.05,   max = 0.95,   a = 1, b = 1, log = TRUE),  # ROS-augmented clearance
    alpha_AMP   = list(min = 0.01,   max = 0.10,   a = 1, b = 1, log = TRUE),  # AMP killing rate

    # Pathogen / commensal dynamics — unchanged (f_Ei = 1-Eh/K_E identical)
    gamma_B_P   = list(min = 0.001,  max = 0.10,   a = 1, b = 1, log = TRUE),  # pathogen translocation
    gamma_B_C   = list(min = 0.001,  max = 0.10,   a = 1, b = 1, log = TRUE),  # commensal translocation

    # Danger and suppression signals
    sigma_D     = list(min = 0.05,     max = 0.95,     a = 1, b = 1, log = TRUE),  # danger signal production
    sigma_S     = list(min = 0.05/100, max = 0.95/100, a = 1, b = 1, log = TRUE),  # suppression signal production (/100)
    tau_D       = list(min = 0.001,    max = 0.50,     a = 1, b = 1, log = TRUE),  # D decay rate
    tau_S       = list(min = 0.001,    max = 0.50,     a = 1, b = 1, log = TRUE),  # S decay rate

    # Bounded dimensionless weight -- LINEAR. Log would pile draws at 0.05 and
    # misrepresent a 0-1 mixing weight, which has no scale to be uniform in.
    w_E_over_M  = list(min = 0.05,     max = 0.95,     a = 1, b = 1, log = FALSE),

    # Knobs — fixed at zero, linear (log is undefined at 0)
    phi_instability = list(min = 0,  max = 0,  a = 1, b = 1, log = FALSE),   # Robustness of Treg phenotype under inflammation
    treg_eff        = list(min = 0,  max = 0,  a = 1, b = 1, log = FALSE)    # Treg efficiency
  )

  # =========================================================================
  # Check for placeholder bounds before proceeding
  # =========================================================================
  na_params = names(Filter(function(s) is.na(s$min) || is.na(s$max), param_specs))
  if (length(na_params) > 0) {
    stop(sprintf("Fill in bounds for: %s", paste(na_params, collapse = ", ")))
  }
  # log sampling needs a strictly positive range
  bad_log = names(Filter(function(s) isTRUE(s$log) && (s$min <= 0 || s$max <= 0), param_specs))
  if (length(bad_log) > 0) {
    stop(sprintf("log = TRUE requires min > 0 and max > 0: %s", paste(bad_log, collapse = ", ")))
  }

  n_params    = length(param_specs)
  param_names = names(param_specs)

  cat(sprintf("Drawing %d samples across %d parameters\n", n_keep, n_params))

  # =========================================================================
  # LHS on all parameters
  # =========================================================================
  U        = lhs_sample(n_keep, n_params)
  tau_grid = data.table(matrix(0, nrow = n_keep, ncol = n_params))
  setnames(tau_grid, param_names)

  for (j in seq_along(param_specs)) {
    spec          = param_specs[[j]]
    beta_vals     = qbeta(U[, j], spec$a, spec$b)
    tau_grid[[j]] = if (isTRUE(spec$log)) {
      # uniform per decade between min and max
      exp(log(spec$min) + beta_vals * (log(spec$max) - log(spec$min)))
    } else {
      spec$min + beta_vals * (spec$max - spec$min)
    }
  }

  # =========================================================================
  # Derived quantities
  # =========================================================================
  tau_grid[, row_ind := .I]

  ## ---- THE GLOBAL KEY, ASSIGNED HERE AND NOWHERE ELSE ---------------------
  ## row_ind restarts at 1 in every grid, so it does not identify a parameter
  ## set on its own. Until now `batch` and `tau_id` were reconstructed in
  ## 00_read_first.R, and three downstream scripts inverted the result by hand
  ## as tau_id %/% TAU_STRIDE against a stride that was nominally a variable.
  ## Stamping both into the grid at the point the rows are created makes the key
  ## a property of the data: the sweep carries the columns through, nothing has
  ## to rebuild them, and a merge that drops `batch` fails loudly instead of
  ## silently returning another batch's parameters.
  ##
  ## TAU_STRIDE must exceed the largest row count any grid will ever have. It is
  ## asserted rather than assumed, because a grid that overflowed it would
  ## produce tau_id collisions across batches -- the exact failure the key
  ## exists to prevent.
  TAU_STRIDE = as.integer(Sys.getenv("TAU_STRIDE", "10000"))
  stopifnot(nrow(tau_grid) < TAU_STRIDE)
  tau_grid[, batch      := as.integer(ind)]
  tau_grid[, tau_stride := TAU_STRIDE]
  tau_grid[, tau_id     := TAU_STRIDE * as.integer(ind) + row_ind]

  # -------------------------------------------------------------------------
  # DERIVED SCALES. The quantities that decide the dynamics are RATIOS of
  # rates, and those ratios depend on several sampled parameters at once --
  # so you cannot choose any single parameter's range to make a ratio come
  # out right. Record the ratios instead, verify their coverage BEFORE
  # spending compute, and stratify the analysis by them rather than by the
  # raw parameter. They need no integration.
  #
  # beta_ros is swept, not sampled, so everything ROS-related is given PER
  # UNIT beta_ros; multiply by the beta_ros value in use. M1 is taken at its
  # ceiling K_M, so these are UPPER BOUNDS (real M1 is often 5-70).
  # -------------------------------------------------------------------------
  K_M = 100; K_E = 100
  tau_grid[, R_per_beta       := K_M / mu_R]                    # R*  = beta_ros * this
  tau_grid[, ros_damage_ratio := tau_r     * R_per_beta / kappa_E]  # tau_r*R / kappa_E, per unit beta_ros
  tau_grid[, ros_clear_per_beta := alpha_1_r * R_per_beta]      # alpha_1_r*R, per unit beta_ros
  tau_grid[, phag_clear       := (alpha_0 + alpha_1_e) * K_M]   # alpha_0*M0 + alpha_1_e*M1, at M1 = K_M
  # Eh at the fixed point, ROS term only, per unit beta_ros:
  #   Eh* = K_E / (1 + ros_damage_ratio * beta_ros)
  # Sanity: with the OLD beta_ros = 0.26 (grid midpoint) this predicts the
  # observed median Eh of ~2; with the corrected 0.0026 it predicts ~65.
  ## Reference beta_ros for this diagnostic = the midpoint of the v6 grid,
  ## seq(0.0015, 0.0095, by = 0.001), i.e. 0.0055. v5 used 0.0026, the midpoint of
  ## its own narrower grid, so the two numbers are not comparable across versions.
  tau_grid[, Eh_star_at := K_E / (1 + ros_damage_ratio * 0.0055)]

  # =========================================================================
  # Summary
  # =========================================================================
  cat("\nParameter ranges:\n")
  for (nm in param_names) {
    cat(sprintf("  %-20s: [%.6g, %.6g]  median=%.6g  %s\n",
                nm, min(tau_grid[[nm]]), max(tau_grid[[nm]]),
                median(tau_grid[[nm]]),
                if (isTRUE(param_specs[[nm]]$log)) "(log)" else "(linear)"))
  }

  cat("\nDerived scales (per unit beta_ros where noted):\n")
  for (nm in c("R_per_beta","ros_damage_ratio","ros_clear_per_beta","phag_clear","Eh_star_at")) {
    q = quantile(tau_grid[[nm]], c(0, .25, .5, .75, 1), na.rm = TRUE)
    cat(sprintf("  %-20s: %s\n", nm, paste(sprintf("%.4g", q), collapse = "  ")))
  }
  cat("\n  Eh_star_at is the predicted epithelial fixed point at beta_ros = 0.0055\n")
  cat("  (midpoint of the CORRECTED grid). If its median is not well above 1,\n")
  cat("  the ROS term still dominates healing and the beta_ros rescaling has\n")
  cat("  not been applied in the sweep scripts.\n")

  saveRDS(tau_grid, paste0('./data/tau_grid_v6_', ind, '.rds'))
  cat(sprintf("\nSaved to ./data/tau_grid_v6_%d.rds\n", ind))

}
