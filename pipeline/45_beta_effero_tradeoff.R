## ============================================================================
## pipeline/45_beta_effero_tradeoff.R -- what raising beta_r does to a Tregs-OFF system,
## read from the stored sweep: the chance that efferocytosis ends overwhelmed and the
## chance that the pathogen is cleared, for every parameter set at all nine beta_r levels.
## ----------------------------------------------------------------------------
## NOTHING IS INTEGRATED HERE. The sweep (pipeline/01_sweep.R, 02_assemble.R) already
## stores, per combination (tau_id, K_g, r_P, beta_r), the two Tregs-OFF attractors in
## full (lo_off_*, hi_off_*) and the initial-condition threshold Eh0_off that separates
## their basins on the E_h(0) line. The attractor reached from E_h(0) is the low-burden
## one if E_h(0) >= Eh0_off and the high-burden one otherwise -- the rule the gated unit
## sample uses (pipeline/43_damage_matched_gated_units.R).
##
## THE RULE IS ONLY VALID WHERE THE THRESHOLD IS ONE CROSSING OF THE LINE, so rows are
## kept only if the Tregs-OFF system is bistable (n_field_stable_off >= 2), its threshold
## direction is "increasing", Eh0_off is finite, the fixed-point map is trusted, and both
## attractors carry finite states. A PARAMETER SET ENTERS ONLY IF ALL NINE beta_r LEVELS
## SURVIVE, so every comparison below is within one parameter set, beta_r being the only
## thing that differs.
##
## Per row: L = E_a + E_n at each attractor, g(L)/a_0 there, "overwhelmed" = L > K_g, and
## "cleared" = P_LP < 1e-3. Chances are taken over the injury line, E_h(0) uniform on
## [0, K_E]: the high-burden attractor takes a share Eh0_off / K_E of it.
##
## in : ./data/fp_big_batch_*.rds
## out: ./data/fig_beta_effero.rds  one row per (parameter set, beta_r)
## RUN: Rscript pipeline/45_beta_effero_tradeoff.R
## ============================================================================
suppressMessages(library(data.table))
SRC = Sys.getenv("SRC", "."); setwd(SRC)
EFF   = list(a0 = 0.0045, a_floor = 0.0002, n_g = 4)      # EFFERO of lib/skeleton_defs.R
KE    = 100                                               # K_E of STRUCTURAL
P_CLR = as.numeric(Sys.getenv("P_CLR", "1e-3"))
OUT   = Sys.getenv("OUT", "./data/fig_beta_effero.rds")

B = rbindlist(lapply(sort(Sys.glob("./data/fp_big_batch_*.rds")), function(f) {
  d = as.data.table(readRDS(f))
  d[r_L_P > 0 & n_field_stable_off >= 2 & dir_off == "increasing" &
    is.finite(Eh0_off) & map_trusted_off == TRUE,
    .(tau_id, Kg, r_L_P, beta = beta_ros, e_star = Eh0_off,
      L_lo = lo_off_Ea + lo_off_En, P_lo = lo_off_PLP, Eh_lo = lo_off_Eh,
      L_hi = hi_off_Ea + hi_off_En, P_hi = hi_off_PLP, Eh_hi = hi_off_Eh)] }))
cat(sprintf("[beta-effero] rows passing the attractor rule: %d in %d parameter sets\n",
            nrow(B), uniqueN(B[, .(tau_id, Kg, r_L_P)])))
B = B[is.finite(L_lo) & is.finite(L_hi) & is.finite(P_lo) & is.finite(P_hi)]
full = B[, .N, by = .(tau_id, Kg, r_L_P)][N == 9, .(tau_id, Kg, r_L_P)]
B = merge(B, full, by = c("tau_id", "Kg", "r_L_P"))
cat(sprintf("[beta-effero] parameter sets with all nine beta_r levels: %d (%d rows)\n", nrow(full), nrow(B)))

gf = function(L, Kg) (EFF$a_floor + (EFF$a0 - EFF$a_floor) * Kg^EFF$n_g / (Kg^EFF$n_g + L^EFF$n_g)) / EFF$a0
B[, `:=`(g_lo = gf(L_lo, Kg), g_hi = gf(L_hi, Kg),
         ov_lo = L_lo > Kg, ov_hi = L_hi > Kg,
         cl_lo = P_lo < P_CLR, cl_hi = P_hi < P_CLR,
         w_hi = pmin(pmax(e_star, 0), KE) / KE)]            # share of the injury line that falls to the high-burden attractor
B[, `:=`(p_over = w_hi * ov_hi + (1 - w_hi) * ov_lo,
         p_clear = w_hi * cl_hi + (1 - w_hi) * cl_lo,
         g_exp = w_hi * g_hi + (1 - w_hi) * g_lo,
         Eh_exp = w_hi * Eh_hi + (1 - w_hi) * Eh_lo)]
saveRDS(B[], OUT)
cat(sprintf("[beta-effero] wrote %s\n", OUT))
print(B[, .(sets = .N, overwhelmed_pct = round(100 * mean(p_over), 1), cleared_pct = round(100 * mean(p_clear), 1),
            g_expected = round(mean(g_exp), 3)), keyby = beta])
M = B[order(beta), .(ov = all(diff(p_over) >= -1e-12), cl = all(diff(p_clear) >= -1e-12)), by = .(tau_id, Kg, r_L_P)]
cat(sprintf("[beta-effero] non-decreasing in beta_r within the parameter set: overwhelm %.1f%%, clearance %.1f%% of sets\n",
            100 * mean(M$ov), 100 * mean(M$cl)))
