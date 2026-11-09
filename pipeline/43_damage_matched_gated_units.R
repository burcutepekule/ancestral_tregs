## ============================================================================
## pipeline/43_damage_matched_gated_units.R -- the units of the damage-matched
## experiment drawn from the WHOLE Tregs-OFF bistable set, with no Treg-effect
## criterion, among the units in which Tregs OFF does not eliminate the pathogen.
## ----------------------------------------------------------------------------
## The stored experiment (pipeline/40_damage_matched.R) draws its cells from the
## substrate of pipeline/05_substrate_cells.R: cells in which Tregs help at one or
## more beta_r levels. Here the only gate is the one of pipeline/03_select.R on the
## Tregs-OFF system (bistable; ./data/selected_big.rds), and a unit is
##   one gated pathogenic combination (kinetic vector, K_g, r_P, beta_r)  x  E_h(0)
## so every unit sits at a beta_r at which its own combination is bistable without
## Tregs. beta_r is therefore one of the nine levels of the main sweep
## (0.0015-0.0095): the levels 0.0125-0.05 were swept only for the substrate cells
## and using them would bring the Treg-effect criterion back in.
##
## PERSIST=1 (default): only the units in which Tregs OFF is predicted NOT to
## eliminate the pathogen -- the attractor reached from E_h(0) (the low-burden one
## if E_h(0) >= the threshold Eh0_off, else the high-burden one) has P_LP > P_TOL.
## Checked against the stored experiment at the two shared levels: 1,735 of 1,737
## units predicted correctly. The experiment applies the exact test again on its own
## Tregs-OFF run (SKIP_OFF_ELIM=1).
## WITH THIS GATE THE SAMPLE CAN SHOW ONLY GAINS OF ELIMINATION, never losses: a
## unit in which Tregs OFF eliminates is not in it.
##
## N_UNIT units are drawn uniformly at random (SEED) and written IN THE ORDER DRAWN:
## an array task runs its units in file order, so whatever part of the file has run
## when a task stops at the walltime is itself a random sample. treg_effect (dead band
## EFFECT_TOL, as 05_substrate_cells.R) is recorded per unit and NOT used.
##
## in : ./data/selected_big.rds
## out: ./data/fig04_gated_units.rds  columns batch, row_ind, tau_id, Kg, r_L_P,
##      b_off, Eh0 (what 40_damage_matched.R reads) + Eh0_off, P_pred, treg_effect
## RUN: Rscript pipeline/43_damage_matched_gated_units.R
## ============================================================================
suppressMessages(library(data.table))
SRC = Sys.getenv("SRC", "."); setwd(SRC)
N_UNIT     = as.integer(Sys.getenv("N_UNIT", "75000"))
SEED       = as.integer(Sys.getenv("SEED", "1"))
EH0        = as.numeric(strsplit(Sys.getenv("EH0", "20,50,80"), ",")[[1]])
PERSIST    = Sys.getenv("PERSIST", "1") == "1"
P_TOL      = as.numeric(Sys.getenv("P_TOL", "1e-3"))
EFFECT_TOL = as.numeric(Sys.getenv("EFFECT_TOL", "1.0"))
OUT        = Sys.getenv("OUT", "./data/fig04_gated_units.rds")

S = as.data.table(readRDS("./data/selected_big.rds"))[r_L_P > 0 & is.finite(Eh0_off)]
cmpf = function(Eh, PLP) sqrt((100 - Eh)^2 + PLP^2)
S[, d_lo := cmpf(lo_on_Eh, lo_on_PLP) - cmpf(lo_off_Eh, lo_off_PLP)]
S[, treg_effect := fifelse(!is.finite(d_lo), "unknown", fifelse(d_lo < -EFFECT_TOL, "help", fifelse(d_lo > EFFECT_TOL, "harm", "neutral")))]
U = CJ(k = seq_len(nrow(S)), Eh0 = EH0)
U = cbind(S[U$k, .(batch = as.integer(batch), row_ind = as.integer(row_ind), tau_id = as.integer(tau_id), Kg = as.integer(Kg),
                   r_L_P = as.numeric(r_L_P), b_off = as.numeric(beta_ros), Eh0_off, lo_off_PLP, hi_off_PLP, treg_effect)], Eh0 = U$Eh0)
U[, P_pred := fifelse(Eh0 >= Eh0_off, lo_off_PLP, hi_off_PLP)][, c("lo_off_PLP", "hi_off_PLP") := NULL]
cat(sprintf("[gated units] %d gated pathogenic combinations in %d cells x %d starts = %d units; Tregs OFF predicted persistent in %d (%.0f%%)\n",
            nrow(S), uniqueN(S[, .(tau_id, Kg, r_L_P)]), length(EH0), nrow(U), U[P_pred > P_TOL, .N], 100 * U[, mean(P_pred > P_TOL)]))
if (PERSIST) U = U[P_pred > P_TOL]
setorder(U, tau_id, Kg, r_L_P, b_off, Eh0)          # a fixed order: the draw depends on it
set.seed(SEED)
U = U[sample.int(.N, min(N_UNIT, .N))]            # kept in the order drawn (see header)
setcolorder(U, c("batch", "row_ind", "tau_id", "Kg", "r_L_P", "b_off", "Eh0"))
saveRDS(U, OUT)
cat(sprintf("[gated units] wrote %s: %d units in %d cells, %d tau rows\n", OUT, nrow(U), uniqueN(U[, .(tau_id, Kg, r_L_P)]), uniqueN(U$tau_id)))
print(dcast(U[, .N, by = .(beta_r = signif(b_off, 3), Eh0)], beta_r ~ Eh0, value.var = "N", fill = 0))
print(U[, .N, by = treg_effect][order(-N)])
