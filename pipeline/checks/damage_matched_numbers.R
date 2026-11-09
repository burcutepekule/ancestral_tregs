## ============================================================================
## pipeline/checks/damage_matched_numbers.R -- every number of the Methods section
## of Fig. 4 (docs/figure_methods.pdf) and of the
## captions of Fig. 4 and Extended Data Fig. 5, recomputed from ./data/.
## Sections are numbered; the text cites them as [NSk].
## in : ./data/fig04_gated_units.rds, ./data/fig04_damage_matched_gated.rds,
##      ./data/fig04b_v3_trajectories.rds, ./data/aux_gated_fourth_arm.rds,
##      ./data/selected_big.rds
## RUN: Rscript pipeline/checks/damage_matched_numbers.R
## ============================================================================
suppressMessages(library(data.table))
SRC = Sys.getenv("SRC", "."); setwd(SRC)
source("./lib/fig_damage_matched.R")
q = function(x, p = c(.25, .5, .75)) paste(signif(quantile(x, p, na.rm = TRUE), 3), collapse = " / ")
pct = function(x) sprintf("%.1f%%", 100 * mean(x))
sec = function(k, t) cat(sprintf("\n== %d. %s ==\n", k, t))

## ---- 1. the population and the sample ------------------------------------------------
sec(1, "population and sample (pipeline/43_damage_matched_gated_units.R)")
S = as.data.table(readRDS("./data/selected_big.rds"))[r_L_P > 0 & is.finite(Eh0_off)]
cat("gated pathogenic combinations:", nrow(S), " cells:", uniqueN(S[, .(tau_id, Kg, r_L_P)]), " tau rows:", uniqueN(S$tau_id), "\n")
cat("beta_r levels:", paste(sort(unique(signif(S$beta_ros, 3))), collapse = " "), "\n")
U0 = CJ(k = seq_len(nrow(S)), Eh0 = c(20, 50, 80))
U0 = cbind(S[U0$k, .(tau_id, Kg, r_L_P, Eh0_off, lo_off_PLP, hi_off_PLP)], Eh0 = U0$Eh0)
U0[, P_pred := fifelse(Eh0 >= Eh0_off, lo_off_PLP, hi_off_PLP)]
cat("units (x 3 starts):", nrow(U0), " predicted persistent without Tregs (P_LP > 1e-3 at the attractor reached):", U0[P_pred > 1e-3, .N],
    sprintf("(%.1f%%)", 100 * mean(U0$P_pred > 1e-3)), " in", uniqueN(U0[P_pred > 1e-3, .(tau_id, Kg, r_L_P)]), "cells\n")
U = as.data.table(readRDS("./data/fig04_gated_units.rds")); U[, i := .I]
cat("drawn:", nrow(U), "units in", uniqueN(U[, .(tau_id, Kg, r_L_P)]), "cells,", uniqueN(U$tau_id), "tau rows; by Treg effect at own beta_r:\n"); print(U[, .N, by = treg_effect][order(-N)])
print(dcast(U[, .N, by = .(beta_r = signif(b_off, 3), Eh0)], beta_r ~ Eh0, value.var = "N"))
print(U[, .N, by = r_L_P][order(r_L_P)]); print(U[, .N, by = Kg][order(Kg)])

## ---- 2. the run and the search ----------------------------------------------------------
sec(2, "the experiment (pipeline/40_damage_matched.R; pipeline/cluster/damage_matched_gated.sh)")
R = as.data.table(readRDS("./data/fig04_damage_matched_gated.rds"))[match == "conv"]
R = merge(R, U[, .(i, treg_effect, P_pred)], by = "i")
cat("units run:", nrow(R), "of", nrow(U), "drawn (tasks stopped at the walltime; the file order is random)\n")
cat("Tregs OFF eliminated after all (status off_eliminates):", R[status == "off_eliminates", .N], "\n")
print(R[, .N, by = status][order(-N)])
el = function(tc, pe) is.finite(tc) | (is.finite(pe) & pe < 1e-3)
R[, `:=`(eo = el(t_clear_off, P_end_off), eb = el(t_clear_onB, P_end_onB), es = el(t_clear_onS, P_end_onS))]
cat("Tregs OFF eliminates by the figure's rule (floor or final P_LP < 1e-3):", sum(R$eo), "\n")
cat("final P_LP without Tregs, quartiles:", q(R$P_end_off), " | predicted vs realised persistent agree:", pct(R$P_pred > 1e-3 & !R$eo), "\n")
M0 = R[grepl("^matched", status)]
cat("matched:", nrow(M0), " | |match_err| quantiles 50/90/99/max:", paste(signif(quantile(abs(M0$match_err), c(.5, .9, .99, 1)), 2), collapse = " "), "\n")
cat("absolute damage gap at the match |D_ON* - D_OFF|, quantiles 50/90/99:", paste(signif(quantile(abs(M0$D_win_onS - M0$D_target), c(.5, .9, .99)), 3), collapse = " "), " | D_target median", signif(median(M0$D_target), 3), "\n")
cat("evaluations per unit (n_eval), quartiles:", q(M0$n_eval), " | beta* > 0.0095 (above the sweep):", M0[b_star > 0.0095, .N], " | max beta*", signif(max(M0$b_star), 3), "\n")
cat("jump at beta* (jump_pts > 10):", M0[jump_pts > 10, .N], " | more than one upward crossing:", M0[n_cross_up > 1, .N], "\n")
R[, brake_saves := D_total_onB < D_total_off]
cat("brake alone (ON at beta_r) saves damage over the run:", sum(R$brake_saves), sprintf("(%.1f%%)", 100 * mean(R$brake_saves)), "\n")
print(R[, .N, by = .(brake_saves, status)][order(brake_saves, -N)])
cat("matched_down pairs:", R[status == "matched_down", .N], " all neutral-labelled:", R[status == "matched_down", all(treg_effect == "neutral")],
    " brake costlier in", R[status == "matched_down" & !brake_saves, .N], " | no_match_lo brake costlier in", R[status == "no_match_lo" & !brake_saves, .N], "of", R[status == "no_match_lo", .N], "\n")

## ---- 3. the master filter ----------------------------------------------------------------
sec(3, "the master filter (lib/fig_damage_matched.R)")
MA = fig_dm_matched(verbose = TRUE)                     # one population for every panel
MB = MA                                                 # (the retired ratio cut kept a second one)
lab = U[, .(i, treg_effect)]
MA = merge(MA[, setdiff(names(MA), "treg_effect"), with = FALSE], lab, by = "i")
MB = merge(MB[, setdiff(names(MB), "treg_effect"), with = FALSE], lab, by = "i")
MB = merge(MB, R[, .(i, eo, eb, es)], by = "i")
cat("panel B (all matched pairs of the filter):", nrow(MA), " | beta* < beta_r:", MA[b_star < b_off, .N], " (min ratio", signif(min(MA$b_star / MA$b_off), 4), ")",
    " | ties (ratio < 1.05):", MA[b_star < 1.05 * b_off, .N], " | ratio >= 1.05:", nrow(MB), "\n")
cat("panel B per level (median beta*, IQR, median ratio):\n")
print(MA[, .(n = .N, med = signif(median(b_star), 3), q1 = signif(quantile(b_star, .25), 3), q3 = signif(quantile(b_star, .75), 3), ratio = round(median(b_star / b_off), 2),
             ratio_q90 = round(quantile(b_star / b_off, .9), 2)), by = .(b_off = signif(b_off, 3))][order(b_off)])
cat("ratio over all pairs: median", round(median(MA$b_star / MA$b_off), 2), " IQR", q(MA$b_star / MA$b_off, c(.25, .75)), " 90th", round(quantile(MA$b_star / MA$b_off, .9), 2), " max", round(max(MA$b_star / MA$b_off), 1), "\n")
cat("by Treg effect at own beta_r (pairs, median ratio, share of ties):\n")
print(MA[, .(pairs = .N, ratio_med = round(median(b_star / b_off), 2), ties = round(100 * mean(b_star < 1.05 * b_off), 1)), by = treg_effect][order(-pairs)])

## ---- 4. the first peak, the window, the trajectories -----------------------------------
sec(4, "the first peak (pipeline/42_first_peak.R), the matching window, the trajectory subset")
P1 = as.data.table(readRDS("./data/fig06_first_peak.rds"))
MA = merge(MA, P1, by = "i")
cat("first local maximum of P_LP, OFF(beta_r): median time", signif(median(MA$t_peak1_off), 3), " quartiles of the value:", q(MA$peak1_off), "\n")
cat("series with no local maximum: OFF", MA[is.na(peak1_off), .N], " ON*", MA[is.na(peak1_onS), .N],
    " | first maximum within 1% of the final load: OFF", MA[no_transient_off == TRUE, .N], " ON*", MA[no_transient_onS == TRUE, .N], "\n")
cat("first maximum lower with Tregs at beta*:", pct(MA$peak1_onS < MA$peak1_off),
    " median reduction", round(100 * median(1 - MA$peak1_onS / MA$peak1_off), 1), "%  max", round(100 * max(1 - MA$peak1_onS / MA$peak1_off), 1), "%\n")
cat("the retired measure (maximum over the run) lower with Tregs:", pct(MA$gmax_onS < MA$gmax_off),
    " | pairs where it is higher:", MA[gmax_onS > gmax_off, .N], " of which the ON* maximum is its attractor:", MA[gmax_onS > gmax_off & no_transient_onS == TRUE, .N], "\n")
cat("window t_conv quartiles:", q(MA$win), " range", signif(range(MA$win), 4), " | T_conv of OFF(beta_r) median", signif(median(MA$tconv_off), 3),
    " | ON(beta*) sets the window in", pct(MA$tconv_onS > MA$tconv_off), "\n")
X = readRDS("./data/fig04b_v3_trajectories.rds"); TR = X$traj; S3 = X$units
## the seed-1 draw of 1,000 pairs is taken from the current filter, so every stored pair is
## one the figure draws; the subsetting below is kept as a guard
n_stored = nrow(S3); S3 = S3[i %in% MA$i]; TR = TR[i %in% S3$i]
cat("trajectory units stored:", n_stored, " drawn by B and C (in the filter):", nrow(S3),
    " | runs complete:", TR[, .N, by = .(i, run)][, all(N == 2001)], "\n")
TP = TR[run == "OFF" & time > 0, { k = which(diff(sign(diff(PLP))) < 0) + 1L
        .(tp = if (length(k)) time[k[1]] else time[which.max(PLP)]) }, by = i]
cat("t_peak recomputed on the 2000-point grid, quartiles:", q(TP$tp), " range", signif(range(TP$tp), 3), "\n")
sup = TR[time > 0, { o = .SD[run == "OFF"]; a = .SD[run == "ON_same"]
        k = which(diff(sign(diff(o$PLP))) < 0) + 1L; tp = if (length(k)) o$time[k[1]] else o$time[which.max(o$PLP)]
        d = o$M1 - a$M1
        .(frac = approx(o$time, d, xout = tp, rule = 2)$y / max(d)) }, by = i]
cat("share of the eventual M1 suppression reached at the first peak: median", round(100 * median(sup$frac), 2),
    "%  90th", round(100 * quantile(sup$frac, .9), 1), "%  above 10% in", sup[frac > 0.1, .N], "of", nrow(sup), "pairs\n")

## ---- 5. the panels' numbers --------------------------------------------------------------
sec(5, "panel C (ratio) and the first-peak claim")
MA[, `:=`(peak_red = 100 * (1 - peak1_onS / peak1_off), win_red = 100 * (1 - P_win_onS / P_win_off))]
cat("pairs:", nrow(MA), " | first maximum lower at beta*:", pct(MA$peak_red > 0), " median", round(median(MA$peak_red), 1), "%\n")
cat("window integral lower:", pct(MA$win_red > 0), " median", round(median(MA$win_red), 1), "%  max", round(max(MA$win_red), 1), "%",
    " | whole-run integral lower:", pct(MA$auc_P_onS < MA$auc_P_off), "\n")
G = exp(seq(log(0.1), log(100), length.out = 400))
rat = merge(TR, TP, by = "i")[time > 0, { o = .SD[run == "OFF"]; a = .SD[run == "ON_same"]; s = .SD[run == "ON_star"]
  x = log(o$time / tp[1]); f = function(y) approx(x, log(pmax(y, 1e-12) / pmax(o$PLP, 1e-12)), xout = log(G), rule = 2)$y
  .(tau = G, same = f(a$PLP), star = f(s$PLP)) }, by = i]
CC = rat[, .(same = exp(median(same)), star = exp(median(star))), by = tau]
for (tt in c(0.5, 1, 2, 5, 10, 50, 100)) { k = which.min(abs(CC$tau - tt)); cat(sprintf("  tau %5.1f: same %.3f  star %.3f\n", tt, CC$same[k], CC$star[k])) }
cat("minimum of the median star ratio:", round(min(CC$star), 3), "at tau", round(CC$tau[which.min(CC$star)], 2), "\n")

## ---- 6. fates and the Sankey -------------------------------------------------------------
sec(6, "fates and routes (Extended Data Fig. 5)")
MA = merge(MA, R[, .(i, eo, eb, es)], by = "i")
MA[, fate := fifelse(eo & es, "both_elim", fifelse(!eo & !es, "both_persist", fifelse(es, "only_on", "only_off")))]
print(MA[, .N, by = fate])
P = MA[fate == "both_persist"]
P[, `:=`(pk = peak1_onS < peak1_off, it = P_win_onS < P_win_off, ss = P_end_onS < P_end_off)]
print(P[, .N, by = .(first_peak_less = pk, int_less = it, ss_less = ss)][order(-N)])
cat("persistent pairs: steady state higher with Tregs in", pct(!P$ss), " median difference OFF - ON", round(median(P$P_end_off - P$P_end_onS), 2),
    " | residual load OFF median", round(median(P$P_end_off), 2), "\n")
cat("effect sizes, all pairs: first peak OFF-ON median", round(median(MA$peak1_off - MA$peak1_onS), 2), " range", round(range(MA$peak1_off - MA$peak1_onS), 2),
    " | integral % median", round(median(MA$win_red), 1), " range", round(range(MA$win_red)),
    " | steady state OFF-ON median", round(median(MA$P_end_off - MA$P_end_onS), 2), " range", round(range(MA$P_end_off - MA$P_end_onS), 1), "\n")
cat("steady state higher with Tregs, all pairs:", pct(MA$P_end_onS > MA$P_end_off), " above one unit in", pct(MA$P_end_onS - MA$P_end_off > 1),
    " | epithelium higher with Tregs:", pct(MA$Eh_end_onS > MA$Eh_end_off), " (median", round(median(MA$Eh_end_onS), 1), "vs", round(median(MA$Eh_end_off), 1), ")\n")
cat("worst integrated excess: ON*/OFF =", round(max(MA$P_win_onS / MA$P_win_off), 1), "\n")

## ---- 7. the eliminations and the fourth arm ------------------------------------------
sec(7, "eliminations and the fourth arm (pipeline/48_gated_fourth_arm.R)")
E = MA[fate == "only_on"]
cat("only Tregs-ON at beta* eliminates:", nrow(E), " in", uniqueN(E[, .(tau_id, Kg, r_L_P)]), "cells | the brake alone eliminates in", MA[eb == TRUE, .N], "\n")
cat("beta*/beta_r of the eliminations, quartiles:", q(E$b_star / E$b_off), " | final P_LP without Tregs", q(E$P_end_off), " | t_clear", q(E$t_clear_onS), "\n")
cat("final E_h, ON* - OFF:", q(E$Eh_end_onS - E$Eh_end_off), "\n")
F4 = as.data.table(readRDS("./data/aux_gated_fourth_arm.rds"))[i %in% E$i]
cat("fourth arm available for", nrow(F4), "of", nrow(E), " | OFF at beta* eliminates too:", sum(F4$offstar_elim), " | cost D(OFF,beta*)/D(OFF,beta_r) quartiles", q(F4$cost_noTreg),
    " | > 1.2:", sum(F4$cost_noTreg > 1.2), " > 1.5:", sum(F4$cost_noTreg > 1.5), " | E_h of OFF at beta* minus OFF at beta_r, median", round(median(F4$Eh_end_OFF_star - F4$Eh_end_OFF), 1), "\n")
print(F4[, .(n = .N, cost_med = round(median(cost_noTreg), 2)), by = treg_effect])

## ---- 8. the demonstration unit ---------------------------------------------------------
sec(8, "the demonstration unit (figure_scripts/Fig_4.R D; lib/fig_demo_case.R)")
DEMO = as.integer(Sys.getenv("UNIT", "4108"))
u = MA[i == DEMO]; f = as.data.table(readRDS("./data/aux_gated_fourth_arm.rds"))[i == DEMO]
cat(sprintf("unit %d: tau %d, K_g %g, r_P %.2f, E_h(0) %g, beta_r %.4g, beta* %.4g (%.2fx), t_peak %.3g, t_conv %s, Treg effect at own beta_r: %s\n",
            u$i, u$tau_id, u$Kg, u$r_L_P, u$Eh0, u$b_off, u$b_star, u$b_star / u$b_off, u$t_peak1_off, format(round(u$win), big.mark = ","), u$treg_effect))
cat(sprintf("stored: P_end OFF %.1f, ON(beta_r) %.1f, ON(beta*) %.3g; t_clear ON(beta*) %.0f; E_h end OFF %.1f ON(beta_r) %.1f ON(beta*) %.1f\n",
            u$P_end_off, u$P_end_onB, u$P_end_onS, u$t_clear_onS, u$Eh_end_off, u$Eh_end_onB, u$Eh_end_onS))
cat(sprintf("fourth arm (400-pt grid): OFF(beta*) eliminates %s at t %.0f, E_h end %.1f, M1 end OFF %.1f ON* %.1f OFF* %.1f, cost %.2f\n",
            f$offstar_elim, f$t_clear_OFF_star, f$Eh_end_OFF_star, f$M1_end_OFF, f$M1_end_ON_star, f$M1_end_OFF_star, f$cost_noTreg))
