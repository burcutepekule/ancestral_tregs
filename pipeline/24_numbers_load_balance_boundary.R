## ============================================================================
## pipeline/24_numbers_load_balance_boundary.R -- every number quoted in the methods for
## Fig. 3C, Extended Data Figs. 2C and 3, and Supplementary Fig. 1 (docs/figure_methods.pdf).
## Reads only ./data/; integrates nothing. Prints each block and writes them long
## (block, group, item, value) to
##   ./source_data/load_balance_numbers.csv   (Fig. 3C, Extended Data Figs. 2C, 3: gate, load-balance table)
##   ./source_data/boundary_numbers.csv       (Supplementary Fig. 1: pool, ledger, statistics)
## in : rescue_<br>.rds, skeleton_<br>.rds, selected_big.rds, load_balance_<br>.rds,
##      corridor_tzero_sterile_round1.rds, corridor_tzero_pathogenic_help.rds,
##      corridor_q_sterile_round1.rds, corridor_q_pathogenic_help.rds,
##      corridor_lo_off_pathogenic_help_disagree.rds, tau grids
## RUN: Rscript pipeline/24_numbers_load_balance_boundary.R
## ============================================================================
suppressMessages(library(data.table))
SRC = Sys.getenv("SRC", "."); setwd(SRC)
CK = c("tau_id","Kg","beta_ros","r_L_P"); BR = c("sterile", "pathogenic_help")
STEP = 10 / 2^ceiling(log2(10 / 0.1))       # rescue-gap lattice step, 0.078125
OUT2 = list(); OUT3 = list()
put2 = function(block, group, item, value) OUT2[[length(OUT2) + 1]] <<- data.table(block, group, item, value = as.character(value))
put3 = function(block, group, item, value) OUT3[[length(OUT3) + 1]] <<- data.table(block, group, item, value = as.character(value))
SB = as.data.table(readRDS("./data/selected_big.rds"))
pct = function(x) sprintf("%.1f", 100 * mean(x))

## ---- 1. load-balance gate, step by step ----------------------------------
for (br in BR) {
  gmin = if (br == "sterile") 5 else 2
  R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", br)))
  S = as.data.table(readRDS(sprintf("./data/skeleton_%s.rds", br)))
  put2("gate", br, "rescue rows", nrow(R))
  for (m in c("bistable", "mono_B")) put2("gate", br, paste("rescue rows", m), R[rescue_mode == m, .N])
  R = merge(R, unique(S[, .(tau_id, Kg, beta_ros, r_L_P, batch, row_ind)]), by = CK)
  put2("gate", br, "joined to skeleton", nrow(R))
  Bg = SB[, .(dev_max_on = max(dev_max_on, na.rm = TRUE), dev_max_off = max(dev_max_off, na.rm = TRUE)), by = CK]
  R = merge(R, Bg, by = CK, all.x = TRUE)
  R[, dev_ok := is.finite(dev_max_on) & dev_max_on <= 1e-2 & is.finite(dev_max_off) & dev_max_off <= 1e-2]
  for (m in c("bistable", "mono_B")) {
    X = R[rescue_mode == m]; put2("gate", br, paste(m, "start"), nrow(X))
    X = X[dev_ok == TRUE]; put2("gate", br, paste(m, "after worst-probe deviation <= 1e-2 on both rays"), nrow(X))
    if (m == "bistable") {
      X = X[dir_off == "increasing" & dir_on == "increasing"]; put2("gate", br, "bistable after both directions increasing", nrow(X))
      X = X[gap_censored == FALSE]; put2("gate", br, "bistable after gap not censored", nrow(X))
      X = X[rescue_gap >= gmin]; put2("gate", br, sprintf("bistable after rescue_gap >= %g", gmin), nrow(X))
    } else {
      X = X[dir_off == "increasing"]; put2("gate", br, "mono_B after OFF direction increasing", nrow(X))
      X = X[is.finite(Eh0_off)]; put2("gate", br, "mono_B after finite Eh0_off", nrow(X))
    }
    put2("gate", br, paste(m, "in the Treg-help class"), X[treg_effect == "help", .N])
    put2("gate", br, paste(m, "drawn"), min(500L, nrow(X)))
  }
}

## ---- 2. the load-balance table (full 0.5 grid summaries) --------------------
TAB = list()
for (br in BR) {
  ## one effect class in both rows of the load-balance figures (Fig. 3C, Extended Data Figs. 2C and 3): the pathogenic
  ## branch is Treg-help by construction, the sterile gate keeps help and neutral
  ## alike, so the summaries below are restricted to help. lib/fig_load_balance.R
  ## filters the same way and the two must stay in step.
  D0 = as.data.table(readRDS(sprintf("./data/load_balance_%s.rds", br)))
  D0[, tau_id := batch * 10000L + row_ind]
  kk = function(X) paste(X$tau_id, X$Kg, X$beta_ros, X$r_L_P)
  HELPK = kk(unique(SB[treg_effect == "help", c(CK), with = FALSE]))
  D = D0[kk(D0) %in% HELPK]
  for (cl in c("bistable", "mono_B"))
    put2("gate", br, paste(cl, "integrated / summarised after the Treg-help filter"),
         sprintf("%d / %d", uniqueN(D0[cls == cl]$cid), uniqueN(D[cls == cl]$cid)))
  D[, `:=`(L0K = L0 / Kg, LminK = L_min / Kg, LendK = L_end / Kg)]
  W = dcast(D, cid + cls ~ cond, value.var = c("L0K","LminK","LendK","frac_removed","dL_turns_pos_in_lowg","ends_lowg","resolves","steady_end",
                                                  "M1_at_off_min","L_at_off_min","M1_max_lowg","L_ever_le_Kg","dL_min","g_end_over_a0",
                                                  "PLP_end_over_K","PLP_max_over_K","t_end","t_Lmin","t95","Kg"))
  setnames(W, function(x) sub("_Tregs OFF$", "_off", sub("_Tregs ON$", "_on", x)))
  for (cl in c("bistable", "mono_B")) {
    w = W[cls == cl]; g = paste(br, ifelse(cl == "mono_B", "ON_mono", "bistable"))
    T = list(
      n = nrow(w),
      reached_TMAX_both = sum(w$t_end_off >= 12000 & w$t_end_on >= 12000),
      starts_L_gt_Kg = sum(w$L0K_off > 1 & w$L0K_on > 1),
      L0_over_Kg_median = sprintf("%.2f", median(w$L0K_off)),
      L0_over_Kg_range = sprintf("%.2f-%.2f", min(w$L0K_off), max(w$L0K_off)),
      OFF_min_L_above_Kg = sum(w$LminK_off > 1),
      OFF_min_L_over_Kg_median = sprintf("%.2f", median(w$LminK_off)),
      OFF_median_load_removed_before_min_pct = sprintf("%.1f", 100 * median(w$frac_removed_off)),
      OFF_time_of_min_median = sprintf("%.1f", median(w$t_Lmin_off)),
      OFF_kill_gt_clear_while_L_gt_Kg = sum(w$dL_turns_pos_in_lowg_off),
      OFF_ends_L_gt_Kg = sum(w$ends_lowg_off),
      OFF_end_L_over_Kg_median = sprintf("%.2f", median(w$LendK_off)),
      OFF_end_L_over_Kg_IQR = sprintf("%.2f-%.2f", quantile(w$LendK_off, .25), quantile(w$LendK_off, .75)),
      OFF_at_rest = sum(w$steady_end_off),
      OFF_resolves_g095_ever = sum(w$resolves_off),
      ON_L_at_OFF_min_over_Kg_median = sprintf("%.2f", median(w$L_at_off_min_on / w$Kg_on)),
      OFF_L_at_OFF_min_over_Kg_median = sprintf("%.2f", median(w$L_at_off_min_off / w$Kg_off)),
      ON_M1_at_OFF_min_median = sprintf("%.1f", median(w$M1_at_off_min_on)),
      OFF_M1_at_OFF_min_median = sprintf("%.1f", median(w$M1_at_off_min_off)),
      ON_lower_M1_at_OFF_min = sum(w$M1_at_off_min_on < w$M1_at_off_min_off),
      ON_lower_peak_M1_while_L_gt_Kg = sum(w$M1_max_lowg_on < w$M1_max_lowg_off),
      peak_M1_while_L_gt_Kg_median_OFF = sprintf("%.1f", median(w$M1_max_lowg_off)),
      peak_M1_while_L_gt_Kg_median_ON = sprintf("%.1f", median(w$M1_max_lowg_on)),
      ON_kill_gt_clear_while_L_gt_Kg = sum(w$dL_turns_pos_in_lowg_on),
      ON_L_reaches_Kg = sum(w$L_ever_le_Kg_on),
      ON_deeper_min_dLdt = sum(w$dL_min_on < w$dL_min_off),
      ratio_min_dLdt_ON_over_OFF_median = sprintf("%.2f", median(w$dL_min_on / w$dL_min_off)),
      ON_resolves_g095_ever = sum(w$resolves_on),
      ON_ends_with_g095 = sum(w$g_end_over_a0_on >= 0.95),
      ON_ends_L_lt_Kg = sum(!w$ends_lowg_on),
      ON_end_L_over_Kg_median = sprintf("%.3f", median(w$LendK_on)),
      ON_at_rest = sum(w$steady_end_on),
      ON_t95_median = sprintf("%.1f", median(w$t95_on, na.rm = TRUE)))
    if (br != "sterile") T = c(T, list(
      OFF_PLP_end_lt_1e4_pct = pct(w$PLP_end_over_K_off < 1e-4), OFF_PLP_end_lt_1e2_pct = pct(w$PLP_end_over_K_off < 1e-2),
      ON_PLP_end_gt_1e2_pct = pct(w$PLP_end_over_K_on > 1e-2), ON_PLP_end_gt_1e1_pct = pct(w$PLP_end_over_K_on > 0.1),
      OFF_PLP_end_median = signif(median(w$PLP_end_over_K_off), 3), ON_PLP_end_median = signif(median(w$PLP_end_over_K_on), 3),
      OFF_PLP_peak_median = sprintf("%.3f", median(w$PLP_max_over_K_off)), ON_PLP_peak_median = sprintf("%.3f", median(w$PLP_max_over_K_on))))
    for (k in names(T)) put2("table", g, k, T[[k]])
    TAB[[g]] = T
  }
  ## the exceptions behind open item 13, listed
  EX = D[(cond == "Tregs ON" & !(resolves & !ends_lowg & g_end_over_a0 >= 0.95)) | (cond == "Tregs OFF" & (resolves | !ends_lowg))]
  if (nrow(EX)) { cat(sprintf("\n[numbers] %s: trajectories that are not 'OFF chronic / ON resolved at rest with g >= 0.95 a0'\n", br))
    print(EX[, .(cid, cls, cond, L0K = round(L0K, 2), LminK = round(LminK, 3), LendK = round(LendK, 3), g_end = round(g_end_over_a0, 3), resolves, ends_lowg, steady_end, t95)])
    for (j in seq_len(nrow(EX))) put2("exceptions", br, sprintf("cid %d %s %s", EX$cid[j], EX$cls[j], EX$cond[j]),
                                     sprintf("L_end/Kg %.3f; min L/Kg %.3f; g_end/a0 %.3f; resolves %s; at rest %s", EX$LendK[j], EX$LminK[j], EX$g_end_over_a0[j], EX$resolves[j], EX$steady_end[j])) }
  ## the drawn subset: points with L/Kg <= 1e-3 are not drawn
  P = D[, .(below = sum(pop[[1]]$L_over_Kg <= 1e-3)), by = .(cid, cls, cond)]
  put2("drawing", br, "trajectories with stored points at L/Kg <= 1e-3 (not drawn)", P[below > 0, .N])
}
TT = rbindlist(lapply(names(TAB), function(g) data.table(group = g, item = names(TAB[[g]]), value = unlist(lapply(TAB[[g]], as.character)))))
cat("\n[numbers] load-balance table\n"); print(dcast(TT, item ~ group, value.var = "value"), nrows = 200)

## ---- 3. ON_mono off the line (open item 14) and the rescue-gap floor (item 15) ----
for (br in BR) {
  R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", br)))
  M = merge(R[rescue_mode == "mono_B", ..CK], SB[, c(CK, "n_field_stable_on"), with = FALSE], by = CK)
  put2("on_mono_field", br, "ON_mono with >= 2 stable Tregs-ON fixed points found (n_field_stable_on >= 2)", sprintf("%d of %d", M[n_field_stable_on >= 2, .N], nrow(M)))
  G = R[rescue_mode == "bistable" & dir_off == "increasing" & dir_on %in% "increasing" & !gap_censored]
  put2("gap_floor", br, "bistable with measured gap, both directions increasing", nrow(G))
  put2("gap_floor", br, "gap exactly 0 (pct)", pct(G$rescue_gap == 0))
  put2("gap_floor", br, "gap < 0 (n)", G[rescue_gap < 0, .N])
  put2("gap_floor", br, sprintf("gap <= one lattice step %g (pct)", STEP), pct(G$rescue_gap <= STEP * 1.01))
  put2("gap_floor", br, "gap max", sprintf("%.3f", max(G$rescue_gap)))
  put2("gap_floor", br, "gap lattice check: max distance of gap/step from an integer", signif(max(abs(G$rescue_gap / STEP - round(G$rescue_gap / STEP))), 3))
}

## ---- 4. Supplementary Fig. 1: corridor pool, round 1, notes ------------------------------
for (br in BR) {
  tag = if (br == "sterile") "_round1" else ""
  R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", br))); U = R[rescue_mode %in% c("bistable", "mono_B")]
  T0 = as.data.table(readRDS(sprintf("./data/corridor_tzero_%s%s.rds", br, tag)))
  K = unique(T0[, ..CK]); K = merge(K, U[, c(CK, "rescue_mode", "rescue_gap", "gap_censored"), with = FALSE], by = CK)
  put3("pool", br, "pool (bistable + mono_B)", nrow(U)); put3("pool", br, "pool bistable", U[rescue_mode == "bistable", .N]); put3("pool", br, "pool mono_B", U[rescue_mode == "mono_B", .N])
  put3("pool", br, "measured combinations", sprintf("%d (%.1f%%)", nrow(K), 100 * nrow(K) / nrow(U)))
  put3("pool", br, "measured bistable", K[rescue_mode == "bistable", .N]); put3("pool", br, "measured mono_B", K[rescue_mode == "mono_B", .N])
  ## is the measured subsample representative of the pool? (it is a prefix of a seeded permutation)
  put3("pool", br, "mono_B share pool / measured (pct)", sprintf("%.2f / %.2f", 100 * mean(U$rescue_mode == "mono_B"), 100 * mean(K$rescue_mode == "mono_B")))
  gp = U[rescue_mode == "bistable" & !gap_censored & is.finite(rescue_gap), rescue_gap]; gk = K[rescue_mode == "bistable" & !gap_censored & is.finite(rescue_gap), rescue_gap]
  put3("pool", br, "bistable rescue gap: pct > one step, pool / measured", sprintf("%.2f / %.2f", 100 * mean(gp > STEP * 1.01), 100 * mean(gk > STEP * 1.01)))
  put3("pool", br, "bistable rescue gap: KS test measured vs pool, D and p", { k = suppressWarnings(ks.test(gk, gp)); sprintf("D = %.4f, p = %.3g", k$statistic, k$p.value) })
  nt = T0[, .(note = note[1]), by = CK][, .N, by = note]
  for (j in seq_len(nrow(nt))) put3("notes", br, paste("per-combination outcome:", nt$note[j]), nt$N[j])
  OK = T0[note == "ok"]
  put3("base", br, "base_frac distribution", paste(OK[, .N, by = base_frac][order(base_frac), sprintf("%g: %d rows", base_frac, N)], collapse = "; "))
  put3("base", br, "combinations with base rung 0.02", sprintf("%d of %d", uniqueN(OK[base_frac == 0.02, ..CK]), uniqueN(OK[, ..CK])))
}

## ---- 5. Supplementary Fig. 1B/C ledger and statistics ------------------------------------
source("./lib/fig_boundary_gain.R")
DIS = as.data.table(readRDS("./data/corridor_lo_off_pathogenic_help_disagree.rds"))
for (br in BR) {
  tag = if (br == "sterile") "_round1" else ""
  gp = gain_panel(br, tag)
  Q = as.data.table(readRDS(sprintf("./data/corridor_q_%s%s.rds", br, tag)))
  put3("ledger", br, "bistable combinations with note ok", uniqueN(Q[rescue_mode == "bistable", ..CK]))
  CNT = gp$counts[order(match(coord, c("Ea","En","M1","R","D")))]
  for (j in seq_len(nrow(CNT))) put3("ledger", br, paste("coord", CNT$coord[j]),
       sprintf("rows %d; both walls (measured) %d; G>0 %d; gap finite, uncensored, >0 %d; at floor %d; plotted %d",
               CNT$rows[j], CNT$measured[j], CNT$gain_pos[j], CNT$gap_ok[j], CNT$at_floor[j], CNT$plotted[j]))
  QK = Q[coord %in% c("Ea","En","M1","R","D") & rescue_mode == "bistable", .N, by = .(coord, q_kind)]
  for (j in seq_len(nrow(QK))) put3("q_kind", br, paste(QK$coord[j], QK$q_kind[j]), QK$N[j])
  QM = Q[coord %in% c("Ea","En","M1","R","D") & rescue_mode == "mono_B", .(n = .N, on_absent = sum(q_kind == "on_absent"), measured = sum(q_kind == "measured")), by = coord]
  for (j in seq_len(nrow(QM))) put3("mono_B", br, paste(QM$coord[j], "rows / no ON wall in range / both walls"), sprintf("%d / %d (%.1f%%) / %d", QM$n[j], QM$on_absent[j], 100 * QM$on_absent[j] / QM$n[j], QM$measured[j]))
  GF = gp$rows
  ST = GF[, {
    ct = suppressWarnings(cor.test(gain, rescue_gap, method = "spearman", exact = FALSE))
    f = lm(log10(rescue_gap) ~ log10(gain))
    .(n = .N, rho = unname(ct$estimate), p = ct$p.value, slope = unname(coef(f)[2]), intercept = unname(coef(f)[1]),
      G_med = median(gain), G_q25 = quantile(gain, .25), G_q75 = quantile(gain, .75),
      gap_med = median(rescue_gap)) }, by = coord][order(-rho)]
  for (j in seq_len(nrow(ST))) put3("stats", br, ST$coord[j],
       sprintf("n %d; rho %.3f; p %.3g; OLS log10(gap) = %.3f + %.3f log10(G); G median %.1f%% (IQR %.1f-%.1f)",
               ST$n[j], ST$rho[j], ST$p[j], ST$intercept[j], ST$slope[j], ST$G_med[j], ST$G_q25[j], ST$G_q75[j]))
  put3("stats", br, "distinct combinations plotted (any coordinate)", uniqueN(GF[, ..CK]))
  put3("stats", br, "rows pooled over coordinates / drawn", sprintf("%d / %d", nrow(GF), min(9000L, nrow(GF))))
  cat(sprintf("\n[numbers] Supplementary Fig. 1 %s\n", br)); print(ST)
  if (br == "pathogenic_help") {
    ## open item 7: the combinations whose PLP span disagrees with the sweep
    GFd = GF[!DIS, on = CK]
    put3("disagree", br, "combinations in the disagreement file", nrow(DIS))
    put3("disagree", br, "plotted rows from those combinations", nrow(GF) - nrow(GFd))
    S2 = GFd[, .(n = .N, rho = suppressWarnings(cor(gain, rescue_gap, method = "spearman"))), by = coord][order(-rho)]
    for (j in seq_len(nrow(S2))) put3("disagree", br, paste("rho without them", S2$coord[j]), sprintf("n %d; rho %.3f", S2$n[j], S2$rho[j]))
    LO = as.data.table(readRDS("./data/corridor_lo_off_pathogenic_help.rds"))
    T0 = as.data.table(readRDS("./data/corridor_tzero_pathogenic_help.rds"))[note == "ok"]
    CH = merge(unique(T0[, c(CK, "coord", "span"), with = FALSE]), LO, by = c(CK, "coord")); CH[, rel := abs(abs(hi_off - lo_off) - span) / pmax(abs(span), 1e-12)]
    put3("disagree", br, "coordinates on which the span differs by > 1%", paste(CH[rel > 0.01, .N, by = coord][, sprintf("%s %d", coord, N)], collapse = "; "))
    put3("disagree", br, "PLP: median |span difference| in those rows", signif(CH[rel > 0.01, median(abs(abs(hi_off - lo_off) - span))], 3))
  }
}
dir.create("./source_data", showWarnings = FALSE); fwrite(rbindlist(OUT2), "./source_data/load_balance_numbers.csv"); fwrite(rbindlist(OUT3), "./source_data/boundary_numbers.csv")
cat("\n[numbers] wrote ./source_data/load_balance_numbers.csv and ./source_data/boundary_numbers.csv\n")
