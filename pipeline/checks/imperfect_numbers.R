## ============================================================================
## pipeline/checks/imperfect_numbers.R -- every number quoted in
## docs/figure_methods.pdf (Fig. 5, Extended Data Figs. 6 and 7),
## recomputed from ./data. No integration.
## in : ./data/rescue_<branch>.rds, ./data/skeleton_<branch>.rds,
##      ./data/selected_big.rds, ./data/fig05_units*.{rds,csv},
##      ./data/fig05_raw_all.rds, ./data/fig05_cells.rds,
##      ./data/fig05_switch_prob_cells.rds, ./data/fig05_switch_prob_units.rds
## RUN: Rscript pipeline/checks/imperfect_numbers.R > imperfect_numbers.txt
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages(library(data.table))
options(width = 200)
sec = function(s) cat("\n==================== ", s, " ====================\n", sep = "")
CK = c("tau_id", "Kg", "beta_ros", "r_L_P")

## ---- 1. eligibility, step by step -----------------------------------------
sec("1. eligibility")
BG = as.data.table(readRDS("./data/selected_big.rds"))[
  , .(dev_max_on = max(dev_max_on, na.rm = TRUE), dev_max_off = max(dev_max_off, na.rm = TRUE)), by = CK]
for (br in c("sterile", "pathogenic_help")) {
  GAP_MIN = if (br == "sterile") 5 else 2
  R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", br)))
  cat(sprintf("\n%s: rescue table %d rows (bistable %d, mono_B %d)\n", br, nrow(R),
              R[rescue_mode == "bistable", .N], R[rescue_mode == "mono_B", .N]))
  S = as.data.table(readRDS(sprintf("./data/skeleton_%s.rds", br)))
  R = merge(R, unique(S[, .(tau_id, Kg, beta_ros, r_L_P, batch, row_ind)]), by = CK)
  cat(sprintf("  after join to skeleton: %d\n", nrow(R)))
  R = merge(R, BG, by = CK, all.x = TRUE)
  R = R[is.finite(dev_max_on) & dev_max_on <= 1e-2 & is.finite(dev_max_off) & dev_max_off <= 1e-2]
  cat(sprintf("  deviation <= 1e-2 in both arms: %d (bistable %d, mono_B %d)\n", nrow(R),
              R[rescue_mode == "bistable", .N], R[rescue_mode == "mono_B", .N]))
  b = R[rescue_mode == "bistable"]
  cat(sprintf("  bistable: gap >= %g %d | + both increasing %d | + not censored %d\n", GAP_MIN,
              b[rescue_gap >= GAP_MIN, .N],
              b[rescue_gap >= GAP_MIN & dir_off == "increasing" & dir_on == "increasing", .N],
              b[rescue_gap >= GAP_MIN & dir_off == "increasing" & dir_on == "increasing" & !gap_censored, .N]))
  m = R[rescue_mode == "mono_B"]
  cat(sprintf("  mono_B: dir_off increasing %d | + Eh0_off finite %d\n",
              m[dir_off == "increasing", .N], m[dir_off == "increasing" & is.finite(Eh0_off), .N]))
  E = rbind(b[rescue_gap >= GAP_MIN & dir_off == "increasing" & dir_on == "increasing" & !gap_censored],
            m[dir_off == "increasing" & is.finite(Eh0_off)], fill = TRUE)
  E[, stratum := paste(Kg, beta_ros, r_L_P)]
  print(E[, .(rows = .N, strata = uniqueN(stratum)), by = rescue_mode])
  cat(sprintf("  strata occupied (modes pooled): %d; rows per stratum min %d median %g max %d\n",
              uniqueN(E$stratum), min(E[, .N, by = stratum]$N), median(E[, .N, by = stratum]$N),
              max(E[, .N, by = stratum]$N)))
  cat(sprintf("  levels: Kg {%s} beta_ros {%s} r_L_P {%s}\n", paste(sort(unique(E$Kg)), collapse = ","),
              paste(sort(unique(E$beta_ros)), collapse = ","), paste(sort(unique(E$r_L_P)), collapse = ",")))
}

## ---- 2. the draw -----------------------------------------------------------
sec("2. drawn units")
U = fread("./data/fig05_units.csv")
print(U[, .(units = .N, strata = uniqueN(stratum), Eh0_q25 = quantile(Eh0, .25), Eh0_med = median(Eh0),
            Eh0_q75 = quantile(Eh0, .75), Eh0_min = min(Eh0), Eh0_max = max(Eh0)), by = .(branch, rescue_mode)][order(branch, rescue_mode)])
print(U[, .(units = .N, strata = uniqueN(stratum), waves = paste(sort(unique(wave)), collapse = ",")), by = branch])
cat("units per stratum, by branch:\n"); print(U[, .N, by = .(branch, stratum)][, .(strata = .N), by = .(branch, N)][order(branch, N)])
cat("bistable units per stratum among strata holding any bistable eligible row (sterile / help):\n")
print(U[, .(bi = sum(rescue_mode == "bistable")), by = .(branch, stratum)][, .N, by = .(branch, bi)][order(branch, bi)])
cat(sprintf("Eh0 inside rescued interval (bistable: Eh0_on < Eh0 < Eh0_off): %d of %d\n",
            U[rescue_mode == "bistable" & Eh0 > Eh0_on & Eh0 < Eh0_off, .N], U[rescue_mode == "bistable", .N]))
cat(sprintf("mono_B Eh0 = Eh0_off / 2 exactly: %d of %d\n",
            U[rescue_mode == "mono_B" & abs(Eh0 - Eh0_off / 2) < 1e-12, .N], U[rescue_mode == "mono_B", .N]))
cat("rescue_gap of drawn bistable units:\n"); print(U[rescue_mode == "bistable", as.list(quantile(rescue_gap, c(0, .25, .5, .75, 1))), by = branch])
cat("units by Kg x branch:\n"); print(dcast(U[, .N, by = .(branch, Kg)], Kg ~ branch, value.var = "N"))
cat("units by r_L_P (help):\n"); print(U[branch == "pathogenic_help", .N, by = r_L_P][order(r_L_P)])

## ---- 3. the sweep: what came back -----------------------------------------
sec("3. sweep coverage")
B = as.data.table(readRDS("./data/fig05_raw_all.rds"))
B = merge(B, U[, .(unit, branch, wave)], by = "unit")
cat(sprintf("rows %d | treated %d | untreated %d | units %d\n", nrow(B), B[theta_treg == 1, .N],
            B[theta_treg == 0, .N], uniqueN(B$unit)))
cat(sprintf("solver_ok FALSE rows present: %d\n", B[solver_ok == FALSE, .N]))
IT = B[theta_treg == 1, .(n = .N), by = .(unit, rep, branch)]
cat(sprintf("(unit, rep) items present %d of %d; runs per present item: %s\n", nrow(IT), nrow(U) * 10L,
            paste(names(table(IT$n)), table(IT$n), sep = ":", collapse = " ")))
ALL = CJ(unit = U$unit, rep = 1:10); ALL = merge(ALL, U[, .(unit, branch, wave)], by = "unit")
MISS = ALL[!IT, on = .(unit, rep)]
cat(sprintf("items absent: %d (%d runs); by branch:\n", nrow(MISS), nrow(MISS) * 121L)); print(MISS[, .N, by = .(branch, wave)])
cat(sprintf("units missing >= 1 item: %d; max items missing in one unit: %d\n",
            uniqueN(MISS$unit), max(MISS[, .N, by = unit]$N)))
OFFU = B[theta_treg == 0, unique(unit)]
cat(sprintf("untreated runs: %d units; absent for %d units (%s), all with replicate-1 item absent: %s\n",
            length(OFFU), nrow(U) - length(OFFU), paste(setdiff(U$unit, OFFU), collapse = ","),
            all(setdiff(U$unit, OFFU) %in% MISS[rep == 1, unit])))
print(U[!unit %in% OFFU, .N, by = .(branch, rescue_mode)])
CG = B[theta_treg == 1, .(n = .N), by = .(unit, phi_instability, treg_eff)][, .(cells = .N), by = unit]
cat(sprintf("units with all 121 cells: %d of %d\n", CG[cells == 121, .N], nrow(CG)))
cat("tail-mean bounds in the stored runs (all must be inside the read filter):\n")
print(B[, .(max_Eh = max(Eh_mean), max_Ei = max(Ei_mean), max_M0 = max(M0_mean), max_M1 = max(M1_mean),
            max_Tr = max(Tr_mean), max_PLP = max(PLP_mean), max_CLP = max(CLP_mean))])
cat("solve seconds per run by eta (treated):\n")
print(B[theta_treg == 1, .(med_s = median(solve_secs), q90 = quantile(solve_secs, .9)), by = treg_eff][order(treg_eff)])
cat(sprintf("total core-hours in stored runs: %.0f\n", sum(B$solve_secs) / 3600))

## ---- 4. cells file consistency -------------------------------------------
sec("4. cells file")
C = as.data.table(readRDS("./data/fig05_cells.rds"))
cat(sprintf("cells %d units %d\n", nrow(C), uniqueN(C$unit)))
print(merge(unique(C[, .(unit)]), U, by = "unit")[, .N, by = .(branch, rescue_mode)])
cat("res_n (replicates per cell) in cells file:\n"); print(table(C$res_n))

## ---- 5. replicate spread ---------------------------------------------------
sec("5. replicate spread")
T1 = B[theta_treg == 1 & unit %in% C$unit]
T1[, L := Ea_mean + En_mean]
T1[, res := sqrt((Ei_mean / 100)^2 + (PLP_mean / 100)^2)]
SP = T1[, .(n = .N, mL = mean(L), sL = sd(L), mr = mean(res), sr = sd(res),
            nopen = sum(L < Kg)), by = .(unit, phi_instability, treg_eff)]
cat(sprintf("eta = 1 cells: max across-replicate sd of L %.3g (deterministic)\n", SP[treg_eff == 1, max(sL)]))
X = SP[treg_eff < 1 & n > 1]
cat(sprintf("eta < 1 cells: %d; CV of res across replicates median %.4g q95 %.4g max %.4g\n", nrow(X),
            median(X$sr / X$mr), quantile(X$sr / X$mr, .95), max(X$sr / X$mr)))
cat(sprintf("eta < 1 cells: CV of L across replicates median %.4g q95 %.4g max %.4g\n",
            median(X$sL / X$mL), quantile(X$sL / X$mL, .95), max(X$sL / X$mL)))
cat(sprintf("split cells (0 < n_open < n): %d of %d (%.2f%%); at eta = 1: %d\n",
            SP[nopen > 0 & nopen < n, .N], nrow(SP), 100 * SP[nopen > 0 & nopen < n, .N] / nrow(SP),
            SP[treg_eff == 1 & nopen > 0 & nopen < n, .N]))
cat(sprintf("among cells NOT split, CV of L median %.4g q95 %.4g\n",
            median(X[nopen == 0 | nopen == n]$sL / X[nopen == 0 | nopen == n]$mL),
            quantile(X[nopen == 0 | nopen == n]$sL / X[nopen == 0 | nopen == n]$mL, .95)))
cat("split cells by eta:\n"); print(SP[, .(cells = .N, split = sum(nopen > 0 & nopen < n)), by = treg_eff][order(treg_eff)])
cat("within-run tail CV of Ei (sqrt(Ei_var)/Ei_mean), treated runs:\n")
print(T1[, as.list(quantile(sqrt(Ei_var) / pmax(Ei_mean, .Machine$double.eps), c(.5, .95, .99, 1)))])

## ---- 6. rescue outcome -----------------------------------------------------
sec("6. rescue outcome")
PC = as.data.table(readRDS("./data/fig05_switch_prob_cells.rds"))
PU = as.data.table(readRDS("./data/fig05_switch_prob_units.rds"))
PU = merge(PU, U[, .(unit, branch)], by = "unit")
cat("n_rep per cell:\n"); print(table(PC$n_rep))
print(PU[, .(units = .N, never = sum(p_mean == 0), always = sum(p_mean == 1),
             crossing = sum(p_mean > 0 & p_mean < 1)), by = .(branch, rescue_mode)][order(branch, rescue_mode)])
print(PU[, .(units = .N, never = sum(p_mean == 0), always = sum(p_mean == 1),
             crossing = sum(p_mean > 0 & p_mean < 1)), by = branch])
OFF = B[theta_treg == 0 & unit %in% PU$unit, .(unit, Kg, L_off = Ea_mean + En_mean, PLP_off = PLP_mean, branch)]
cat(sprintf("untreated arm rescued (L < Kg): %d of %d units\n", OFF[L_off < Kg, .N], nrow(OFF)))
cat("untreated L/Kg quantiles:\n"); print(quantile(OFF$L_off / OFF$Kg, c(0, .05, .5, .95, 1)))
PU2 = merge(PU, PC[abs(phi) < 1e-9 & abs(eta - 1) < 1e-9, .(unit, p_perfect = p_open)], by = "unit")
cat("perfect Treg (phi = 0, eta = 1): units rescued\n"); print(PU2[, .(units = .N, rescued = sum(p_perfect == 1), not = sum(p_perfect == 0)), by = branch])
PU3 = merge(PU, PC[abs(phi - 1) < 1e-9 & abs(eta) < 1e-9, .(unit, p_worst = p_open)], by = "unit")
cat("worst Treg (phi = 1, eta = 0): mean p by branch, all units\n"); print(PU3[, .(mean_p = mean(p_worst)), by = branch])
cat("crossing units not rescued at the perfect Treg:\n"); print(PU2[p_mean > 0 & p_mean < 1 & p_perfect < 1, .N, by = branch])
## pathogen in rescued runs (the criterion ignores P_LP)
RR = T1[L < Kg]
cat("tail-mean PLP among rescued treated runs, pathogenic_help:\n")
print(RR[branch == "pathogenic_help", as.list(quantile(PLP_mean, c(.5, .9, .99, 1)))])
cat(sprintf("rescued help runs with tail PLP > 1: %d of %d\n", RR[branch == "pathogenic_help" & PLP_mean > 1, .N],
            RR[branch == "pathogenic_help", .N]))
cat("tail L/Kg in rescued runs and non-rescued runs (quantiles):\n")
print(quantile(T1[L < Kg]$L / T1[L < Kg]$Kg, c(0, .5, .95, .99, 1)))
print(quantile(T1[L >= Kg]$L / T1[L >= Kg]$Kg, c(0, .01, .05, .5, 1)))

## ---- 7. the plotted surfaces ----------------------------------------------
sec("7. plotted quantities")
for (br in c("", "sterile", "pathogenic_help")) {
  su = if (nzchar(br)) U[branch == br, unit] else U$unit
  cr = PU[unit %in% su & p_mean > 0 & p_mean < 1, unit]
  H = PC[unit %in% cr, .(p = mean(p_open)), by = .(stab = round(stab, 1), eta = round(eta, 1))]
  cat(sprintf("\n--- %s: %d crossing units ---\n", if (nzchar(br)) br else "pooled", length(cr)))
  cat("P(rescue) surface (rows 1-phi from 1.0 down, columns eta 0..1):\n")
  W = dcast(H, stab ~ eta, value.var = "p"); setorder(W, -stab); print(W, digits = 3)
  setorder(H, stab, eta)
  DE = H[, .(eta = head(eta, -1) + 0.05, d = diff(p)), by = stab]
  DS = H[, .(stab = head(stab, -1) + 0.05, d = diff(p)), by = eta]
  cat(sprintf("steps along eta: min %.4f median %.4f max %.4f | along stab: min %.4f median %.4f max %.4f | negative steps %d / %d\n",
              min(DE$d), median(DE$d), max(DE$d), min(DS$d), median(DS$d), max(DS$d), sum(DE$d < 0), sum(DS$d < 0)))
  cat(sprintf("stab step 0.9 -> 1.0 by eta: %s\n", paste(sprintf("%.3f", DS[abs(stab - 0.95) < 1e-6][order(eta)]$d), collapse = " ")))
  cat(sprintf("stab step 0.9 -> 1.0 share of the total rise 0 -> 1, by eta: %s\n",
              paste(sprintf("%.2f", DS[abs(stab - 0.95) < 1e-6][order(eta)]$d /
                              (H[stab == 1][order(eta)]$p - H[stab == 0][order(eta)]$p)), collapse = " ")))
  CC = PC[unit %in% cr]
  MS = CC[, .(p = mean(p_open)), by = .(unit, x = round(stab, 1))]
  ME = CC[, .(p = mean(p_open)), by = .(unit, x = round(eta, 1))]
  q = function(M) M[, .(mean = mean(p), q05 = quantile(p, .05), q25 = quantile(p, .25), med = median(p),
                        q75 = quantile(p, .75), q95 = quantile(p, .95)), by = x][order(x)]
  cat("marginal along 1 - phi (per-unit mean over eta):\n"); print(q(MS), digits = 3)
  cat("marginal along eta (per-unit mean over 1 - phi):\n"); print(q(ME), digits = 3)
  ## trapezoid area of the mean curve minus the identity's 0.5
  ar = function(Q) { Q = Q[order(x)]; sum(diff(Q$x) * (head(Q$mean, -1) + tail(Q$mean, -1)) / 2) }
  cat(sprintf("area under mean curve: stability %.4f | eta %.4f (identity 0.5)\n", ar(q(MS)), ar(q(ME))))
}
