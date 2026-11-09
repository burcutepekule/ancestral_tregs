SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## census_core.R -- every count the Methods (docs/methods.pdf)
## quotes, recomputed from the pipeline products. Integrates nothing.
## ----------------------------------------------------------------------------
## Reads the assembled sweep batches (main and extension, trusted and untrusted),
## the tau grids, and the selection products; prints each count under a tag that
## core.tex cites in its "% from:" comments.
##
## The assembled main-sweep batches are ~750 MB each and are not copied into
## ./data; FP_DIR points at the folder that holds them (default ./data). Only the
## scalar columns are kept after reading, one batch at a time.
##
## in : <FP_DIR>/fp_big_batch_<b>.rds (+ _untrusted); ./data/tau_grid_v6_*.rds;
##      SEL_BIG (default ./data/selected_big.rds)
## out: console only
## RUN: FP_DIR=<folder with fp_big_batch_*.rds> Rscript pipeline/census_core.R
## ============================================================================
suppressMessages(library(data.table))
FP_DIR  = Sys.getenv("FP_DIR", "./data")
SEL_BIG = Sys.getenv("SEL_BIG", "./data/selected_big.rds")
GAP_CUT = as.numeric(Sys.getenv("GAP_CUT", "20"))
EFFECT_TOL = 1
STEP = 10 / 2^7                                   # bisection lattice, PROBE_BY / 2^7
pc = function(a, b) sprintf("%s of %s (%.2f%%)", format(a, big.mark = ","), format(b, big.mark = ","), 100 * a / b)
hdr = function(x) cat(sprintf("\n==================== %s ====================\n", x))

KEEP = c("tau_id","batch","row_ind","Kg","beta_ros","r_L_P","hi_lo_gap_off","n_cross_off","n_cross_on",
  "dir_off","dir_on","n_ray_off","n_ray_on","note_off","note_on","lo_off","hi_off","lo_on","hi_on",
  "lo_off_Eh","lo_off_PLP","lo_on_Eh","lo_on_PLP","hi_off_Eh","hi_off_PLP","hi_on_Eh","hi_on_PLP",
  "rescue_gap","rescue_gap_kind","Eh0_off","Eh0_on","n_field_stable_off","n_field_stable_on",
  "seed_mode_off","seed_mode_on","dev_med_off","dev_med_on","n_fail_off","n_fail_on",
  "n_unmatched_off","n_unmatched_on","map_trusted")
rd = function(f) { X = readRDS(f); setDT(X); Y = X[, intersect(KEEP, names(X)), with = FALSE]; rm(X); gc(); Y }
load_set = function(files) rbindlist(lapply(files, function(f) {
  Y = rd(f); Y[, trusted := !grepl("_untrusted\\.rds$", f)]; Y }), fill = TRUE)

## ---- [grid] the parameter grids -------------------------------------------------
hdr("[grid] tau grids")
tf = sort(Sys.glob("./data/tau_grid_v6_*.rds"))
TG = rbindlist(lapply(tf, readRDS), fill = TRUE)
print(TG[, .(rows = .N, tau_id_min = min(tau_id), tau_id_max = max(tau_id)), by = batch])
cat("kinetic parameter ranges over all grids (min / max):\n")
KIN = c("kappa_E","tau_r","tau_P_B","mu_Ea","mu_R","tau_01","tau_10","tau_1_d","tau_eff",
        "tau_t_a","tau_t_d","alpha_0","alpha_1_e","alpha_1_r","alpha_AMP","gamma_B_P",
        "gamma_B_C","sigma_D","sigma_S","tau_D","tau_S","w_E_over_M","phi_instability","treg_eff")
print(data.table(par = KIN, min = sapply(KIN, function(v) signif(min(TG[[v]]), 4)),
                 max = sapply(KIN, function(v) signif(max(TG[[v]]), 4))))

## ---- [sweep] what came back -------------------------------------------------------
hdr("[sweep] main sweep, assembled batches")
mf = sort(Sys.glob(file.path(FP_DIR, "fp_big_batch_*.rds")))
if (!length(mf)) stop("no fp_big_batch_*.rds under FP_DIR = ", FP_DIR)
A = load_set(mf)
A[, patho := r_L_P > 0]
cat(sprintf("combinations returned: %s (trusted %s, untrusted %s)\n", format(nrow(A), big.mark = ","),
            format(A[trusted == TRUE, .N], big.mark = ","), format(A[trusted == FALSE, .N], big.mark = ",")))
cat(sprintf("key unique: %s\n", uniqueN(A, by = c("tau_id","Kg","beta_ros","r_L_P")) == nrow(A)))
sw_b = sort(unique(A$batch))
cat("batches swept:", sw_b, "\n")
per_row = A[, .N, by = .(batch, row_ind)]
cat("combinations per returned tau row:\n"); print(per_row[, .N, by = .(combos = N)])
miss = TG[batch %in% sw_b][!per_row, on = .(batch, row_ind)][, .(batch, row_ind)]
cat(sprintf("tau rows in the swept grids: %d | returned: %d | missing: %d\n",
            TG[batch %in% sw_b, .N], nrow(per_row), nrow(miss)))
print(miss)
cat(sprintf("expected combinations %d x 441 = %s; returned %s\n", TG[batch %in% sw_b, .N],
            format(TG[batch %in% sw_b, .N] * 441, big.mark = ","), format(nrow(A), big.mark = ",")))
cat("levels: Kg", sort(unique(A$Kg)), "| beta_ros", sort(unique(A$beta_ros)), "| r_L_P", sort(unique(A$r_L_P)), "\n")

hdr("[trust] convergence exclusion (median relative deviation <= 1e-2, both arms)")
cat("untrusted:", pc(A[trusted == FALSE, .N], nrow(A)), "\n")
cat("check: untrusted <=> !(dev_med_off <= 1e-2 & dev_med_on <= 1e-2):",
    A[, all(trusted == (is.finite(dev_med_off) & dev_med_off <= 1e-2 & is.finite(dev_med_on) & dev_med_on <= 1e-2))], "\n")
cat("untrusted by arm: OFF fails", A[trusted == FALSE & !(is.finite(dev_med_off) & dev_med_off <= 1e-2), .N],
    "| ON fails", A[trusted == FALSE & !(is.finite(dev_med_on) & dev_med_on <= 1e-2), .N],
    "| dev_med NA (no stable assignment) OFF", A[trusted == FALSE & !is.finite(dev_med_off), .N],
    "ON", A[trusted == FALSE & !is.finite(dev_med_on), .N], "\n")
for (v in c("r_L_P", "Kg", "beta_ros"))
  print(A[, .(N = .N, excluded = sum(!trusted), pct = round(100 * mean(!trusted), 2)), by = v][order(get(v))])
T = A[trusted == TRUE]
cat(sprintf("retained %s: sterile %s, pathogenic %s\n", format(nrow(T), big.mark = ","),
            format(T[patho == FALSE, .N], big.mark = ","), format(T[patho == TRUE, .N], big.mark = ",")))

hdr("[seed] how the endpoint-seeded fixed-point set was built (retained combinations)")
print(dcast(T[, .N, by = .(patho, seed_mode_off)], seed_mode_off ~ patho, value.var = "N", fill = 0))
print(dcast(T[, .(N = .N), by = .(patho, seed_mode_on)], seed_mode_on ~ patho, value.var = "N", fill = 0))
cat("integrations failed (n_fail) > 0: OFF", T[n_fail_off > 0, .N], "ON", T[n_fail_on > 0, .N], "\n")
cat("probed Eh0 unassigned (n_unmatched) > 0: OFF", T[n_unmatched_off > 0, .N], "ON", T[n_unmatched_on > 0, .N], "\n")
cat("n_ray_off distribution:\n"); print(T[, .N, by = .(patho, n_ray_off)][order(patho, n_ray_off)])
cat("note_off:\n"); print(dcast(T[, .N, by = .(patho, note_off)], note_off ~ patho, value.var = "N", fill = 0))
cat("n_ray == 1 <=> n_cross == 0 violations: OFF", T[(n_ray_off == 1L) != (n_cross_off == 0L), .N],
    "ON", T[(n_ray_on == 1L) != (n_cross_on == 0L), .N], "\n")

## ---- [gate] ----------------------------------------------------------------------
hdr("[gate] the bistability gate, Tregs OFF, per base branch")
T[, `:=`(g1 = is.finite(hi_lo_gap_off) & hi_lo_gap_off > GAP_CUT,
         g2 = is.finite(n_cross_off) & n_cross_off == 1L,
         g3 = !is.na(dir_off) & dir_off == "increasing")]
G = T[, .(start = .N, gap_gt_20 = sum(g1), plus_one_crossing = sum(g1 & g2), plus_increasing = sum(g1 & g2 & g3)), by = patho]
print(G)
cat("hi_lo_gap_off finite (two attractors reached, Tregs OFF):\n"); print(T[, .(finite = sum(is.finite(hi_lo_gap_off)), n_ray_off_ge2 = sum(n_ray_off >= 2)), by = patho])
cat("removed by (2), by n_cross_off:\n"); print(T[g1 & !g2, .N, by = .(patho, n_cross_off)][order(patho, n_cross_off)])
cat("removed by (3), by dir_off:\n"); print(T[g1 & g2 & !g3, .N, by = .(patho, dir_off)])
cat("pathogenic removed by (3), note_off:\n"); print(T[g1 & g2 & !g3, .N, by = .(patho, note_off)])
cat("pathogenic decreasing by r_L_P:\n"); print(T[g1 & g2 & patho == TRUE, .(n = .N, decreasing = sum(dir_off == "decreasing", na.rm = TRUE)), by = r_L_P][order(r_L_P)])
SEL = T[g1 & g2 & g3]
## the same gate steps split by the Treg effect each combination would be assigned
## (computable wherever both low-burden attractors exist; the branches themselves
## are formed only after the gate)
T[, eff := { d = sqrt((100 - lo_on_Eh)^2 + lo_on_PLP^2) - sqrt((100 - lo_off_Eh)^2 + lo_off_PLP^2)
             fifelse(!is.finite(d), "unknown", fifelse(d < -1, "help", fifelse(d > 1, "harm", "neutral"))) }]
cat("gate steps by base branch x Treg effect:\n")
print(T[, .(start = .N, gap_gt_20 = sum(g1), plus_one_crossing = sum(g1 & g2), plus_increasing = sum(g1 & g2 & g3)),
        by = .(patho, eff)][order(patho, eff)])
cat(sprintf("gate yields %s\n", format(nrow(SEL), big.mark = ",")))
## the Chebyshev separation the gate does not use, for the record
SELc = T[is.finite(hi_lo_gap_off)]
SELc[, cheb := pmax(abs(hi_off_Eh - lo_off_Eh), abs(hi_off_PLP - lo_off_PLP))]
cat("composite vs Chebyshev separation at 20 (rows with two OFF attractors): agree",
    pc(SELc[, sum((hi_lo_gap_off > GAP_CUT) == (cheb > GAP_CUT))], nrow(SELc)), "\n")
print(SELc[, .(n = .N, comp_gt20 = sum(hi_lo_gap_off > GAP_CUT), cheb_gt20 = sum(cheb > GAP_CUT),
               comp_only = sum(hi_lo_gap_off > GAP_CUT & !(cheb > GAP_CUT)), cheb_only = sum(cheb > GAP_CUT & !(hi_lo_gap_off > GAP_CUT))), by = patho])
cat("gate (1) per beta_ros, pathogenic, pass fraction:\n")
print(T[patho == TRUE, .(N = .N, pass = sum(g1 & g2 & g3), pct = round(100 * mean(g1 & g2 & g3), 2)), by = beta_ros][order(beta_ros)])
cat("gate per beta_ros, sterile, pass fraction:\n")
print(T[patho == FALSE, .(N = .N, pass = sum(g1 & g2 & g3), pct = round(100 * mean(g1 & g2 & g3), 2)), by = beta_ros][order(beta_ros)])
cat("gate per Kg:\n")
print(T[, .(N = .N, pass = sum(g1 & g2 & g3), pct = round(100 * mean(g1 & g2 & g3), 2)), by = .(patho, Kg)][order(patho, Kg)])
cat("gate per r_L_P:\n")
print(T[, .(N = .N, pass = sum(g1 & g2 & g3), pct = round(100 * mean(g1 & g2 & g3), 2)), by = r_L_P][order(r_L_P)])
cat("distinct tau rows in the gated set: sterile", SEL[patho == FALSE, uniqueN(tau_id)], "pathogenic", SEL[patho == TRUE, uniqueN(tau_id)], "all", uniqueN(SEL$tau_id), "\n")

## ---- [class] rescue_mode and treg_effect, exactly as 03_select.R --------------------
hdr("[class] what Tregs do, on the gated set")
SEL[, rescue_mode := fifelse(!is.finite(n_ray_on), "unresolved_on",
                     fifelse(n_ray_on >= 2L, "bistable",
                      fifelse(note_on == "lo_basin_throughout", "mono_B",
                       fifelse(is.finite(lo_on) & is.finite(lo_off) & is.finite(hi_off) &
                               abs(lo_on - lo_off) <= abs(lo_on - hi_off), "mono_B", "mono_A"))))]
cat("how the mono_B label was reached: note_on == lo_basin_throughout",
    SEL[is.finite(n_ray_on) & n_ray_on < 2 & note_on == "lo_basin_throughout", .N],
    "| by the composite rule", SEL[is.finite(n_ray_on) & n_ray_on < 2 & note_on != "lo_basin_throughout" & rescue_mode == "mono_B", .N], "\n")
cat("note_on among n_ray_on == 1:\n"); print(SEL[n_ray_on == 1, .N, by = .(patho, note_on)])
cat("n_ray_on distribution:\n"); print(SEL[, .N, by = .(patho, n_ray_on)][order(patho, n_ray_on)])
cmpf = function(Eh, PLP) sqrt((100 - Eh)^2 + PLP^2)
SEL[, d_lo_sel := cmpf(lo_on_Eh, lo_on_PLP) - cmpf(lo_off_Eh, lo_off_PLP)]
SEL[, treg_effect := fifelse(!is.finite(d_lo_sel), "unknown",
                      fifelse(d_lo_sel < -EFFECT_TOL, "help", fifelse(d_lo_sel > EFFECT_TOL, "harm", "neutral")))]
print(dcast(SEL[, .N, by = .(patho, rescue_mode, treg_effect)], patho + rescue_mode ~ treg_effect, value.var = "N", fill = 0))
cat("lo_on vs component composite disagreement > 1e-6:", SEL[is.finite(lo_on) & abs(lo_on - cmpf(lo_on_Eh, lo_on_PLP)) > 1e-6, .N], "\n")

## the three branches
SEL[, branch := fifelse(!patho, "sterile", fifelse(treg_effect == "help", "pathogenic_help",
                 fifelse(treg_effect == "harm", "pathogenic_harm", paste0("pathogenic_", treg_effect))))]
cat("\nbranches:\n"); print(SEL[, .N, by = branch])

hdr("[gap] rescue gap per branch and class")
for (br in c("sterile", "pathogenic_help", "pathogenic_harm")) {
  B = SEL[branch == br]
  cat(sprintf("\n-- %s (n = %s) --\n", br, format(nrow(B), big.mark = ",")))
  print(B[, .N, by = .(rescue_mode, rescue_gap_kind, dir_on)][order(rescue_mode, rescue_gap_kind)])
  m = B[rescue_mode == "bistable" & rescue_gap_kind == "measured" & is.finite(rescue_gap), rescue_gap]
  cat(sprintf("bistable measured gaps: %d | exactly 0: %.1f%% | <= one step (%.6f): %.1f%% | > 0: %.1f%% | < 0: %d\n",
              length(m), 100 * mean(m == 0), STEP, 100 * mean(abs(m) <= STEP + 1e-12), 100 * mean(m > 0), sum(m < 0)))
  if (length(m)) { cat("quantiles 50/75/90/95/99/max: "); cat(round(quantile(m, c(.5, .75, .9, .95, .99, 1)), 4), "\n")
    cat("on lattice (multiple of step within 1e-9, or of 9/128 in [90,99]):",
        mean(abs(m / STEP - round(m / STEP)) < 1e-6), "\n") }
  bi_fig = B[rescue_mode == "bistable" & is.finite(d_lo_sel) & is.finite(hi_on) & is.finite(rescue_gap)]
  cat("bistable with finite d_hi and finite rescue_gap (S-figure population B/C):", nrow(bi_fig), "\n")
  mb = B[rescue_mode == "mono_B"]
  cat("mono_B:", nrow(mb), "| finite lo_on & hi_off:", mb[is.finite(lo_on) & is.finite(hi_off), .N],
      "| gap kind:", paste(mb[, .N, by = rescue_gap_kind][, sprintf("%s %d", rescue_gap_kind, N)], collapse = "; "),
      "| second stable ON fixed point off the line (n_field_stable_on >= 2):", mb[n_field_stable_on >= 2, .N], "\n")
  ma = B[rescue_mode == "mono_A"]
  cat("mono_A:", nrow(ma), "| finite lo_on & hi_off:", ma[is.finite(lo_on) & is.finite(hi_off), .N], "\n")
  cat("Eh0_off on the line: median", round(median(B$Eh0_off), 3), "range", round(range(B$Eh0_off), 3), "\n")
}

## ---- [check] the stored selection is what this recomputation gives ------------------
hdr("[check] stored selection products")
if (file.exists(SEL_BIG)) {
  S0 = as.data.table(readRDS(SEL_BIG))
  CK = c("tau_id", "Kg", "beta_ros", "r_L_P")
  cat(SEL_BIG, nrow(S0), "rows | same keys as recomputed gate:",
      nrow(S0) == nrow(SEL) && nrow(merge(S0[, ..CK], SEL[, ..CK], by = CK)) == nrow(SEL), "\n")
  m = merge(S0[, c(CK, "rescue_mode", "treg_effect"), with = FALSE], SEL[, c(CK, "rescue_mode", "treg_effect"), with = FALSE], by = CK)
  cat("rescue_mode identical:", m[, all(rescue_mode.x == rescue_mode.y)], "| treg_effect identical:", m[, all(treg_effect.x == treg_effect.y)], "\n")
}
for (br in c("sterile", "pathogenic_help", "pathogenic_harm")) for (k in c("skeleton", "rescue")) {
  f = sprintf("./data/%s_%s.rds", k, br)
  if (file.exists(f)) cat(f, nrow(readRDS(f)), "rows (recomputed branch:", SEL[branch == br, .N], ")\n")
}

## ---- [skel] the numbers behind Figs. 3A, 3B, Extended Data Figs. 1, 2A, 2B, Supplementary Fig. 2 ------------------
hdr("[skel] populations, attractor shifts and rescue gaps drawn in Fig. 3 / Extended Data Figs. 1, 2 / Supplementary Fig. 2")
CKs = c("tau_id","Kg","beta_ros","r_L_P")
q3 = function(x) { x = x[is.finite(x)]; sprintf("%.2f [%.2f, %.2f] (min %.2f, max %.2f)", median(x), quantile(x, .25), quantile(x, .75), min(x), max(x)) }
for (br in c("sterile", "pathogenic_help", "pathogenic_harm")) {
  f = sprintf("./data/skeleton_%s.rds", br); if (!file.exists(f)) next
  K = as.data.table(readRDS(f))
  BI = K[rescue_mode == "bistable" & is.finite(d_hi) & is.finite(rescue_gap)]
  MB = K[rescue_mode == "mono_B" & is.finite(lo_on) & is.finite(hi_off)]
  cat(sprintf("\n-- %s: skeleton rows %d | bistable %d -> drawn %d | ON_mono %d -> drawn %d | mono_A %d\n", br, nrow(K),
              K[rescue_mode == "bistable", .N], nrow(BI), K[rescue_mode == "mono_B", .N], nrow(MB), K[rescue_mode == "mono_A", .N]))
  cat("  bistable dropped: d_hi not finite", K[rescue_mode == "bistable" & !is.finite(d_hi), .N],
      "| rescue_gap not finite", K[rescue_mode == "bistable" & is.finite(d_hi) & !is.finite(rescue_gap), .N], "\n")
  sh = list(`comp L` = BI$d_lo, `comp H` = BI$d_hi, `dEh L` = BI$d_lo_Eh, `dEh H` = BI$d_hi_Eh, `dPLP L` = BI$d_lo_PLP, `dPLP H` = BI$d_hi_PLP)
  for (k in names(sh)) cat(sprintf("  B %-7s median [IQR] %s | <0 %.1f%% | >0 %.1f%%\n", k, q3(sh[[k]]),
                                   100 * mean(sh[[k]] < 0, na.rm = TRUE), 100 * mean(sh[[k]] > 0, na.rm = TRUE)))
  sh = list(`comp L` = MB$lo_on - MB$lo_off, `comp H` = MB$lo_on - MB$hi_off, `dEh L` = MB$lo_Eh_on - MB$lo_Eh_off,
            `dEh H` = MB$lo_Eh_on - MB$hi_Eh_off, `dPLP L` = MB$lo_PLP_on - MB$lo_PLP_off, `dPLP H` = MB$lo_PLP_on - MB$hi_PLP_off)
  for (k in names(sh)) cat(sprintf("  E %-7s median [IQR] %s | <0 %.1f%% | >0 %.1f%%\n", k, q3(sh[[k]]),
                                   100 * mean(sh[[k]] < 0, na.rm = TRUE), 100 * mean(sh[[k]] > 0, na.rm = TRUE)))
  g = BI$rescue_gap
  cat(sprintf("  C gaps: n %d | exactly 0 %.1f%% | > 0 %.1f%% | < 0 %d | <= one step %.1f%% | max %.4f | median of >0 %.4f\n",
              length(g), 100 * mean(g == 0), 100 * mean(g > 0), sum(g < 0), 100 * mean(abs(g) <= STEP + 1e-12), max(g), median(g[g > 0])))
  tb = table(round(g / STEP) * STEP); cat("  C largest bars:", paste(sprintf("%s:%.1f%%", names(tb), 100 * as.numeric(tb) / length(g))[order(-tb)][1:4], collapse = " "), "\n")
  cat(sprintf("  F ON_mono: censored %d of %d (gap kind on_lo_basin_throughout); recorded extension Eh0_off - 0: median %.2f, range %.2f-%.2f\n",
              MB[gap_censored %in% TRUE, .N], nrow(MB), median(MB$rescue_gap), min(MB$rescue_gap), max(MB$rescue_gap)))
  ex = BI[which.max(rescue_gap)]; em = MB[which.max(lo_off - lo_on)]
  for (e in list(list("A", ex), list("D", em))) { z = e[[2]]
    cat(sprintf("  %s example: tau_id %d (batch %d row %d) Kg %g beta_ros %g r_L_P %g | gap %.4f | Eh0_off %.4f Eh0_on %s\n", e[[1]], z$tau_id, z$batch, z$row_ind,
                z$Kg, z$beta_ros, z$r_L_P, z$rescue_gap, z$Eh0_off, format(z$Eh0_on)))
    cat(sprintf("     OFF lo (Eh %.2f, M1 %.2f, PLP %.3f, comp %.2f) hi (Eh %.2f, M1 %.2f, PLP %.3f, comp %.2f)\n",
                z$lo_Eh_off, z$lo_M1_off, z$lo_PLP_off, z$lo_off, z$hi_Eh_off, z$hi_M1_off, z$hi_PLP_off, z$hi_off))
    cat(sprintf("     ON  lo (Eh %.2f, M1 %.2f, PLP %.3f, comp %.2f) hi (Eh %s, M1 %s, PLP %s, comp %s)\n",
                z$lo_Eh_on, z$lo_M1_on, z$lo_PLP_on, z$lo_on, format(round(z$hi_Eh_on, 2)), format(round(z$hi_M1_on, 2)),
                format(round(z$hi_PLP_on, 3)), format(round(z$hi_on, 2)))) }
}
cat("\n[census] done\n")
