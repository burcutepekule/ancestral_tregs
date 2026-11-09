## ============================================================================
## pipeline/checks/harm_methods_numbers.R -- every count and summary statistic
## quoted in docs/figure_methods.pdf that is not printed by
## figure_scripts/Supplementary_Fig_3.R or by the two harm_check_* scripts.
## Reads only ./data; integrates nothing.
## RUN: Rscript pipeline/checks/harm_methods_numbers.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages(library(data.table))
KK = c("batch", "row_ind", "Kg", "beta_ros", "r_L_P")
q3 = function(x) sprintf("median %.4g (IQR %.4g-%.4g; range %.4g-%.4g; n=%d)",
                         median(x), quantile(x, .25), quantile(x, .75), min(x), max(x), length(x))

cat("== 1. the gate on skeleton_pathogenic_harm.rds\n")
S = as.data.table(readRDS("./data/skeleton_pathogenic_harm.rds"))
cat("start", nrow(S), "| rescue_mode:", paste(names(table(S$rescue_mode)), table(S$rescue_mode), collapse = " "), "\n")
S1 = S[is.finite(Eh0_off) & Eh0_off < 95]; cat("Eh0_off < 95:", nrow(S1), "\n")
S2 = S1[is.na(dir_on) | dir_on != "decreasing"]; cat("dir_on != decreasing:", nrow(S2), " (removed", nrow(S1) - nrow(S2), ")\n")
B  = S2[rescue_mode == "bistable"]
cat("bistable after dir:", nrow(B), "| rescue_gap == 0:", sum(B$rescue_gap == 0, na.rm = TRUE),
    "| < 0:", sum(B$rescue_gap < 0, na.rm = TRUE), "| non-finite:", sum(!is.finite(B$rescue_gap)), "\n")
S3 = S2[rescue_mode != "bistable" | (is.finite(rescue_gap) & rescue_gap >= 1e-9)]
cat("GAP_MIN 1e-9:", nrow(S3), "| mode:", paste(names(table(S3$rescue_mode)), table(S3$rescue_mode), collapse = " "), "\n")
cat("tau rows", uniqueN(S3[, .(batch, row_ind)]), "| Kg", paste(sort(unique(S3$Kg)), collapse = ","),
    "| r_L_P", paste(sort(unique(S3$r_L_P)), collapse = ","), "| beta_ros", paste(sort(unique(S3$beta_ros)), collapse = ","), "\n")
cat("d_lo (sweep harm size):", q3(S3$d_lo), "\n")
cat("Eh0_off:", q3(S3$Eh0_off), "\n")

cat("\n== 2. the integrated file\n")
Z = readRDS("./data/harm_rescue_data.rds"); SM = as.data.table(Z$sm); TH = as.data.table(Z$thresh)
cat("sm rows", nrow(SM), "combinations", uniqueN(SM$cid), "arms per combination:",
    paste(unique(SM[, .N, by = cid]$N), collapse = ","), "\n")
cat("cids present: max", max(SM$cid), " missing from 1..1773:", paste(setdiff(1:1773, SM$cid), collapse = ","), "\n")
cat("output times", length(unique(round(Z$tr$time, 8))), " tr rows", nrow(Z$tr), "\n")
TOL = Z$settings$EFFECT_TOL
SM[, effect := fifelse(dcomp < -TOL, "help", fifelse(dcomp > TOL, "harm", "neutral"))]
u = SM[treat == "delay" & level == 0]
cat("untreated class:", paste(names(table(u$effect)), table(u$effect), collapse = " "), "\n")
a = merge(u[, c(KK, "dcomp"), with = FALSE], SM[treat == "strength" & level == 1, c(KK, "dcomp"), with = FALSE], by = KK)
b = merge(u[, c(KK, "dcomp"), with = FALSE], SM[treat == "deact" & level == 1, c(KK, "dcomp"), with = FALSE], by = KK)
d = merge(u[, c(KK, "dcomp"), with = FALSE], SM[treat == "induce" & level == 1, c(KK, "dcomp"), with = FALSE], by = KK)
cat("untreated arm four ways, max |diff|:", max(abs(c(a$dcomp.x - a$dcomp.y, b$dcomp.x - b$dcomp.y, d$dcomp.x - d$dcomp.y))), "\n")
cat("untreated dcomp among harmed:", q3(u[effect == "harm", dcomp]), "\n")

cat("\n== 3. dip landmarks on the Tregs OFF arm\n")
O = SM[treat == "off"]; H = merge(O, u[effect == "harm", ..KK], by = KK)
for (nm in c("all 1771" , "harmed")) { X = if (nm == "harmed") H else O
  cat("--", nm, "n =", nrow(X), "\n")
  cat("t_PLP_peak", q3(X$t_PLP_peak), "\nPLP_peak", q3(X$PLP_peak), "\n")
  cat("t_dip", q3(X$t_dip), "\nPLP_dip", q3(X$PLP_dip), "\n")
  cat("t_ctrl defined", sum(!is.na(X$t_ctrl)), ":", q3(X$t_ctrl[!is.na(X$t_ctrl)]), "\n")
  cat("t_zero defined", sum(!is.na(X$t_zero)), ":", q3(X$t_zero[!is.na(X$t_zero)]), "\n")
  cat("t_dip == t_zero where defined:", X[!is.na(t_zero), all(t_dip == t_zero)], "\n")
  cat("PLP_dip == 0 count", sum(X$PLP_dip == 0), "; PL_end<EXT & PLP_end<EXT", sum(X$PL_end < 1e-11 & X$PLP_end < 1e-11), "\n")
  cat("t_dip == last output time (25000) count", sum(X$t_dip >= 25000 - 1e-6), "\n")
  cat("PLP_off (end)", q3(X$PLP_off), "\nEh_off", q3(X$Eh_off), "\n")
}

cat("\n== 4. thresholds\n")
print(TH[, .(n = .N, nonmono = sum(nonmono), help = sum(class_at_thr == "Harm to help", na.rm = TRUE),
             neutral = sum(class_at_thr == "Harm to neutral", na.rm = TRUE),
             harm = sum(class_at_thr == "Harm to harm", na.rm = TRUE), na = sum(is.na(thr))), by = .(treat, bound)][order(treat, bound)])
cat("Treat with no bracket (still harmed at the far end of the ladder): per treatment, harmed untreated 1745 minus neutral brackets\n")
TH[, fold := fifelse(treat %in% c("strength", "induce"), 1 / thr, thr)]
print(TH[, .(n = .N, frac_of_1745 = round(.N / 1745, 4),
             med = signif(median(fold), 4), q25 = signif(quantile(fold, .25), 4),
             q75 = signif(quantile(fold, .75), 4), min = signif(min(fold), 4), max = signif(max(fold), 4)),
         by = .(treat, bound)][order(treat, bound)])
W = dcast(TH, batch + row_ind + Kg + beta_ros + r_L_P + treat ~ bound, value.var = "thr")
print(W[is.finite(neutral) & is.finite(help), .(n = .N, same = sum(abs(log(help / neutral)) < 1e-6),
        sep_med_log10 = round(median(abs(log10(help / neutral))), 4)), by = treat])
## bracket ends: the ladder levels a and b
print(TH[treat == "delay" & bound == "neutral", .N, by = .(a, b)][order(a)])
cat("delay d_star:", q3(TH[treat == "delay" & bound == "neutral", thr]), "\n")
cat("bracket widths respected: delay |thr-a| > 0 all:", TH[treat == "delay", all(thr > a & thr <= b)], "\n")

cat("\n== 5. help/harm explore cache\n")
E = readRDS("./data/harm_explore_cache.rds"); str(E$key); TRx = as.data.table(E$tr)
print(TRx[, .(combinations = uniqueN(id)), by = grp])
cat("output times", uniqueN(TRx$time), " range", range(TRx$time), "\n")
cat("at last time, PLP <= 1: OFF", TRx[theta == 0 & time == max(time), mean(PLP <= 1)],
    " ON", TRx[theta == 1 & time == max(time), mean(PLP <= 1)], "\n")
print(TRx[time == max(time), .(frac_extinct = mean(PLP == 0), med_PLP = median(PLP)), by = .(grp, theta)])
