## ============================================================================
## Fig_4.R -- Fig. 4: the
## damage-matched comparison of the Tregs OFF and Tregs ON systems. With a brake on
## inflammation the host can run a stronger ROS response at the same cumulative
## epithelial damage, and in some combinations that clears an infection which
## persists without Tregs.
## ----------------------------------------------------------------------------
## Reads the damage-matched experiment on the gated units (pipeline/43_damage_matched_
## gated_units.R, pipeline/40_damage_matched.R via pipeline/cluster/damage_matched_gated.sh;
## cost matched over the convergence window, MATCH = conv), the trajectories of 1,000 of
## its pairs (pipeline/41_damage_matched_trajectories.R) and the demonstration unit
## (lib/fig_demo_case.R). Pairs: the filter of lib/fig_damage_matched.R. Panels:
##   A  per reference level beta_r, the beta_r* Tregs ON needs for equal damage (median,
##      IQR over all pairs of the filter; label: median beta_r*/beta_r), against the diagonal.
##   B  the order of events across the trajectory pairs, time in units of each pair's
##      t_peak: P_LP of Tregs OFF and of Tregs ON at the same beta_r, T_r of
##      Tregs ON and the suppression M1(OFF) - M1(ON), each pair scaled to its own
##      range over the drawn span; median and IQR over pairs. NOTHING IN B IS A MAGNITUDE.
##   C  the same pairs: P_LP of Tregs ON at beta_r and at beta_r*, relative to Tregs OFF
##      at beta_r (median, IQR over pairs of the log ratio).
##   D  the demonstration unit (UNIT; default 38614, one of the 15 pairs in which only
##      Tregs-ON at beta_r* eliminates): M1, P_LP and E_h of the four arms
##      OFF(beta_r), ON(beta_r), ON(beta_r*), OFF(beta_r*).
## t_peak: the time of the first peak of the pathogen of Tregs OFF at beta_r -- the first
## time the load turns over -- read off each pair's own trajectory. It needs no threshold,
## and it puts that peak at 1 in B and C, so the suppression that follows it is placed
## against the event it follows rather than against a chosen divergence time.
##
## in : ./data/fig04_damage_matched_gated.rds, ./data/fig04_gated_units.rds,
##      ./data/fig04b_v3_trajectories.rds, ./lib/*
## out: ./figures/main/Fig_4.png
## RUN: Rscript figure_scripts/Fig_4.R
## ============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
source("./lib/fig_damage_matched.R")
BASE_SZ  = as.numeric(Sys.getenv("BASE_SIZE", "42"))
F_TRAJ   = Sys.getenv("F_TRAJ", "./data/fig04b_v3_trajectories.rds")
UNIT     = as.integer(Sys.getenv("UNIT", "38614"))
OUT_FIG  = Sys.getenv("OUT_FIG", "./figures/main/Fig_4.png")
G_N      = as.integer(Sys.getenv("G_N", "400"))

## one population for every panel: the pairs whose brake saves at least DAM_MIN of the
## cumulative damage (lib/fig_damage_matched.R). B and C use those of them whose
## trajectories were stored.
M = fig_dm_matched()

## ---- the stored trajectory pairs (pipeline/41_damage_matched_trajectories.R) ------------
## Only the stored pairs that the filter of lib/fig_damage_matched.R keeps are drawn; the
## file written by pipeline/41 with its defaults (N_TRAJ = 1000, SEED = 1) is drawn from that
## filter, so all of its pairs are used. The count is reported.
if (!file.exists(F_TRAJ)) stop(F_TRAJ, " not found -- run: Rscript pipeline/41_damage_matched_trajectories.R")
X = readRDS(F_TRAJ); S = X$units; TR = X$traj
n_stored = uniqueN(S$i); S = S[i %in% M$i]; TR = TR[i %in% S$i]
if (!nrow(S)) stop("no stored trajectory pair survives the damage floor -- re-run pipeline/41_damage_matched_trajectories.R")
TR = merge(TR, S[, .(i, b_off, b_star, Eh0, r_L_P)], by = "i")
## each pair's own clock: the time of the FIRST peak of the Tregs-OFF pathogen -- the first
## strict local maximum on the fine grid, the same quantity Extended Data Fig. 5 compares (pipeline/42_first_peak.R).
## Not the maximum over the run: for a pair that climbs to a high-burden attractor that maximum
## is the attractor, which would put the clock hundreds of time units after the acute phase.
TP = TR[run == "OFF" & time > 0, { k = which(diff(sign(diff(PLP))) < 0) + 1L
          .(t_peak = if (length(k)) time[k[1]] else time[which.max(PLP)]) }, by = i]
stopifnot(all(is.finite(TP$t_peak)), all(TP$t_peak > 0))
TR = merge(TR, TP, by = "i")
cat(sprintf("[Fig. 4] trajectories for %d of the %d stored pairs (the rest are outside the damage floor)\n",
            uniqueN(TR$i), n_stored))

## colours: teal for the Tregs-OFF pathogen (A), grey for the beta_r condition and pink
## for the beta_r* condition (B, C, D)
COL = c(OFF = "#2b6b7f", same = "#555555", star = "#D81B7A")
th = theme_minimal(base_size = BASE_SZ) +
  theme(panel.grid.minor = element_blank(), legend.title = element_blank(),
        plot.title = element_text(size = BASE_SZ, face = "bold"))
## the span every pair covers with real output on the t_peak clock, so nothing is extrapolated
G = exp(seq(log(0.1), log(100), length.out = G_N))
TB = scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = function(z) sapply(z, format),
                   minor_breaks = as.vector(outer(1:9, 10^(-2:2))))
T_PEAK_LINE = geom_vline(xintercept = 1, colour = "grey30", linetype = "dotted", linewidth = 1.1)

## ---- the order of events (drawn second, tagged B) ------------------------------------------
## EACH PAIR IS SCALED TO ITS OWN RANGE over the drawn time span, then the median and IQR
## are taken over pairs at each tau: the two pathogen arms together (so OFF and ON stay
## comparable within a pair), Tregs and suppression each on their own. Pairs differ in
## scale by orders of magnitude; scaling pooled quantiles would let a few large pairs set
## the range.
ser = TR[time > 0, {
  o = .SD[run == "OFF"]; a = .SD[run == "ON_same"]
  x = log(o$time / t_peak[1]); f = function(y) approx(x, y, xout = log(G), rule = 2)$y
  sc = function(y, lo, hi) (y - lo) / pmax(hi - lo, .Machine$double.eps)
  po = f(o$PLP); pa = f(a$PLP); tr = f(a$Tr); dm = f(o$M1 - a$M1)
  lo = min(po, pa); hi = max(po, pa)
  list(tau = G, PLP_off = sc(po, lo, hi), PLP_on = sc(pa, lo, hi), Tr = sc(tr, min(tr), max(tr)), dM1 = sc(dm, min(dm), max(dm)))
}, by = i]
Z = melt(ser, id.vars = c("i", "tau"), variable.name = "var")
Z = Z[, .(med = median(value), q1 = quantile(value, .25), q3 = quantile(value, .75)), by = .(var, tau)]
LV_A = c("PLP_off", "PLP_on", "Tr", "dM1"); Z[, var := factor(var, levels = LV_A)]
COL_A = setNames(c(COL[["OFF"]], COL[["OFF"]], "#009E73", "#E07B00"), LV_A)
LTY_A = setNames(c("dashed", "solid", "solid", "solid"), LV_A)
## every series in this panel is at beta_r, so the rate is stated in the caption, not four times in the key
LAB_A = expression(paste("Tregs-OFF, ", P[LP]), paste("Tregs-ON, ", P[LP]),
                   paste("Tregs-ON, ", T[r]), paste("Tregs (OFF " - "ON), ", M[1]))
pnA = ggplot(Z, aes(tau, med, colour = var, fill = var, linetype = var)) +
  T_PEAK_LINE +
  geom_ribbon(aes(ymin = q1, ymax = q3), colour = NA, alpha = 0.13, show.legend = FALSE) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = COL_A, labels = LAB_A) + scale_fill_manual(values = COL_A, labels = LAB_A) +
  scale_linetype_manual(values = LTY_A, labels = LAB_A) +
  TB + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression("Time / " * t[peak]), y = "Scaled to its own range") + th +
  theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.4),
        panel.grid.minor = element_line(colour = "grey93", linewidth = 0.25),
        ## the key sits above the panel, two rows of two: every corner inside the panel is taken
        legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, -BASE_SZ * 2.0),
        legend.box.spacing = unit(2, "pt"), legend.text = element_text(size = BASE_SZ * 0.76),
        legend.key.width = unit(BASE_SZ * 0.06, "cm"), legend.key.spacing.y = unit(0, "pt"),
        legend.key.spacing.x = unit(BASE_SZ * 0.2, "pt")) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE))

## ---- the headroom, per reference level (drawn first, tagged A) -----------------------------
B = M[, .(med = median(b_star), q1 = quantile(b_star, .25), q3 = quantile(b_star, .75),
          ratio = median(b_star / b_off), n = .N), by = b_off]
lim = range(c(B$b_off, B$q1, B$q3))
pnB = ggplot(B, aes(b_off, med)) +
  geom_abline(slope = 1, intercept = 0, colour = COL[["OFF"]], linetype = "dashed") +
  geom_linerange(aes(ymin = q1, ymax = q3), colour = COL[["star"]], linewidth = 1.2) +
  geom_point(colour = COL[["star"]], size = 3) +
  ## nine closely spaced levels: labels repelled along y in a column to the left, each joined to its point
  ggrepel::geom_text_repel(aes(label = sprintf("%.2f×", ratio)), nudge_x = -0.35, hjust = 1, direction = "y", size = BASE_SZ * 0.24,
                           colour = "grey25", segment.colour = "grey70", segment.size = 0.3, min.segment.length = 0, seed = 1) +
  scale_x_log10(limits = lim * c(0.3, 1.1), breaks = c(0.003, 0.03)) + scale_y_log10(limits = lim) +
  labs(x = expression(beta[r] * ", Tregs-OFF"), y = expression(beta[r]^"*" * ", Tregs-ON")) + th

## ---- C: the pathogen relative to Tregs OFF, over time -----------------------------------
rat = TR[time > 0, {
  o = .SD[run == "OFF"]; a = .SD[run == "ON_same"]; s = .SD[run == "ON_star"]
  x = log(o$time / t_peak[1])
  f = function(y) approx(x, log(pmax(y, 1e-12) / pmax(o$PLP, 1e-12)), xout = log(G), rule = 2)$y
  list(tau = G, same = f(a$PLP), star = f(s$PLP))
}, by = i]
C = melt(rat, id.vars = c("i", "tau"), measure.vars = c("same", "star"), variable.name = "run", value.name = "lr")
C = C[, .(med = exp(median(lr)), q1 = exp(quantile(lr, .25)), q3 = exp(quantile(lr, .75))), by = .(run, tau)]
C[, run := factor(run, levels = c("same", "star"))]
LAB_C = c(expression("Tregs-ON " * beta[r] * ", Tregs-OFF " * beta[r]),
          expression("Tregs-ON " * beta[r]^"*" * ", Tregs-OFF " * beta[r]))
pnC = ggplot(C, aes(tau, med, colour = run, fill = run)) +
  geom_hline(yintercept = 1, colour = COL[["OFF"]], linewidth = 0.8) +
  T_PEAK_LINE +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.18, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_manual(values = COL[c("same", "star")], labels = LAB_C) + scale_fill_manual(values = COL[c("same", "star")], labels = LAB_C) +
  TB + scale_y_log10(breaks = c(0.5, 0.7, 1, 1.5, 2), minor_breaks = seq(0.35, 2.4, by = 0.05)) +
  labs(x = expression("Time / " * t[peak]), y = expression(P[LP] * ", Tregs-ON / Tregs-OFF")) + th +
  theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.4),
        panel.grid.minor = element_line(colour = "grey93", linewidth = 0.25),
        ## the key above the panel, as in A: the row already carries A's two-line key, so
        ## the white space is there and the panel keeps its full height
        legend.position = "top", legend.justification = "left",
        ## box.spacing well above A's 2 pt: C's top tick label (2.0) reaches the panel edge,
        ## so the key needs clearance under it
        legend.margin = margin(0, 0, 0, -BASE_SZ * 2.0), legend.box.spacing = unit(22, "pt"),
        legend.text = element_text(size = BASE_SZ * 0.80), legend.key.width = unit(BASE_SZ * 0.06, "cm"),
        legend.key.spacing.y = unit(0, "pt")) +
  guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))

## ---- D: the demonstration unit, four arms (lib/fig_demo_case.R) ---------------------------
assign("BRANCH", Sys.getenv("BRANCH", "pathogenic_harm"), envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R")); suppressMessages(source("./lib/extinction_floor.R"))
source("./lib/fig_demo_case.R")
DEMO = fig_demo_case(UNIT, BASE_SZ)

## ---- layout: A, B, C in one row; D (three facets under one tag) below -------------------
## Each row is wrapped on its own, so D's wide right margin (its keys) does not push the
## panels of the first row; the tags are therefore set on the panels themselves.
TAG = function(g, l) g + labs(tag = l) + theme(plot.tag = element_text(size = BASE_SZ * 1.2, face = "bold"))
## the story starts with the headroom, so the beta_r* panel is A and the order of events B
ROW1 = (TAG(pnB, "A") | TAG(pnA, "B") | TAG(pnC, "C")) + plot_layout(widths = c(0.85, 1.15, 1.0))
ROW2 = (TAG(DEMO$pC, "D") | DEMO$pA | DEMO$pB) + plot_layout(guides = "collect") &
       theme(legend.position = "right", legend.spacing.y = unit(BASE_SZ * 0.9, "pt"), legend.box.just = "left")
Gp = wrap_elements(full = ROW1, ignore_tag = TRUE) / wrap_elements(full = ROW2, ignore_tag = TRUE) + plot_layout(heights = c(1, 1))
dir.create(dirname(OUT_FIG), showWarnings = FALSE, recursive = TRUE)
png(OUT_FIG, width = 6400, height = 4400, res = 200); print(Gp); invisible(dev.off())
cat(sprintf("[Fig. 4] wrote %s\n", OUT_FIG))
print(B[order(b_off), .(b_off, beta_star_median = signif(med, 3), ratio = round(ratio, 2), n)])
for (tt in c(0.5, 1, 2, 5, 10, 50)) { k = which.min(abs(C[run == "star"]$tau - tt))
  cat(sprintf("[Fig. 4] C at tau %4.1f: same %.3f  star %.3f\n", tt, C[run == "same"]$med[k], C[run == "star"]$med[k])) }
