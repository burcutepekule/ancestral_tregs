## ============================================================================
## Supplementary_Fig_1.R -- Supplementary Fig. 1: displacement of the basin boundary.
##   A  (left column, full height) one sterile bistable combination followed in
##      time: the basin boundaries of Tregs OFF and Tregs ON along En re-located at
##      sampled times along the trajectory, and three initial conditions (below,
##      inside, above the corridor), with brackets under the trajectory ends
##      grouping them by the attractor each one reaches
##   B  boundary displacement at t = 0 against the rescue gap, sterile
##   C  the same, pathogenic help
## A: lib/corridor_defs.R. B and C: lib/fig_boundary_gain.R.
##
## The attractor of each trajectory is the classifier's label (eng$lab), not its
## position at t = T_END.
##
## WHICH CORRIDOR FILE EACH PANEL READS. The sterile corridor exists only as the
## round-1 subsample (corridor_q_sterile_round1.rds: 16,990 of the 89,328 pool
## combinations); the pathogenic-help corridor is complete (corridor_q_pathogenic_help.rds).
## Per branch: the final file when present, the round-1 file otherwise;
## CORR_TAG_<BRANCH> forces one. The file read is printed.
##
## PANEL A's BOUNDARY SWEEP (60 times x 6 trajectories, each a 40-point scan plus 20
## bisections of 12,000-unit integrations) is cached in
##   ./data/cache_bnd_v6_v6_10621_15_0.0015_0_En_60.46609_76.79843_107.11795_60_1000_40_20_12000_sterile_cheb1.rds
## and read from there. Without the cache it is recomputed (hours; set N_CORES).
## Choosing the example (rule "gap0") and labelling the attractors: about a minute.
##
## in : ./data/corridor_q_sterile_round1.rds, ./data/corridor_q_pathogenic_help.rds,
##      ./data/tau_grid_v6_<n>.rds, the boundary cache above
## out: ./figures/supplementary_information/Supplementary_Fig_1.png, ./source_data/Supplementary_Fig_1_example.csv
## RUN: N_CORES=6 Rscript figure_scripts/Supplementary_Fig_1.R
## ============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
Sys.setenv(COMPOSITE = "1")
CT_BOTH = Sys.getenv("CORR_TAG")
ctag = function(br) {
  v  = Sys.getenv(paste0("CORR_TAG_", toupper(br)))
  tg = if (nzchar(v)) v else if (nzchar(CT_BOTH)) CT_BOTH else
       if (file.exists(sprintf("./data/corridor_q_%s.rds", br))) "" else "_round1"
  f = sprintf("./data/corridor_q_%s%s.rds", br, tg)
  if (!file.exists(f)) stop(sprintf("[Supp. Fig. 1] %s: %s not found -- run pipeline/23_corridor_position.R for it", br, f))
  cat(sprintf("[Supp. Fig. 1] %-16s reads %s\n", br, f)); tg }
CTAG = vapply(c(sterile = "sterile", pathogenic_help = "pathogenic_help"), ctag, "")
if (!nzchar(Sys.getenv("RAY_CLUST_TOL"))) Sys.setenv(RAY_CLUST_TOL = "1")
if (!nzchar(Sys.getenv("T_END")))   Sys.setenv(T_END = "1000")
if (!nzchar(Sys.getenv("N_CORES"))) Sys.setenv(N_CORES = "6")
if (!nzchar(Sys.getenv("BASE_SIZE"))) Sys.setenv(BASE_SIZE = "30")
BS = as.numeric(Sys.getenv("BASE_SIZE"))
WIN = Sys.getenv("WIN", "En"); EX_PICK = Sys.getenv("EX_PICK", "gap0")
tagged = function(p, tg, ttl = NULL) p + labs(tag = tg, subtitle = ttl) + theme(plot.tag = element_text(size = 1.15 * BS, face = "bold"), plot.tag.position = c(0.01, 0.99), plot.subtitle = element_text(size = 0.9 * BS, face = "bold", hjust = 0))

## ---- B and C ----------------------------------------------------------------
source("./lib/fig_boundary_gain.R")
built_y = function(g) ggplot_build(g)$layout$panel_params[[1]]$y.range   # log10 units
YL = 10^range(c(built_y(gain_panel("sterile", CTAG[["sterile"]], base_size = BS)$plot),
                built_y(gain_panel("pathogenic_help", CTAG[["pathogenic_help"]], base_size = BS)$plot)))
GB = gain_panel("sterile", CTAG[["sterile"]], y_lim = YL, base_size = BS)
GC = gain_panel("pathogenic_help", CTAG[["pathogenic_help"]], y_lim = YL, base_size = BS)
pB = tagged(GB$plot, "B", "Sterile")
pC = tagged(GC$plot, "C", "Pathogenic")

## ---- A ----------------------------------------------------------------------
Sys.setenv(BRANCH = "sterile", CORR_TAG = CTAG[["sterile"]]); BRANCH = "sterile"; source("./lib/corridor_defs.R")
COL_SEED  = c(below = "#261CC1", corridor = "#c2189a", above = "#FF7070")
COL_BASIN = c(resolved = "#261CC1", chronic = "#FF7070")
Lc = load_corridor(BRANCH); pick = pick_combo(Lc, WIN, EX_PICK); eng = combo_engine(pick); sds = make_seeds(pick)
TR = run_traj(eng, WIN, sds); BND = sweep_bnd(eng, WIN, pick, sds, TR)
pA = panel_corridor(WIN, TR, BND)
## the numbers behind A, for the methods
dir.create("./source_data", showWarnings = FALSE)
fwrite(data.table(tau_id = pick$tau_id, batch = pick$tau_id %/% 10000L, row_ind = pick$tau_id %% 10000L, Kg = pick$Kg, beta_ros = pick$beta_ros,
                  rescue_gap = pick$rescue_gap, e0_base = pick$e0_base, base = pick$base, x_off = pick$x_off, x_on = pick$x_on,
                  seed_below = sds[["below"]], seed_corridor = sds[["corridor"]], seed_above = sds[["above"]],
                  bnd_t0_off = BND[t == min(t) & theta == 0 & is.finite(b), b][1], bnd_t0_on = BND[t == min(t) & theta == 1 & is.finite(b), b][1]),
       "./source_data/Supplementary_Fig_1_example.csv")

## each trajectory's attractor, and where it sits along WIN at t = T_END
ENDS = TR[time <= T_END, .SD[.N], by = .(seed_key, theta, system)]
ENDS[, basin := vapply(seq_len(.N), function(i)
  eng$lab(eng$SETC(eng$ybase, WIN, sds[[as.character(seed_key[i])]]), theta[i]), integer(1))]
ENDS[, x := get(WIN)]
print(ENDS[, .(seed_key, system, time, x, basin)])
if (anyNA(ENDS$basin)) stop("[Supp. Fig. 1] a trajectory has no attractor label")

## the brackets sit in a strip added below t = T_END (the y axis is reversed);
## each is centred on its group's end points and spans them, widened to MIN_W of the
## x range when narrower; the low-burden ends sit on the axis, so its bracket and label
## are centred there and may run past the panel's left edge (clip is off)
XR = ggplot_build(pA)$layout$panel_params[[1]]$x.range
MIN_W = 0.05 * diff(XR)
BRK = ENDS[, .(lo = min(x), hi = max(x)), by = basin]
BRK[, `:=`(mid = (lo + hi) / 2, w = pmax(hi - lo, MIN_W))]
BRK[, `:=`(lo = mid - w / 2, hi = pmin(mid + w / 2, XR[2] - 0.01 * diff(XR)))]
BRK[, `:=`(mid = (lo + hi) / 2, txt = fifelse(basin == 1L, "Low-burden\nattractor", "High-burden\nattractor"))]
## strip heights in time units on the pseudo-log axis: ticks drop from y_top to y_bar
fwd = function(t) asinh(t / 2) / log(10); inv = function(u) 2 * sinh(u * log(10))
U = fwd(T_END); y_top = inv(U + 0.04); y_bar = inv(U + 0.12); y_txt = inv(U + 0.18); y_lim = inv(U + 0.55)
pA = suppressMessages(pA + coord_cartesian(xlim = XR, ylim = c(y_lim, 0), expand = FALSE, clip = "off")) +
  annotate("rect", xmin = XR[1], xmax = XR[2], ymin = T_END, ymax = y_lim, fill = "white") +
  geom_segment(data = BRK, aes(x = lo, xend = hi, y = y_bar, yend = y_bar), linewidth = 0.9) +
  geom_segment(data = BRK, aes(x = lo, xend = lo, y = y_top, yend = y_bar), linewidth = 0.9) +
  geom_segment(data = BRK, aes(x = hi, xend = hi, y = y_top, yend = y_bar), linewidth = 0.9) +
  geom_text(data = BRK, aes(x = mid, y = y_txt, label = txt), vjust = 1, lineheight = 0.9, size = 0.85 * BS / .pt) +
  theme(panel.grid.minor.y = element_blank())
pA = tagged(pA, "A", "Sterile, one combination")

p = pA | (pB / pC)
OUT = Sys.getenv("OUT", "./figures/supplementary_information/Supplementary_Fig_1.png"); dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
png(OUT, width = 3600, height = 3024, res = 200); print(p); invisible(dev.off()); cat("[Supp. Fig. 1] wrote", OUT, "\n")
