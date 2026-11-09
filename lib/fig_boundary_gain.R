## ============================================================================
## lib/fig_boundary_gain.R -- Supplementary Fig. 1B/C: displacement of the basin boundary at
## t = 0 against the rescue gap, one branch per panel.
##
## THE METRIC (per combination and coordinate X, from pipeline/23_corridor_position.R)
##      G_X = 100 * (x_on - x_off) / (x_off - base)          (%)
## how much further the state can be displaced along X before the outcome flips with
## Tregs, as a percentage of how far it could be displaced without them. The
## denominator is the Tregs-OFF wall distance from the base state, not the
## compartment capacity: dividing by K_X would multiply the Treg effect by the
## system's baseline fragility (x_off/K_X).
##
## ROWS PLOTTED: treg_effect == "help"; coordinate in (Ea, En, R, D, M1, PLP --
## the last on the pathogenic branch only, where it exists);
## q_kind == "measured" (both walls
## located); G_X finite and > 0; bistable; rescue gap finite, uncensored and above
## the bisection resolution floor STEP = PROBE_BY / 2^ceiling(log2(PROBE_BY/EH0_TOL))
## = 0.078125 (rows with rescue_gap <= 1.01 STEP, including exact zeros, are
## excluded: a log axis cannot show 0, and at the floor the bisection could not
## separate the two thresholds).
##
## THE HELP CLASS. The two corridor pools are not matched in Treg effect: the
## sterile branch is defined by r_P = 0 alone and so carries its help and neutral
## combinations together, while the pathogenic branches are defined by the effect
## and the pathogenic neutral class belongs to no branch at all. Keeping
## treg_effect == "help" makes both panels one class, as every other figure of the
## paper is (figure_scripts/Fig_3.R, figure_scripts/Extended_Data_Fig_1.R). It is a no-op on the pathogenic
## pool, which is help by construction, and on the sterile pool it drops the
## neutral combinations: 3,900 -> 2,805 plotted combinations, with no Spearman rho
## moving by more than 0.004 (help 0.937-0.942 against 0.935-0.939 pooled).
## treg_effect is not in the corridor file; it is joined from selected_big.rds on
## (tau_id, Kg, beta_ros, r_L_P), the same key 23_corridor_position.R uses.
## Per coordinate: Spearman's rho (cor.test, exact = FALSE, asymptotic p) of G_X
## against the rescue gap, and an OLS line of log10(gap) on log10(G_X)
## (geom_smooth(method = "lm") on log-log scales). Drawn points: NPT_C = 9000 rows
## sampled (seed 1) from the rows pooled over coordinates; fits and rho use all rows.
##
## use: source("./lib/fig_boundary_gain.R")
##      gp = gain_panel("sterile", "_round1", y_lim = NULL, base_size = 20)
##      gp$plot (ggplot), gp$stat (rho, p, n per coordinate), gp$counts (filter ledger)
## in : ./data/corridor_q_<branch><corr_tag>.rds
## ============================================================================
suppressMessages({library(data.table); library(ggplot2)})

gain_panel = function(branch, corr_tag = "", y_lim = NULL, base_size = 13,
                      npt = 9000L, probe_by = 10, eh0_tol = 0.1) {
  COORDS = c("Ea","En","R","D","M1","PLP")
  ## five cool hues for the host coordinates; Ea and En take the adjacent pair
  ## (teal, green). PLP, the pathogen and the one coordinate that exists on a single
  ## branch, is warm so that it cannot be mistaken for one of them.
  PAL_COORD = c(M1 = "#2C1B5E", D = "#A05EC7", R = "#2469C9", Ea = "#1DA6A0", En = "#5AC85A",
                PLP = "#D55E00")
  IN = sprintf("./data/corridor_q_%s%s.rds", branch, corr_tag)
  if (!file.exists(IN)) stop("no ", IN, " -- run pipeline/23_corridor_position.R first")
  L = as.data.table(readRDS(IN))
  ## ---- one Treg-effect class, joined on the corridor key ---------------------
  CK = c("tau_id", "Kg", "beta_ros", "r_L_P")
  SELF = "./data/selected_big.rds"
  if (!file.exists(SELF)) stop("no ", SELF, " -- treg_effect cannot be joined")
  EFF = unique(as.data.table(readRDS(SELF))[, c(CK, "treg_effect"), with = FALSE])
  if (uniqueN(EFF, by = CK) < nrow(EFF)) stop("treg_effect is not unique on the corridor key")
  n0 = nrow(L)
  L = merge(L, EFF, by = CK, all.x = TRUE)
  if (nrow(L) != n0) stop("the treg_effect join changed the row count: ", n0, " -> ", nrow(L))
  if (anyNA(L$treg_effect)) stop(sum(is.na(L$treg_effect)), " corridor rows have no treg_effect")
  cat(sprintf("[gain_panel] %s | classes measured: %s\n", branch,
              paste(sprintf("%s %s", names(table(L$treg_effect)),
                            format(as.integer(table(L$treg_effect)), big.mark = ",")), collapse = " | ")))
  L = L[treg_effect == "help"]
  if (!nrow(L)) stop("no help-class rows in ", IN)
  ## joined to the kinetic parameters of each vector; the join also fixes the row order
  ## (sorted by tau_id) from which the drawn points are sampled
  tf = sort(Sys.glob("./data/tau_grid_v6_*.rds"))
  TG = rbindlist(lapply(tf, function(p) as.data.table(readRDS(p))), fill = TRUE)
  D = merge(L[coord %in% COORDS], unique(TG[, .(tau_id, mu_R, sigma_D, tau_D)]), by = "tau_id", all.x = TRUE)
  D[, gain := 100 * (x_on - x_off) / (x_off - base)]
  D[, class := fifelse(rescue_mode == "bistable", "bistable", "mono_B")]
  P = D[q_kind == "measured" & rescue_mode %in% c("bistable","mono_B") & is.finite(gain) & gain > 0]
  G = P[class == "bistable" & is.finite(rescue_gap) & !gap_censored & rescue_gap > 0]
  if (!nrow(G)) stop("no bistable rows with a usable rescue gap")
  STEP = probe_by / 2^ceiling(log2(probe_by / eh0_tol))
  G[, at_floor := rescue_gap <= STEP * 1.01]
  GF = G[at_floor == FALSE]
  ## filter ledger per coordinate (bistable rows)
  CNT = D[class == "bistable", .(rows = .N, measured = sum(q_kind == "measured"),
                                 gain_pos = sum(q_kind == "measured" & is.finite(gain) & gain > 0)), by = coord]
  CNT = merge(CNT, G[, .(gap_ok = .N, at_floor = sum(at_floor), plotted = sum(!at_floor)), by = coord], by = "coord", all.x = TRUE)
  ## the coordinate order of the label block (decreasing rho) is what the figure shows
  STAT = GF[, {
    ct = suppressWarnings(cor.test(gain, rescue_gap, method = "spearman", exact = FALSE))
    .(rho = unname(ct$estimate), p = ct$p.value, n = .N) }, by = coord]
  STAT = STAT[order(-rho)]
  if (any(STAT$p >= 2.2e-16)) stop("[gain_panel] a p value is not < 2.2e-16; keep p in the labels")
  SYM06 = c(En = "E[n]", Ea = "E[a]", M1 = "M[1]", R = "R", D = "D", PLP = "P[LP]")
  STAT[, sym := unname(SYM06[as.character(coord)])]
  STAT[, lab := sprintf("rho == '%+.3f'", rho)]
  set.seed(1)
  GS = GF[sample(.N, min(npt, .N))]
  xr = range(log10(GF$gain), na.rm = TRUE); yr = range(log10(GF$rescue_gap), na.rm = TRUE)
  if (!is.null(y_lim)) yr = log10(y_lim)
  STAT[, lx := 10^(xr[1] + 0.03 * diff(xr))]
  STAT[, ly := 10^(yr[2] - ((if (is.null(y_lim)) 0.02 else 0.05) + 0.085 * (seq_len(.N) - 1)) * diff(yr))]
  TH = theme_classic(base_size = base_size) +
    theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = base_size),
          axis.text = element_text(size = base_size * 0.85, colour = "grey20"),
          plot.tag = element_text(size = base_size * 1.4, face = "plain"),
          legend.position = "none", plot.title = element_blank())
  pC = ggplot() +
    geom_point(data = GS, aes(gain, rescue_gap), colour = "#B8BCC6", size = 0.6, alpha = 1, shape = 16) +
    geom_smooth(data = GF, aes(gain, rescue_gap, colour = coord), method = "lm",
                se = FALSE, linewidth = 0.8, formula = y ~ x) +
    geom_text(data = STAT, aes(x = lx, y = ly, label = sym, colour = coord), parse = TRUE,
              hjust = 0, size = base_size * 0.27, fontface = "bold", inherit.aes = FALSE, show.legend = FALSE) +
    geom_text(data = STAT, aes(x = 10^(log10(lx) + 0.10 * diff(xr)), y = ly, label = lab, colour = coord), parse = TRUE,
              hjust = 0, size = base_size * 0.27, fontface = "bold", inherit.aes = FALSE, show.legend = FALSE) +
    scale_colour_manual(values = PAL_COORD) +
    ## plain decimal tick labels: the default switches to 1e-01 notation once the
    ## range is wide enough, which it is on the pathogenic panel but not the sterile
    ## one, and the two panels would then be labelled differently
    scale_x_log10(labels = function(z) format(z, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)) +
    (if (is.null(y_lim)) scale_y_log10() else scale_y_log10(limits = y_lim, breaks = 10^(-2:2), expand = c(0, 0))) +
    labs(x = "Additional displacement tolerated\nat t = 0 (%)",
         y = "Rescue gap") + TH
  list(plot = pC, stat = STAT, counts = CNT, rows = GF, step = STEP, file = IN)
}
