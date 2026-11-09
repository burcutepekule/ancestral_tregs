## ============================================================================
## lib/fig_imperfect_treg.R -- drawing code shared by Fig. 5 and
## Extended Data Figs. 6 and 7 (the imperfect Treg).
## ----------------------------------------------------------------------------
## imperfect_treg_panels(branch) builds six panels from the rescue-probability
## tables written by pipeline/63_imperfect_rescue_prob.R:
##   pA  heatmap of P(rescue), the mean over crossing units of p_open, in each
##       cell of stability (1 - phi) against eta
##   pB  one-step forward difference of that surface along eta,
##       P(1 - phi, eta + 0.1) - P(1 - phi, eta), at the step midpoint
##   pC  the same along stability, P(1 - phi + 0.1, eta) - P(1 - phi, eta)
##       (B and C share one grey scale, 0 to the largest difference in either)
##   pD  per unit, p_open averaged over eta, against 1 - phi: mean over units
##       (line) and N_BANDS nested central quantile bands (5-95% ... 45-55%)
##   pE  the same against eta, averaged over 1 - phi
##   pF  the two mean curves of D and E together
## all with the identity line on D, E and F. The CROSSING units are those whose
## p_open is neither 0 in every cell nor 1 in every cell (0 < p_mean < 1). Two
## three further filters, all on by default, reduce the crossing set to the units
## the panels are drawn from: DROP_NEUTRAL=0 keeps the units whose Treg effect is
## not "help", MIN_RANGE=0 keeps the units whose marginals both span less than
## 0.5 in P(rescue), and DROP_NONMONO=0 keeps the units whose per-unit marginal
## P(rescue) falls somewhere along either axis (see the notes at the filters).
## branch = "" pools both branches; "sterile" or "pathogenic_help" restricts
## to that branch's units (from ./data/fig05_units.csv).
##
## Fig. 5 = the aligned A and F; Extended Data Fig. 7 = B, C and single combinations
## of this analysis. Panel letters are reassigned per figure. Passing a branch
## restricts the population to it, which the Methods tables report as numbers.
## ============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

IMP_BASE_SZ = as.numeric(Sys.getenv("BASE_SIZE", "13"))
IMP_FIG_W   = as.numeric(Sys.getenv("FIG_W", "5000"))
IMP_FIG_H   = as.numeric(Sys.getenv("FIG_H", "3300"))
IMP_FIG_RES = as.numeric(Sys.getenv("FIG_RES", "200"))

## panel titles are "<bold letter>   <body>"; the letter is set per figure
imp_ttl = function(L, body) bquote(bold(.(L)) * "   " * .(body))
IMP_T = list(
  A = quote("P(rescue)"),
  B = quote("P(rescue), one step along " * eta),
  C = quote("P(rescue), one step along " * (1 - phi)),
  D = quote("P(rescue) along " * (1 - phi)),
  E = quote("P(rescue) along " * eta),
  F = quote("P(rescue)"))

imperfect_treg_panels = function(branch = "") {
  BASE_SZ = IMP_BASE_SZ
  TXT     = as.numeric(Sys.getenv("TILE_TXT", "0.42"))   # tile-number size, x BASE_SZ
  ID_LW   = as.numeric(Sys.getenv("ID_LW", "1.3"))        # identity line weight
  ID_COL  = Sys.getenv("ID_COL", "#B80000")               # identity line colour
  ID_FAN  = Sys.getenv("ID_FAN", ID_COL)
  N_BANDS = as.integer(Sys.getenv("N_BANDS", "9"))        # nested central intervals
  FAN_A   = as.numeric(Sys.getenv("FAN_ALPHA", "0.06"))   # alpha of each band
  FAN_LW  = as.numeric(Sys.getenv("FAN_LW", "0"))         # band outline weight, 0 = none

  CELL = as.data.table(readRDS("./data/fig05_switch_prob_cells.rds"))
  U    = as.data.table(readRDS("./data/fig05_switch_prob_units.rds"))
  if (nzchar(branch)) {
    UN = fread("./data/fig05_units.csv")
    SU = UN$unit[UN$branch == branch]
    U = U[unit %in% SU]; CELL = CELL[unit %in% SU]
    cat(sprintf("[imperfect] branch = %s: %d units\n", branch, nrow(U)))
  }
  cross = U[p_mean > 0 & p_mean < 1, unit]
  cat(sprintf("[imperfect] %s: %d units, %d never rescued, %d always rescued, %d crossing\n",
              if (nzchar(branch)) branch else "pooled", nrow(U), U[p_mean == 0, .N],
              U[p_mean == 1, .N], length(cross)))
  ## ---- the two population filters ------------------------------------------
  ## DROP_NEUTRAL: units whose Tregs-ON low-burden attractor sits within the
  ## help/harm tolerance of the Tregs-OFF one (treg_effect == "neutral" in
  ## data/rescue_<branch>.rds). All of them are on the sterile branch.
  ## DROP_NONMONO: units whose per-unit marginal P(rescue) -- the mean over the
  ## other axis, the quantity panels D, E and F draw -- falls at any step of
  ## either axis. The perfect Treg is then not the best Treg of that unit's own
  ## grid, so the unit cannot speak to the shape of the decline away from it.
  ## NOTE ON THE SECOND FILTER. The cell probabilities rest on at most ten
  ## replicates, so a marginal can fall by one replicate through the noise
  ## alone; under a monotone null, re-drawing ten replicates flags 162 of the
  ## 717 crossing units on average against the 148 observed. The false-positive
  ## rate also rises with the unit's overall rescue rate (0.04 below p = 0.2,
  ## 0.62 above p = 0.8), so the filter removes the readily rescued units
  ## preferentially and the whole surface sits lower than the unfiltered one.
  if (!identical(Sys.getenv("DROP_NEUTRAL", "1"), "0")) {
    ## the keys are doubles, and ./data/fig05_units.csv does not round-trip
    ## beta_ros exactly (0.0045 comes back one ulp low, and 220 of the 977 units
    ## then fail to match), so the per-branch units files are read instead
    BRS = c("sterile", "pathogenic_help")
    CKU = c("tau_id", "Kg", "beta_ros", "r_L_P")
    UB = rbindlist(lapply(BRS, function(b) {
           x = as.data.table(readRDS(sprintf("./data/fig05_units_%s.rds", b)))
           x[, c("unit", CKU), with = FALSE][, branch := b] }), fill = TRUE)
    RE = rbindlist(lapply(BRS, function(b) {
           x = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", b))); x[, branch := b]
           x[, c(CKU, "branch", "treg_effect"), with = FALSE] }), fill = TRUE)
    EFF = merge(UB, RE, by = c(CKU, "branch"))
    if (nrow(EFF) != nrow(UB))
      stop("the Treg-effect merge lost ", nrow(UB) - nrow(EFF), " of ", nrow(UB), " units")
    drop = intersect(cross, EFF[treg_effect != "help", unit])
    cross = setdiff(cross, drop)
    cat(sprintf("[imperfect] DROP_NEUTRAL: %d units of effect other than help removed, %d left\n",
                length(drop), length(cross)))
  }
  ## MIN_RANGE: a unit whose marginals barely move carries no shape to read. The
  ## test is on the per-unit marginals, the curves the panels draw, and not on
  ## the 121 cells: a unit rescued in every cell but one has cells spanning the
  ## whole of [0,1] while its marginals move by 1/11, so a cell-range test
  ## cannot see it. A unit is kept when the larger of its two marginal ranges
  ## reaches MIN_RANGE, so a unit whose stability curve spans the axis is kept
  ## however flat its eta curve is; a per-axis minimum would instead bind almost
  ## only on eta, whose per-unit range has median 0.15 against 0.98 along
  ## stability, and would keep only the units in which eta matters.
  MIN_RANGE = as.numeric(Sys.getenv("MIN_RANGE", "0.5"))
  if (MIN_RANGE > 0) {
    CC = CELL[unit %in% cross]
    rng = function(ax) CC[, .(p = mean(p_open)), by = c("unit", ax)][
      , .(r = max(p) - min(p)), by = unit]
    G = merge(rng("stab"), rng("eta"), by = "unit")[, .(r = max(r.x, r.y)), by = unit]
    drop = G[r < MIN_RANGE, unit]
    cross = setdiff(cross, drop)
    cat(sprintf("[imperfect] MIN_RANGE = %.2f: %d units whose marginals span less removed, %d left\n",
                MIN_RANGE, length(drop), length(cross)))
  }
  if (!identical(Sys.getenv("DROP_NONMONO", "1"), "0")) {
    CC = CELL[unit %in% cross]
    fell = function(ax) CC[, .(p = mean(p_open)), by = c("unit", ax)][
      order(unit, get(ax))][, .(nm = any(diff(p) < -1e-12)), by = unit][nm == TRUE, unit]
    drop = union(fell("stab"), fell("eta"))
    cross = setdiff(cross, drop)
    cat(sprintf("[imperfect] DROP_NONMONO: %d units with a falling marginal removed, %d left\n",
                length(drop), length(cross)))
  }
  H = CELL[unit %in% cross, .(p = mean(p_open)), by = .(stab, eta)]
  H[, `:=`(stab = round(stab, 1), eta = round(eta, 1))]
  setorder(H, stab, eta)

  ## finite differences, placed at the midpoint of the step
  DE = H[, .(eta = head(eta, -1) + 0.05, d = diff(p)), by = stab]    # along eta, per row
  DS = H[, .(stab = head(stab, -1) + 0.05, d = diff(p)), by = eta]   # along stability, per column
  lim = c(0, max(c(DE$d, DS$d)))
  cat(sprintf("[imperfect] one-step change in P(rescue): along eta median %.3f max %.3f | along stability median %.3f max %.3f\n",
              median(DE$d), max(DE$d), median(DS$d), max(DS$d)))

  AXT = theme(axis.title   = element_text(size = BASE_SZ * 1.9),
              axis.text    = element_text(size = BASE_SZ * 1.45),
              legend.title = element_text(size = BASE_SZ * 1.5),
              legend.text  = element_text(size = BASE_SZ * 1.3))
  thm = theme_minimal(base_size = BASE_SZ) + AXT +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = BASE_SZ * 1.7),
          plot.subtitle = element_text(size = BASE_SZ * 1.25, colour = "grey25"))
  pA = ggplot(H, aes(eta, stab, fill = p)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", p), colour = p < 0.55), size = BASE_SZ * TXT) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
    scale_fill_viridis_c(limits = c(0, 1), name = "P(rescue)") +
    guides(fill = guide_colourbar(barheight = unit(12, "lines"), barwidth = unit(1.6, "lines"))) +
    scale_x_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
    labs(title = imp_ttl("A", IMP_T$A), x = expression(eta), y = expression(1 - phi)) +
    thm + theme(legend.position = "right")
  gp = function(D, xv, yv, ttl, sub, xlab, ylab, leg = "right") {
    ggplot(D, aes(get(xv), get(yv), fill = d)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%.2f", d), colour = d > 0.6 * lim[2]), size = BASE_SZ * TXT) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
      scale_fill_gradient(low = "white", high = "black", limits = lim, name = "change in\nP(rescue)\nper step") +
      scale_x_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
      guides(fill = guide_colourbar(barheight = unit(12, "lines"), barwidth = unit(1.6, "lines"))) +
      labs(title = ttl, subtitle = sub, x = xlab, y = ylab) + thm + theme(legend.position = leg)
  }
  pB = gp(DE, "eta", "stab", imp_ttl("B", IMP_T$B),
          expression(P(1 - phi, eta + 0.1) - P(1 - phi, eta)),
          expression(eta * "  (step midpoint)"), expression(1 - phi), leg = "none")
  pC = gp(DS, "eta", "stab", imp_ttl("C", IMP_T$C),
          expression(P(1 - phi + 0.1, eta) - P(1 - phi, eta)),
          expression(eta), expression((1 - phi) * "  (step midpoint)"), leg = "right")

  ## per-unit marginals of p_open: mean over the other axis
  CC = CELL[unit %in% cross]
  MS = CC[, .(p = mean(p_open)), by = .(unit, x = stab)]
  ME = CC[, .(p = mean(p_open)), by = .(unit, x = eta)]
  LO = seq(0.05, 0.45, length.out = N_BANDS)               # lower quantile of each band
  qb = function(M) M[, .(mean = mean(p), med = median(p)), by = x][order(x)]
  fan = function(M) rbindlist(lapply(LO, function(q)
    M[, .(lo = quantile(p, q), hi = quantile(p, 1 - q), band = q), by = x]))[order(band, x)]
  QS = qb(MS); QE = qb(ME); FS = fan(MS); FE = fan(ME)
  band = function(Q, Fn, col, ttl, xlab, idcol = ID_FAN) {
    ggplot(Q, aes(x)) +
      geom_ribbon(data = Fn, aes(x, ymin = lo, ymax = hi, group = band), fill = col, alpha = FAN_A,
                  colour = if (FAN_LW > 0) col else NA, linetype = "dashed",
                  linewidth = FAN_LW, outline.type = "both") +
      geom_abline(slope = 1, intercept = 0, colour = idcol, linetype = "dotdash", linewidth = ID_LW) +
      geom_line(aes(y = mean), colour = "black", linewidth = 1.6) +
      geom_point(aes(y = mean), colour = "black", size = 2) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      labs(title = ttl, x = xlab, y = "P(rescue), mean over the other axis") +
      theme_minimal(base_size = BASE_SZ) + AXT +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(size = BASE_SZ * 1.7),
            plot.subtitle = element_text(size = BASE_SZ * 1.25, colour = "grey25"))
  }
  pD = band(QS, FS, "#4481B5", imp_ttl("D", IMP_T$D), expression(1 - phi))
  pE = band(QE, FE, "#8A21A2", imp_ttl("E", IMP_T$E), expression(eta))
  BOTH = rbind(QS[, .(x, mean, axis = "1 - phi")], QE[, .(x, mean, axis = "eta")])
  pF = ggplot(BOTH, aes(x, mean, colour = axis)) +
    geom_abline(slope = 1, intercept = 0, colour = ID_COL, linetype = "dotdash", linewidth = ID_LW) +
    geom_line(linewidth = 1.6) + geom_point(size = 2) +
    scale_colour_manual(values = c(`1 - phi` = "#4481B5", eta = "#8A21A2"),
                        labels = c(expression(1 - phi), expression(eta)), name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = imp_ttl("F", IMP_T$F),
         x = "axis value", y = "P(rescue), mean over the other axis") +
    theme_minimal(base_size = BASE_SZ) + AXT +
    theme(panel.grid.minor = element_blank(), legend.position = c(0.18, 0.82),
          legend.text = element_text(size = BASE_SZ * 1.5),
          legend.title = element_text(size = BASE_SZ * 1.5),
          legend.key.size = unit(2.2, "lines"),
          legend.text.align = 0,
          plot.title = element_text(size = BASE_SZ * 1.7),
          plot.subtitle = element_text(size = BASE_SZ * 1.25, colour = "grey25"))
  list(pA = pA, pB = pB, pC = pC, pD = pD, pE = pE, pF = pF,
       H = H, DE = DE, DS = DS, QS = QS, QE = QE, n_units = nrow(U), n_cross = length(cross),
       units = cross)
}

## the supplementary layout: B, C over E, D, relettered A-D, so that each column
## is one axis: left eta (step along eta over the marginal along eta), right
## stability (step along 1 - phi over the marginal along 1 - phi)
imperfect_treg_supp = function(P)
  (P$pB + labs(title = imp_ttl("A", IMP_T$B)) | P$pC + labs(title = imp_ttl("B", IMP_T$C))) /
  (P$pE + labs(title = imp_ttl("C", IMP_T$E)) | P$pD + labs(title = imp_ttl("D", IMP_T$D)))

## the main layout: A and F, F relettered B
imperfect_treg_main = function(P)
  P$pA + labs(title = imp_ttl("A", IMP_T$A)) | P$pF + labs(title = imp_ttl("B", IMP_T$F))

imperfect_treg_write = function(fig, f, w, h) {
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  png(f, width = w, height = h, res = IMP_FIG_RES); print(fig + plot_annotation()); invisible(dev.off())
  cat(sprintf("[imperfect] wrote %s\n", f))
}

## ============================================================================
## imperfect_treg_carpets(branch) -- every retained unit, one raster row each,
## in place of a summary across units: the per-unit marginal P(rescue) (the mean
## over the other axis, the quantity the marginal panels average) as a colour,
## with the units ordered by the area under their own curve. A unit that crosses
## once at its own threshold shows as a dark-to-bright row with a sharp edge; a
## unit that gains gradually shows as a graded row. The population is the one
## imperfect_treg_panels() draws, under the same filters.
## ============================================================================
imperfect_treg_carpet_data = function(branch = "") {
  P = imperfect_treg_panels(branch)
  CELL = as.data.table(readRDS("./data/fig05_switch_prob_cells.rds"))
  CELL = CELL[unit %in% P$units, .(unit, stab = round(stab, 1), eta = round(eta, 1), p = p_open)]
  mk = function(ax) {
    M = CELL[, .(p = mean(p)), by = c("unit", ax)]
    setnames(M, ax, "x")
    A = M[, .(a = sum(diff(sort(unique(x))) * (head(p[order(x)], -1) + tail(p[order(x)], -1)) / 2)),
          by = unit][order(a)][, r := .I]
    merge(M, A[, .(unit, r)], by = "unit")
  }
  list(S = mk("stab"), E = mk("eta"), n = length(P$units))
}

imperfect_treg_carpets = function(branch = "", letters = c("A", "B")) {
  D = imperfect_treg_carpet_data(branch)
  one = function(M, L, ttl, xlab) ggplot(M, aes(x, r, fill = p)) +
    geom_raster() +
    scale_fill_viridis_c(limits = c(0, 1), name = "P(rescue)") +
    guides(fill = guide_colourbar(barheight = unit(12, "lines"), barwidth = unit(1.6, "lines"))) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = imp_ttl(L, ttl), x = xlab,
         y = sprintf("unit (%d), ordered by area under its own curve", D$n)) +
    theme_minimal(base_size = IMP_BASE_SZ) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = IMP_BASE_SZ * 1.7),
          axis.title = element_text(size = IMP_BASE_SZ * 1.55),
          axis.text = element_text(size = IMP_BASE_SZ * 1.3),
          legend.title = element_text(size = IMP_BASE_SZ * 1.5),
          legend.text = element_text(size = IMP_BASE_SZ * 1.3))
  list(pS = one(D$S, letters[1], quote("P(rescue) per unit along " * (1 - phi)), expression(1 - phi)),
       pE = one(D$E, letters[2], quote("P(rescue) per unit along " * eta), expression(eta)))
}

## ============================================================================
## THE SHAPE OF THE RESPONSE, over combinations rather than averaged into one
## curve. imperfect_treg_shape_data() builds the per-combination marginals and
## their aligned versions; imperfect_treg_fourpanel() draws an earlier four-panel version from them
## and imperfect_treg_composition() the band-composition panels (not in the paper).
## ----------------------------------------------------------------------------
## WHY ALIGNING. The mean marginal (Extended Data Fig. 6) cannot distinguish a population
## of smooth curves from a population of steps whose thresholds differ, because
## averaging over the thresholds blurs the step. Shifting each curve onto its
## own midpoint first removes exactly that blur.
## WHY TRIMMING. A curve covers a different stretch of the shifted axis
## depending on where its midpoint sits, so at the extremes the mean and the
## quantiles would be taken over a shrinking and self-selected subset -- the
## combinations whose midpoint is early, which are not representative of the
## rest. The aligned panels are therefore drawn only where at least TRIM of the
## curves are defined.
## ============================================================================
IMP_AX_COL = c(stability = "#4481B5", specificity = "#8A21A2")

imperfect_treg_shape_data = function(branch = "", trim = 0.6, step = 0.05) {
  P = imperfect_treg_panels(branch)
  CELL = as.data.table(readRDS("./data/fig05_switch_prob_cells.rds"))
  X = CELL[unit %in% P$units, .(unit, stab = round(stab, 1), eta = round(eta, 1), p = p_open)]
  M = rbind(X[, .(p = mean(p)), by = .(unit, x = stab)][, axis := "stability"],
            X[, .(p = mean(p)), by = .(unit, x = eta )][, axis := "specificity"])
  M[, axis := factor(axis, levels = c("stability", "specificity"))]
  RNG = M[, .(rng = max(p) - min(p)), by = .(unit, axis)]
  M = merge(M, RNG, by = c("unit", "axis"))
  ## the midpoint of each curve's own rise, by linear interpolation
  cross = function(xx, pp, lev) {
    if (all(pp < lev) || all(pp >= lev)) return(NA_real_)
    i = which(pp >= lev)[1]
    if (i == 1L) return(xx[1])
    xx[i-1] + (lev - pp[i-1]) * (xx[i] - xx[i-1]) / (pp[i] - pp[i-1])
  }
  W = M[order(unit, axis, x), {
    lo = min(p); r = max(p) - lo
    .(x50 = cross(x, p, lo + 0.5 * r)) }, by = .(unit, axis)][is.finite(x50)]
  A = merge(M, W, by = c("unit", "axis"))[, dx := x - x50]
  GRID = seq(-1, 1, by = step)
  IN = A[, { f = approx(dx, p, xout = GRID, rule = 1); .(dbin = GRID, p = f$y) },
         by = .(unit, axis)][is.finite(p)]
  IN = merge(IN, RNG, by = c("unit", "axis"))
  NAX = IN[, .(n_all = uniqueN(unit)), by = axis]
  KEEP = merge(IN[, .(n = .N), by = .(axis, dbin)], NAX, by = "axis")[n >= trim * n_all]
  list(units = P$units, M = M, A = A, IN = IN, KEEP = KEEP, GRID = GRID, trim = trim,
       lims = KEEP[, .(lo = min(dbin), hi = max(dbin), n_min = min(n), n_all = n_all[1]), by = axis])
}

## nine nested central quantile bands, as Extended Data Fig. 6 draws them
imp_fan = function(D, xcol, lo = seq(0.05, 0.45, length.out = 9))
  rbindlist(lapply(lo, function(q)
    D[, .(lo = quantile(p, q), hi = quantile(p, 1 - q), band = q), by = c("axis", xcol)]))

## imp_shape_panels() builds the four marginal panels once and returns them
## named, so that the main figure can take the aligned pair, the supplement the
## unaligned pair, and the archived version0 script all four, without any of the
## three re-implementing the panel. `letters` is recycled in the order
## raw_stability, raw_specificity, aligned_stability, aligned_specificity.
imp_shape_panels = function(S, letters = c("A", "B", "C", "D")) {
  BS = IMP_BASE_SZ
  thm = theme_minimal(base_size = BS) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = BS * 1.7),
          axis.title = element_text(size = BS * 1.75),
          axis.text  = element_text(size = BS * 1.35))
  FAN_RAW = imp_fan(S$M, "x")
  ALN_K = merge(S$IN, S$KEEP[, .(axis, dbin)], by = c("axis", "dbin"))
  FAN_ALN = imp_fan(ALN_K, "dbin"); setnames(FAN_ALN, "dbin", "dx")
  MU_RAW = S$M[, .(p = mean(p)), by = .(axis, x)]
  MU_ALN = ALN_K[, .(p = mean(p)), by = .(axis, dbin)]; setnames(MU_ALN, "dbin", "dx")
  pan = function(Fn, mu, xcol, L, ttl, xlab, xlim, ax, ref) {
    ggplot(Fn, aes(get(xcol))) +
      geom_ribbon(aes(ymin = lo, ymax = hi, group = band),
                  fill = IMP_AX_COL[[ax]], alpha = 0.06, colour = NA) +
      (if (ref == "identity")
         geom_abline(slope = 1, intercept = 0, colour = "#B80000",
                     linetype = "dotdash", linewidth = 1.1)
       else geom_vline(xintercept = 0, colour = "#B80000",
                       linetype = "dotdash", linewidth = 1.1)) +
      geom_line(data = mu, aes(get(xcol), p), colour = "black", linewidth = 1.6) +
      scale_x_continuous(limits = xlim, breaks = pretty(xlim, 4)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = imp_ttl(L, ttl), x = xlab, y = "P(rescue)") + thm
  }
  lim_of = function(ax) c(S$lims[axis == ax, lo], S$lims[axis == ax, hi])
  X_ALN = "Distance from own midpoint"   # the caption gives the full definition
  list(
    raw_stab = pan(FAN_RAW[axis == "stability"], MU_RAW[axis == "stability"], "x", letters[1],
                   quote("P(rescue) along " * (1 - phi)), expression(1 - phi),
                   c(0, 1), "stability", "identity"),
    raw_spec = pan(FAN_RAW[axis == "specificity"], MU_RAW[axis == "specificity"], "x", letters[2],
                   quote("P(rescue) along " * eta), expression(eta),
                   c(0, 1), "specificity", "identity"),
    aln_stab = pan(FAN_ALN[axis == "stability"], MU_ALN[axis == "stability"], "dx", letters[3],
                   quote("P(rescue) along " * (1 - phi) * ", aligned"), X_ALN,
                   lim_of("stability"), "stability", "zero"),
    aln_spec = pan(FAN_ALN[axis == "specificity"], MU_ALN[axis == "specificity"], "dx", letters[4],
                   quote("P(rescue) along " * eta * ", aligned"), X_ALN,
                   lim_of("specificity"), "specificity", "zero"))
}

## the archived version0 layout: raw row over aligned row, lettered A-D
imperfect_treg_fourpanel = function(S, letters = c("A", "B", "C", "D")) {
  P = imp_shape_panels(S, letters)
  (P$raw_stab | P$raw_spec) / (P$aln_stab | P$aln_spec)
}

## the unaligned pair alone, for Extended Data Fig. 6
imperfect_treg_tworow_raw = function(S, letters = c("A", "B")) {
  P = imp_shape_panels(S, c(letters[1], letters[2], "C", "D"))
  P$raw_stab | P$raw_spec
}

## the aligned pair alone, the top row of Fig. 5
imperfect_treg_tworow_aligned = function(S, letters = c("A", "B")) {
  P = imp_shape_panels(S, c("A", "B", letters[1], letters[2]))
  P$aln_stab | P$aln_spec
}

## at each grid value, the fraction of combinations in each band of P(rescue):
## a threshold response holds them at 0 or 1, a graded one fills the middle
imperfect_treg_composition = function(S, letters = c("A", "B")) {
  BS = IMP_BASE_SZ
  B = copy(S$M)[, band := cut(p, c(-1e-9, .001, .25, .5, .75, .999, 1.001),
                              labels = c("0", "(0, .25]", "(.25, .5]",
                                         "(.5, .75]", "(.75, 1)", "1"))]
  CMP = B[, .N, by = .(axis, x, band)][, frac := N / sum(N), by = .(axis, x)]
  FILL = c("0" = "#1d3557", "(0, .25]" = "#457b9d", "(.25, .5]" = "#a8dadc",
           "(.5, .75]" = "#f4a261", "(.75, 1)" = "#e76f51", "1" = "#ffd166")
  pan = function(ax, L, ttl, xlab) ggplot(CMP[axis == ax], aes(x, frac, fill = band)) +
    geom_col(width = 0.085) +
    scale_fill_manual(values = FILL, name = "P(rescue)") +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(title = imp_ttl(L, ttl), x = xlab,
         y = sprintf("fraction of the %d combinations", length(S$units))) +
    theme_minimal(base_size = BS) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = BS * 1.7),
          axis.title = element_text(size = BS * 1.4),
          axis.text = element_text(size = BS * 1.2),
          legend.title = element_text(size = BS * 1.3),
          legend.text = element_text(size = BS * 1.15))
  pan("stability", letters[1], quote("P(rescue) by band, along " * (1 - phi)),
      expression(1 - phi)) |
  pan("specificity", letters[2], quote("P(rescue) by band, along " * eta),
      expression(eta))
}

## ----------------------------------------------------------------------------
## imperfect_treg_examples(S, ids) -- single combinations, to show that the
## midpoint of the response sits in a different place in each. Both marginals
## are drawn against the axis value, with a vertical rule at each curve's own
## midpoint; the panels are ordered by the stability midpoint, so the quantity
## the aligned panels of Fig. 5 remove reads left to right.
## ----------------------------------------------------------------------------
imperfect_treg_examples = function(S, ids, letters = LETTERS[3:6]) {
  BS = IMP_BASE_SZ
  M50 = unique(S$A[, .(unit, axis, x50)])
  ids = M50[axis == "stability" & unit %in% ids][order(x50), unit]
  one = function(u, L) {
    D = S$M[unit == u]; V = M50[unit == u]
    ggplot(D, aes(x, p, colour = axis)) +
      geom_vline(data = V, aes(xintercept = x50, colour = axis),
                 linetype = "dotdash", linewidth = 0.9, show.legend = FALSE) +
      geom_line(linewidth = 1.3) + geom_point(size = 1.6) +
      scale_colour_manual(values = IMP_AX_COL,
                          labels = c(expression(1 - phi), expression(eta)), name = NULL) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
      ## the panel letter alone; which combinations these are is in the legend
      labs(title = bquote(bold(.(L))), x = "axis value", y = "P(rescue)") +
      theme_minimal(base_size = BS) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(size = BS * 1.5),
            axis.title = element_text(size = BS * 1.3),
            axis.text = element_text(size = BS * 1.1),
            legend.position = "none")
  }
  P = Map(one, ids, letters[seq_along(ids)])
  ## the key goes to the right of the last panel, so it lands at the figure's
  ## right-hand edge, under the grey scale of the row above
  P[[length(P)]] = P[[length(P)]] +
    theme(legend.position = "right",
          legend.text = element_text(size = BS * 1.4),
          legend.text.align = 0, legend.key.size = unit(1.6, "lines"))
  Reduce(`|`, P)
}
