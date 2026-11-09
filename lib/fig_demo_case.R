## ============================================================================
## lib/fig_demo_case.R -- the four-arm time courses of one damage-matched unit of
## the gated experiment (the "affordable aggression" demonstration), shared by
## figure_scripts/Fig_4.R (panel D).
## ----------------------------------------------------------------------------
##   fig_demo_case(unit, base_sz)  re-integrates the four arms of unit `unit` (its
##     index i in data/fig04_gated_units.rds) with the experiment's code path on a
##     fine grid, prints the fates and the cumulative damage of each arm over the
##     matching window, and returns list(pA, pB, pC, u, S): the panels P_LP (log),
##     E_h and M1, the unit's row and the summary table. Arms:
##       OFF(beta_r)  ON(beta_r)  ON(beta_r*)  OFF(beta_r*)
##     Dotted vertical line: t_peak, the Tregs-OFF pathogen crest, as in panels B and C;
##     solid grey: t_match (the end of the matching window).
## The caller must have set the working directory and sourced lib/skeleton_defs.R
## and lib/extinction_floor.R (BRANCH = pathogenic_harm).
## ============================================================================
suppressMessages({library(data.table); library(ggplot2)})

fig_demo_case = function(unit = as.integer(Sys.getenv("UNIT", "11773")), base_sz = 26,
                         data = "./data/fig04_damage_matched_gated.rds", tmax = 12000, ke = 100, p_min = 1e-3) {
  fmt = function(z) format(z, big.mark = ",", scientific = FALSE, trim = TRUE)
  trap = function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  TT = c(0, exp(seq(log(0.05), log(tmax), length.out = 3000)))
  u = as.data.table(readRDS(data))[match == "conv" & i == unit]
  if (nrow(u) != 1) stop("unit ", unit, " not in ", data)
  trow = tg[.(u$batch, u$row_ind), nomatch = NULL]
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = u$Kg; p$r_L_P = u$r_L_P; p$r_LP_P = u$r_L_P
  run = function(beta, th, lab) { pp = p; pp$beta_ros = beta
    as.data.table(ode_ext(seed(u$Eh0), TT, rhs_app, c(pp, list(theta_treg = th)), maxsteps = 1e5))[, .(time, PLP, Eh, M1, R)][, run := lab] }
  LV = c("OFF", "ON_same", "ON_star", "OFF_star")
  X = rbindlist(list(run(u$b_off, 0, "OFF"), run(u$b_off, 1, "ON_same"), run(u$b_star, 1, "ON_star"), run(u$b_star, 0, "OFF_star")))[time > 0]
  X[, run := factor(run, levels = LV)]
  S = X[, .(t_clear = time[which(PLP <= EXT_TOL)[1]], P_end = PLP[.N], Eh_end = Eh[.N], M1_end = M1[.N],
            D_win = trap(time[time <= u$win], ke - Eh[time <= u$win])), by = run]
  S[, cost := D_win / S[run == "OFF"]$D_win]
  ## the first peak of the Tregs-OFF pathogen, as in panels B and C and in Extended Data Fig. 5: the first
  ## local maximum, not the maximum over the run (which, for an arm that climbs to a
  ## high-burden attractor, is the attractor)
  off = X[run == "OFF"]; .k = which(diff(sign(diff(off$PLP))) < 0) + 1L
  t_pk = if (length(.k)) off$time[.k[1]] else off$time[which.max(off$PLP)]
  cat(sprintf("[demo] unit %d: tau %d, K_g %g, r_P %.2f, E_h(0) %g, beta_r %.4g -> beta_r* %.4g (%.1fx), t_peak %.3g, t_match %s\n",
              u$i, u$tau_id, u$Kg, u$r_L_P, u$Eh0, u$b_off, u$b_star, u$b_star / u$b_off, t_pk, fmt(round(u$win))))
  print(S[, .(run, t_clear = round(t_clear), P_end = signif(P_end, 3), Eh_end = round(Eh_end, 1), M1_end = round(M1_end, 1), damage_vs_ref = round(cost, 2))])

  ## colour = the beta_r level (dark grey beta_r, pink beta_r*, as in the panels above);
  ## linetype = Tregs (ON solid, OFF dashed, as in the rest of the paper)
  COL = c(OFF = "#555555", ON_same = "#555555", ON_star = "#D81B7A", OFF_star = "#D81B7A")
  LTY = c(OFF = "dashed", ON_same = "solid", ON_star = "solid", OFF_star = "dashed")
  LWD = c(OFF = 1.3, ON_same = 1.3, ON_star = 1.4, OFF_star = 1.4)
  LAB = c(OFF = expression("Tregs-OFF, " * beta[r]), ON_same = expression("Tregs-ON, " * beta[r]),
          ON_star = expression("Tregs-ON, " * beta[r]^"*"), OFF_star = expression("Tregs-OFF, " * beta[r]^"*"))
  th = theme_minimal(base_size = base_sz) +
    theme(panel.grid.minor = element_blank(), legend.title = element_blank(), legend.position = "right",
          legend.text = element_text(size = base_sz * 0.8), legend.key.width = unit(base_sz * 0.1, "cm"),
          legend.key.spacing.y = unit(base_sz * 0.25, "pt"))
  XL = c(min(X$time), tmax)
  base = function(y, ylab, logy = FALSE) {
    D = X; if (logy) D = D[get(y) > p_min]
    g = ggplot(D, aes(time, .data[[y]], colour = run, linetype = run, linewidth = run)) +
      ## t_peak as in panels B and C (grey, dotted); t_match (end of the matching window) solid grey
      ## the time markers carry their own aesthetic (alpha, constant) so they get a key of
      ## their own, separate from the arms'
      geom_vline(data = data.table(x = t_pk, mk = factor("t_peak", levels = c("t_peak", "t_match"))),
                 aes(xintercept = x, alpha = mk), colour = "grey30", linetype = "dotted", linewidth = 1.1, inherit.aes = FALSE, show.legend = c(alpha = TRUE, colour = FALSE, linetype = FALSE, linewidth = FALSE), key_glyph = "vline") +
      geom_vline(data = data.table(x = u$win, mk = factor("t_match", levels = c("t_peak", "t_match"))),
                 aes(xintercept = x, alpha = mk), colour = "grey65", linetype = "solid", linewidth = 1.2, inherit.aes = FALSE, show.legend = c(alpha = TRUE, colour = FALSE, linetype = FALSE, linewidth = FALSE), key_glyph = "vline") +
      scale_alpha_manual(values = c(t_peak = 1, t_match = 1), name = NULL,
                         labels = c(t_peak = expression(t[peak]), t_match = expression(t[match])),
                         guide = guide_legend(order = 2, override.aes = list(linetype = c("dotted", "solid"), colour = c("grey30", "grey65"), linewidth = c(1.1, 1.2)),
                                              theme = theme(legend.text = element_text(size = base_sz * 0.95), legend.key.size = unit(base_sz * 0.07, "cm")))) +
      geom_line() +
      ## the arms' key draws only the geom_line layer (the vline layers would add a tick)
      scale_colour_manual(values = COL, labels = LAB, guide = guide_legend(order = 1, override.aes = list(alpha = 1))) +
      scale_linetype_manual(values = LTY, labels = LAB, guide = guide_legend(order = 1)) +
      scale_linewidth_manual(values = LWD, labels = LAB, guide = guide_legend(order = 1)) +
      scale_x_log10(breaks = c(0.1, 10, 1000), labels = function(z) sapply(z, format)) +
      labs(x = "Time", y = ylab) + th
    if (logy) g + scale_y_log10(labels = function(z) sapply(z, format)) + coord_cartesian(xlim = XL, ylim = c(p_min, NA))
    else g + coord_cartesian(xlim = XL)
  }
  ## A: the eliminated runs labelled where they leave the axis, the persistent runs' final loads at the right
  TX = X[run %in% S[is.finite(t_clear)]$run & PLP > p_min, .(t = max(time)), by = run]
  ## labels left of the drop (the panel may be narrow), stacked at the bottom of the axis
  TX = merge(TX, S[, .(run, t_clear)], by = "run"); TX[, y := p_min * c(3.2, 1)[seq_len(.N)]]
  EP = S[!is.finite(t_clear)]
  NM = c(OFF = "Tregs-OFF", ON_same = "Tregs-ON", ON_star = "Tregs-ON", OFF_star = "Tregs-OFF")
  ## the final load of every arm at the far right, level with its line; the eliminated arms
  ## read 0 at the bottom of the axis
  ## labels whose values coincide (the two eliminated arms at 0; two arms ending at the same
  ## E_h) are offset vertically so both stay legible
  spread = function(v, vj = rep(0.5, length(v)), tol) { o = order(v); d = c(Inf, diff(v[o])); g = cumsum(d > tol)
    for (k in unique(g)) { idx = o[g == k]; if (length(idx) > 1) vj[idx] = seq(1.3, -0.3, length.out = length(idx)) }
    vj }
  EN = copy(S)[, y := fifelse(is.finite(t_clear), p_min * 1.6, P_end)][, lab := fifelse(is.finite(t_clear), "0", sprintf("%.0f", P_end))]
  EN[, vj := spread(log10(y), tol = 0.12)]
  pA = base("PLP", expression(P[LP]), logy = TRUE) +
    geom_text(data = EN, aes(x = tmax * 1.6, y = y, label = lab, vjust = vj), inherit.aes = FALSE,
              colour = COL[as.character(EN$run)], size = base_sz * 0.28, hjust = 0, fontface = "bold") +
    coord_cartesian(xlim = c(min(X$time), tmax * 2.2), ylim = c(p_min, NA), clip = "off")
  endlab = function(g, y) { E2 = S[, .(run, v = get(y))]; E2[, vj := spread(v, tol = 0.04 * diff(range(X[[sub("_end", "", y)]])))]
    g + geom_text(data = E2, aes(x = tmax * 1.6, y = v, label = sprintf("%.0f", v), vjust = vj), inherit.aes = FALSE,
                  colour = COL[as.character(E2$run)], size = base_sz * 0.28, hjust = 0, fontface = "bold") +
      coord_cartesian(xlim = c(min(X$time), tmax * 2.2), clip = "off") }
  list(pA = pA, pB = endlab(base("Eh", expression(E[h])), "Eh_end"), pC = endlab(base("M1", expression(M[1])), "M1_end"), u = u, S = S)
}
