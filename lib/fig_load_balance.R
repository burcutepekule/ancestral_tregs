## ============================================================================
## lib/fig_load_balance.R -- drawing code shared by Fig. 3C and Extended Data Fig. 2C (the
## bistable class) and Extended Data Fig. 3 (the ON_mono class).
##
## One figure per class: rows sterile (A, B) and pathogenic help (C, D), columns
## Tregs OFF and Tregs ON. N_SHOW combinations per class and branch, drawn at random
## (seed 1), the same combinations in both columns, each integrated from the midpoint
## of its rescued interval (pipeline/20_load_balance.R). Paths in the plane
## (Ea + En)/Kg (log) against dL/dt = killing - clearance; circle = initial condition,
## triangle = last stored point (filled if at rest at the horizon). Sterile paths and
## markers coloured by M1 (linear, 0-100); pathogenic by PLP/K_LP,P (log, 1e-4 to 1,
## floored at 1e-4); one scale per branch, shared by Tregs OFF and ON. Shaded:
## L > Kg (overwhelmed efferocytosis). Points with L/Kg <= X_FLOOR are not drawn.
##
## in : ./data/load_balance_sterile.rds, ./data/load_balance_pathogenic_help.rds
## use: source("./lib/fig_load_balance.R"); g = fig_load_balance("Bistable")   # or "ON_mono"
##      grid::grid.draw(g)  on a 4300 x 3700 px device at res 200
## `branches` and `letters` draw ONE row on its own, for the merged figures of the
## The overlaid form (figure_scripts/Fig_3.R, figure_scripts/Extended_Data_Fig_2.R):
##      fig_load_balance("Bistable", branches = "Sterile", letters = c("C", "D"))
## Everything up to the row grobs is computed over BOTH branches whatever is asked
## for, so the shared y range YL and the per-panel x ranges, grids and region
## arrows of a one-row figure are the ones it would have had inside the two-row
## figure; `branches` only selects which rows are drawn and stacked. The marker
## key goes on the LAST drawn row, so a one-row figure still carries it.
## ============================================================================
suppressMessages({library(data.table); library(ggplot2); library(ggnewscale); library(ggh4x); library(gtable)})

## The device every OVERLAID panel is rendered on, defined once and used by all of
## them. One overlaid row is 4300 x 1850 at res 200 -- 21.5 x 9.25
## physical inches -- rescaled to LB_W px wide, res and height together, so the
## panel keeps that physical size, and therefore its type, markers and legend keep
## theirs, whatever the pixel width. A figure is n_rows of these.
## n_rows is the ONLY thing that differs between the three figures that carry
## overlaid panels, and it differs only because Extended Data Fig. 3 has no
## phase-portrait half above its panels while Fig. 3 and Extended Data Fig. 2 do.
LB_W = 6400
lb_png = function(path, n_rows = 1L) {
  sc = LB_W / 4300
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  png(path, width = LB_W, height = round(n_rows * 1850 * sc), res = 200 * sc)
}

## One effect class in both rows. The pathogenic branch is the Treg-help class by
## construction; the sterile gate of pipeline/20_load_balance.R keeps every bistable
## combination with a rescue gap of at least 5, help and neutral alike, so the sterile
## combinations are restricted to help here. pipeline/24_numbers_load_balance_boundary.R applies the
## same filter, and the two must stay in step.
HELP_ONLY = function(D, br) {
  EFF = unique(as.data.table(readRDS("./data/selected_big.rds"))[, .(tau_id, Kg, beta_ros, r_L_P, treg_effect)])
  k = function(X) paste(X$tau_id, X$Kg, X$beta_ros, X$r_L_P)
  D = copy(D)[, tau_id := batch * 10000L + row_ind]
  if (!all(k(D) %in% k(EFF))) stop("[fig_load_balance] a combination is missing from selected_big.rds")
  out = D[k(D) %in% k(EFF[treg_effect == "help"])]
  cat(sprintf("[fig_load_balance] %-16s Treg-help filter: %d of %d combination-condition rows kept\n",
              br, nrow(out), nrow(D)))
  out }

fig_load_balance = function(CLS, n_show = as.integer(Sys.getenv("N_SHOW", "20")),
                            branches = c("Sterile", "Pathogenic"), letters = NULL,
                            overlay = FALSE, short_ylab = overlay || length(branches) < 2L) {
  branches = match.arg(branches, c("Sterile", "Pathogenic"), several.ok = TRUE)
  P_FLOOR = 1e-4; X_FLOOR = 1e-3
  COL_OFF = COL_ON = viridisLite::viridis(6, direction = -1, end = 0.95)
  hex_of = function(v, cols) { r = colorRamp(cols)((v - log10(P_FLOOR)) / (0 - log10(P_FLOOR))); rgb(r[, 1], r[, 2], r[, 3], maxColorValue = 255) }
  contour_of = function(hex) { r = col2rgb(hex) / 255; ifelse(0.299 * r[1, ] + 0.587 * r[2, ] + 0.114 * r[3, ] < 0.5, "white", "black") }
  ## overlay puts the two Tregs conditions of a branch in ONE panel, so the panel
  ## key loses its OFF/ON half and everything downstream -- XR, REG, GX, the facet
  ## -- collapses to one panel per branch with no further special-casing
  PANELS = if (overlay) c("Sterile", "Pathogenic") else c("Sterile OFF", "Sterile ON", "Pathogenic OFF", "Pathogenic ON")
  pkey = function(branch, cond) if (overlay) branch else paste(branch, ifelse(cond == "Tregs OFF", "OFF", "ON"))
  ## dL/dt = killing - clearance exactly (M0 + M1 conserved)
  YLAB = expression(atop("Epithelial killing - clearance  (cells per time unit)", d*(E[a]+E[n])/dt == E[h]*(tau[r]*R + tau[P*","*B]*P[LP]/K[LP*","*P]) - g(E[a]+E[n])*(M[0]+M[1])*(E[a]+E[n])))
  ## Set on one line that label is 10.5 in long and the equation under it 15.1 in;
  ## 15.1 in is what two stacked rows of panels give it to span. A one-row figure
  ## gives about 7.1 in, so it drops the equation -- which the caption states
  ## anyway -- and breaks the words over two lines, 6.0 in at the same 36 pt.
  ## Nothing is scaled down: both labels are set at YLAB_FS. `short_ylab` defaults
  ## to the short form wherever the equation cannot fit or the figure belongs to the
  ## merged set, and is an argument so a caller can ask for either.
  YLAB_1ROW = expression(atop("Epithelial killing - clearance", "(cells per time unit)"))
  YLAB_FS = 36
  CLS_OF = function(x) factor(ifelse(x == "mono_B", "ON_mono", "Bistable"), levels = c("Bistable", "ON_mono"))

  load_branch = function(br, lab) {
    D = HELP_ONLY(as.data.table(readRDS(sprintf("./data/load_balance_%s.rds", br))), br)
    P = D[, rbindlist(lapply(seq_len(.N), function(j) data.table(cid = cid[j], cls = cls[j], cond = cond[j], pop[[j]])))]
    P[, cls := CLS_OF(cls)]
    P[, `:=`(branch = lab, x = L_over_Kg, y = kill_minus_clear, plog = log10(pmax(PLP_over_K, P_FLOOR)))]
    P = P[is.finite(x) & x > X_FLOOR & is.finite(y)]
    ## the drawn combinations: n_show per class, seed 1, drawn from both classes in order
    set.seed(1); show = P[, .(cid = sample(unique(cid), min(n_show, uniqueN(cid)))), by = cls]
    IND = P[show, on = .(cls, cid)][order(cid, cond, time)]
    ST = IND[time == 0]; EN = IND[, .SD[which.max(time)], by = .(cid, cls, cond)]
    EN = merge(EN, D[, .(cid, cls = CLS_OF(cls), cond, steady_end)], by = c("cid", "cls", "cond"))
    ## shp is the triangle the Tregs-OFF end marker uses everywhere; shp_on is the
    ## circle the Tregs-ON end marker uses under overlay = TRUE (see row_of). Both
    ## keep the filled/open distinction: filled when at rest at the horizon.
    EN[, `:=`(shp = ifelse(steady_end, 24L, 2L), shp_on = ifelse(steady_end, 21L, 1L))]; ST[, shp := 21L]
    list(IND = IND, ST = ST, EN = EN) }
  S = load_branch("sterile", "Sterile"); H = load_branch("pathogenic_help", "Pathogenic")
  IND_ALL = rbind(S$IND, H$IND); MK_ALL = rbind(S$ST, S$EN, H$ST, H$EN, fill = TRUE)
  IND_ALL[, panel := factor(pkey(branch, cond), levels = PANELS)]
  MK_ALL[, panel := factor(pkey(branch, cond), levels = PANELS)]
  ## marker fills: pathogenic by pathogen load, sterile by M1 (linear over the pool of 100)
  MK_ALL[, fillhex := ifelse(branch == "Sterile", NA_character_, ifelse(cond == "Tregs OFF", hex_of(plog, COL_OFF), hex_of(plog, COL_ON)))]
  MK_ALL[, contour := ifelse(branch == "Sterile", NA_character_, contour_of(fillhex))]
  brks = log10(c(1e-4, 1e-3, 1e-2, 1e-1, 1)); lbls = c("≤ 1e-4", "1e-3", "1e-2", "1e-1", "1")
  GY = data.table(y = seq(-10, 10, 5))
  ## These are the published values and they stay. The ramp encodes a variable; the
  ## panel shading is decoration, so if the pale end of the ramp is hard to see on
  ## the shading it is the shading that has to give way, not the scale. See the note
  ## on the shaded rect in row_of().
  COL_M1 = c("#fde0dd", "#fa9fb5", "#f768a1", "#c51b8a", "#7a0177", "#49006a")   # light pink to dark purple
  M1_MAX = 100
  hex_m1 = function(v) { r = colorRamp(COL_M1)(pmin(pmax(v, 0), M1_MAX) / M1_MAX); rgb(r[, 1], r[, 2], r[, 3], maxColorValue = 255) }
  MK_ALL[branch == "Sterile", fillhex := hex_m1(M1)]
  MK_ALL[branch == "Sterile", contour := contour_of(fillhex)]
  LWD_PATH = 1.2
  ## overlay marker sizes. The end markers are larger than the 6.4 of the separate
  ## panels: one overlaid panel carries twice as many of them, and the shape is now
  ## what tells the two arms apart, so it has to be legible at a glance.
  ## overlay marker hierarchy. The end markers are larger than the 6.4 of the
  ## separate panels and carry a heavy contour: one overlaid panel holds twice as
  ## many of them and the shape is now what tells the two arms apart, so it has to
  ## read at a glance. The initial-condition cross sits clearly below them on both
  ## counts -- smaller, and a contour thin enough not to compete.
  ## On an unfilled glyph a "contour" is not a rim around a fill: it is a thicker
  ## dark cross drawn under a thinner coloured one, so what the eye weighs is
  ## MK_IC_HALO, not MK_IC_HALO - MK_IC_STROKE. The cross therefore has to come in
  ## under MK_STROKE on that number, not merely on the difference, or a cluster of
  ## twenty overlapping crosses reads as a dark knot.
  MK_END = 8.4; MK_STROKE = 1.7          # steady-state circle and triangle
  MK_IC  = 4.8; MK_IC_STROKE = 1.0       # initial-condition cross: 57% of MK_END
  MK_IC_HALO = 1.5                       # its contour: under MK_STROKE, as asked

  ## ---- this class ----------------------------------------------------------
  IND = IND_ALL[cls == CLS]; MK = MK_ALL[cls == CLS]
  XR = IND[, .(xmin = min(x), xmax = max(x)), by = panel]; xmax_of = function(p) XR[match(p, panel), xmax]
  LAB = data.table(panel = factor(PANELS, levels = PANELS))
  LAB[, letter := LETTERS[seq_len(.N)]]
  ## DRAWN is a subset of PANELS in PANELS order, so `letters` lines up row by row
  DRAWN = PANELS[sub(" .*", "", PANELS) %in% branches]
  if (!is.null(letters)) {
    if (length(letters) != length(DRAWN))
      stop(sprintf("[fig_load_balance] %d letters for %d drawn panels", length(letters), length(DRAWN)))
    LAB[panel %in% DRAWN, letter := letters]
  }
  ## overlay: one title for two arms, so it names the branch and class only and
  ## the arms are carried by the line-type and end-marker keys
  LAB[, text := if (overlay) sprintf("%s   %s, %s", letter, panel, ifelse(CLS == "Bistable", "bistable", "monostable"))
                else sprintf("%s   %s, %s, Tregs-%s", letter, sub(" .*", "", panel), ifelse(CLS == "Bistable", "bistable", "monostable"), sub(".* ", "", panel))]
  LAB[, title := factor(text, levels = text)]   # facet strip titles, A-D in order
  IND = merge(IND, LAB[, .(panel, title)], by = "panel"); MK = merge(MK, LAB[, .(panel, title)], by = "panel")
  ## The shared y range: computed over BOTH branches, so every panel of a class gets
  ## the same one. This is what the published figures do and it is left alone.
  yr = range(IND$y); pad = 0.04 * diff(yr); YL = c(yr[1] - pad, yr[2] + 8 * pad)

  ## Under overlay the range is computed PER PANEL instead, by yl_of() below, from
  ## that panel's own trajectories. The shared range gave three of the four overlaid
  ## panels headroom driven by data they do not contain -- the sterile panel peaks at
  ## about +3.5 but inherited a top of +10.3 from the pathogenic one -- and no
  ## multiplier fixes that, because the baseline it multiplies is the wrong number.
  ## The range is then built downward from the top of the panel, out of quantities
  ## measured off the render rather than guessed:
  ##   TOP_M   margin above the label block
  ##   LBL_H   the label block itself, two lines at size 10
  ##   LBL_GAP the gap between the label's bottom and its arrow
  ## so that the label's BOTTOM lands CLEAR data units above the highest path, and
  ## YARROW_OV is where the arrow then falls as a fraction of the finished range.
  ## The bottom keeps the published 0.04 pad.
  TOP_M = 0.02; LBL_H = 0.117; LBL_GAP = 0.024; CLEAR = 0.75
  ## TOP_M has two jobs that are now separated. yl_of() still uses TOP_M to SIZE
  ## the range, so YLp and the y breaks do not move. TOP_M_DRAW only PLACES the
  ## band inside that range, so shrinking it lifts the arrows and their labels
  ## toward the top of the panel and away from the trajectory peak, by exactly
  ## (TOP_M - TOP_M_DRAW) * diff(YLp). The axis is untouched by design.
  TOP_M_DRAW = 0.004
  YARROW_OV = TOP_M_DRAW + LBL_H + LBL_GAP
  yl_of = function(I) {
    r = range(I$y); span = diff(r)
    lo = r[1] - 0.04 * span
    c(lo, lo + (1.04 * span + CLEAR) / (1 - TOP_M - LBL_H))
  }
  REG = rbind(
    XR[, .(panel, side = "right", x0 = 1.4, x1 = pmin(xmax * 0.9, 22), xt = pmin(xmax * 0.95, 24), hj = 1, label = "Overwhelmed\nefferocytosis")],
    XR[xmin < 0.3, .(panel, side = "left", x0 = 0.72, x1 = pmax(xmin * 1.5, 0.04), xt = 0.72, hj = 1, label = "Functional\nefferocytosis")])
  REG[, vj := -0.35]
  ## Without overlay, a panel carrying BOTH labels drops the right one under its
  ## arrow to keep the two apart. Under overlay they both sit above their arrows --
  ## the two are far apart on the x axis anyway -- so the rule is skipped and the
  ## annotation collapses into one band at the top of the panel.
  if (!overlay)
    REG[side == "right" & panel %in% REG[side == "left", panel], `:=`(xt = xmax_of(panel), hj = 1, vj = 1.4)]
  REG = merge(REG, LAB[, .(panel, title)], by = "panel")
  ## grid drawn as layers: above the shaded rect, below the paths and markers, inside the data range only
  ## overlay spans the union of the two arms' ranges -- about 0.018 to 20, three
  ## decades -- over which breaks_log() returns decades only (0.1, 1, 10 inside the
  ## range). 1-and-3 breaks keep the resolution the separate panels had and still
  ## put a label on the x = 1 efferocytosis boundary; the axis uses the same set,
  ## so the labels sit on the grid lines.
  XBRK = if (overlay) sort(c(10^(-4:4), 3 * 10^(-4:4))) else NULL
  ## and one label format for all of them. scale_x_log10() picks its own format from
  ## the range, so a panel spanning 4.5 decades (sterile ON_mono) came out as
  ## 1e-03, 3e-03, ... while one spanning 3 came out as 0.03, 0.10, ... -- the same
  ## breaks labelled two different ways in panels that have to match.
  XBRK_LAB = function(b) trimws(formatC(b, format = "fg", digits = 3))
  GX = XR[, .(x = if (overlay) XBRK else scales::breaks_log()(c(xmin, xmax))), by = .(panel, xmin, xmax)][x >= xmin & x <= xmax, .(panel, x)]
  GX = merge(GX, LAB[, .(panel, title)], by = "panel"); XR = merge(XR, LAB[, .(panel, title)], by = "panel")

  ## grid cannot dash a line whose colour varies along it, and every path here is
  ## coloured by M1 or by P_LP, which change as the trajectory runs. So the dashes
  ## of the Tregs-OFF arm are cut into the data instead of asked of the renderer:
  ## arc length is accumulated in panel-normalised coordinates (log10 x over the
  ## panel's x span, y over YL, x weighted by DASH_ASP because the overlay panel is
  ## about twice as wide as it is tall), and alternate runs of DASH_LEN are dropped.
  ## The gaps are then even on the page whatever the spacing of the stored points.
  DASH_LEN = 0.011; DASH_ASP = 2.0; DASH_STEP = DASH_LEN / 8
  dash_cut = function(D, xr, yl) {
    D = copy(D)[order(cid, cond, time)]
    lx0 = log10(xr[1]); lxs = diff(log10(xr)); ys0 = yl[1]; ysp = diff(yl)
    ## Resample each path at a fixed step in panel-normalised arc length BEFORE
    ## cutting it into dashes. Without this the run index is assigned to the STORED
    ## vertices, so wherever the trajectory moves fast and two stored points are more
    ## than one DASH_LEN apart, a single segment spans several runs: the pattern
    ## degenerates into an over-long straight dash or a missing chunk. A solid line
    ## never shows this, which is why the same path looks smooth undashed.
    VARS = intersect(c("M1", "plog", "time"), names(D))
    D = D[, {
      sx = DASH_ASP * (log10(x) - lx0) / lxs; sy = (y - ys0) / ysp
      s  = cumsum(c(0, sqrt(diff(sx)^2 + diff(sy)^2)))
      k  = c(TRUE, diff(s) > 0)          # approx() needs a strictly increasing index
      if (sum(k) < 2L) {
        out = list(s = 0, x = x[1], y = y[1])
        for (v in VARS) out[[v]] = .SD[[v]][1]
      } else {
        sk = s[k]; tgt = sort(unique(c(sk, seq(0, max(sk), by = DASH_STEP))))
        out = list(s = tgt,
                   x = 10^(lx0 + approx(sk, sx[k], tgt)$y / DASH_ASP * lxs),
                   y = ys0 + approx(sk, sy[k], tgt)$y * ysp)
        for (v in VARS) out[[v]] = approx(sk, .SD[[v]][k], tgt)$y
      }
      as.data.table(out)
    }, by = .(cid, cond)]
    D[, run := floor(s / DASH_LEN)]
    ## each kept run carries the first point of the run after it, so a dash ends
    ## where the gap begins instead of a point short of it
    D = rbind(D[run %% 2L == 0L], D[run %% 2L == 1L, .SD[1], by = .(cid, cond, run)][, run := run - 1L], fill = TRUE)
    D[order(cid, cond, s)][, grp := paste(cid, cond, run)][] }

  row_of = function(br) {
    I = IND[branch == br]; M = MK[branch == br]; R2 = REG[title %in% I$title]; G2 = GX[title %in% I$title]; X2 = XR[title %in% I$title]
    ## per-panel y range under overlay; R2 is a fresh data.table, so re-placing the
    ## arrows on it does not touch REG or the other panel
    YLp = if (overlay) yl_of(I) else YL
    if (overlay) {
      R2[, yl := YLp[2] - YARROW_OV * diff(YLp)]
      bk = seq(-10, 10, 5); bk = bk[bk >= YLp[1] & bk <= YLp[2]]
      cat(sprintf("[fig_load_balance] %-10s %-8s peak %6.2f | YL [%6.2f, %6.2f] | arrow %5.2f | label bottom %5.2f | clearance %4.2f | y breaks %s\n",
                  br, CLS, max(I$y), YLp[1], YLp[2], R2$yl[1], R2$yl[1] + LBL_GAP * diff(YLp),
                  R2$yl[1] + LBL_GAP * diff(YLp) - max(I$y), paste(bk, collapse = ",")))
    }
    g = ggplot() +
      geom_blank(data = X2[, .(x = c(xmin, xmax)), by = title], aes(x = x)) +
      ## The fill was #EAECF0. A thin path reads against a background by LIGHTNESS,
      ## and in CIE L* that fill is 93.4 against 91.4 for the pale end of the M1 ramp
      ## -- a margin of 1.9, which is why the first stretch of every trajectory, where
      ## M1 is still near 0, was invisible inside the shading and perfectly legible on
      ## the white half of the same panel (white is L* 100, a margin of 8.6).
      ## The ramp encodes a variable and this fill is decoration, so the fill gives
      ## way. It cannot go DARKER: inside this region the drawn points occupy L* 91.4
      ## down to 21.1 continuously, so any fill in that band collides with some path.
      ## Lighter is the only direction, and it trades against how much region is left:
      ## #F3F5F8 is L* 96.5, a margin of 5.0 from the pale ramp end and 3.5 from white.
      ## The x = 1 rule and the two labelled arrows carry the boundary regardless.
      annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F3F5F8", alpha = 1) +
      geom_vline(data = G2, aes(xintercept = x), colour = "grey82", linewidth = 0.4) +
      geom_hline(data = GY, aes(yintercept = y), colour = "grey82", linewidth = 0.4) +
      geom_vline(xintercept = 1, colour = "grey40", linewidth = 0.4) +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4, linetype = "33")
    if (overlay) {
      ## one panel, both arms. Colour still carries M1 (sterile) or P_LP
      ## (pathogenic); the arm is carried by line type, solid ON / dashed OFF, the
      ## convention Fig. 2 and Figs. 3A, 3B already use. The OFF arm is drawn
      ## solid from pre-cut dashes (see dash_cut); the line-type scale exists for
      ## the legend only, which is why both its values are drawn, not mapped.
      IOFF = dash_cut(I[cond == "Tregs OFF"], c(X2$xmin[1], X2$xmax[1]), YLp)
      ION  = I[cond == "Tregs ON"]
      cvar = if (br == "Sterile") quote(M1) else quote(plog)
      g = g + geom_path(data = IOFF, aes(x, y, group = grp, colour = !!cvar), linewidth = LWD_PATH, alpha = 0.9) +
        geom_path(data = ION, aes(x, y, group = cid, colour = !!cvar), linewidth = LWD_PATH, alpha = 0.9) +
        (if (br == "Sterile")
           scale_colour_gradientn(colours = COL_M1, limits = c(0, M1_MAX), breaks = c(0, 25, 50, 75, 100), name = expression(M[1]), guide = guide_colourbar(order = 1, barheight = unit(5.5, "cm")))
         else
           scale_colour_gradientn(colours = COL_OFF, limits = log10(c(P_FLOOR, 1)), breaks = brks, labels = lbls, name = expression(P[LP] / K[LP*","*P]), guide = guide_colourbar(order = 1, barheight = unit(5.5, "cm")))) +
        new_scale_colour()
    } else if (br == "Sterile") {
      g = g + geom_path(data = I, aes(x, y, group = interaction(cid, cond), colour = M1), linewidth = LWD_PATH, alpha = 0.9) +
        scale_colour_gradientn(colours = COL_M1, limits = c(0, M1_MAX), breaks = c(0, 25, 50, 75, 100), name = expression(M[1]), guide = guide_colourbar(order = 1, barheight = unit(5.5, "cm"))) +
        new_scale_colour()
    } else {
      g = g + geom_path(data = I[cond == "Tregs OFF"], aes(x, y, group = cid, colour = plog), linewidth = LWD_PATH, alpha = 0.9) +
        scale_colour_gradientn(colours = COL_OFF, limits = log10(c(P_FLOOR, 1)), breaks = brks, labels = lbls, name = expression(P[LP] / K[LP*","*P]), guide = guide_colourbar(order = 1, barheight = unit(5.5, "cm"))) +
        new_scale_colour() +
        geom_path(data = I[cond == "Tregs ON"], aes(x, y, group = cid, colour = plog), linewidth = LWD_PATH, alpha = 0.9) +
        ## the ON paths: a second layer on the SAME ramp, so a scale of their own (ggnewscale) but no second bar
        scale_colour_gradientn(colours = COL_ON, limits = log10(c(P_FLOOR, 1)), breaks = brks, labels = lbls, guide = "none") +
        new_scale_colour()
    }
    if (overlay) {
      ## The two arms of a pair start from the same state, verified bitwise equal in
      ## x and y for all 20 pairs, so the start marker is drawn once, off the OFF arm.
      ## It is a cross, pch 4, the glyph the initial conditions already carry on the
      ## E_h axis of the phase portraits, and it is coloured by the panel's own scale
      ## like everything else in the panel. A cross has no fill and no outline, so it
      ## cannot get its contrast from contour_of() the way the filled markers do; it
      ## has no fill and no outline of its own, so it gets its contour the way the
      ## filled markers do -- contour_of() -- by being drawn TWICE: a slightly larger,
      ## heavier cross in the contour colour, then the scale colour over it. The rim
      ## that leaves is (MK_IC_HALO - MK_IC_STROKE)/2, deliberately well under the
      ## MK_STROKE the steady-state markers carry, so the two do not read as equals.
      ## The end markers keep the fill and the filled/open "at rest" distinction and
      ## carry the arm by shape: triangle Tregs-OFF, circle Tregs-ON.
      g = g + geom_point(data = M[shp != 21L & cond == "Tregs OFF"], aes(x, y, fill = I(fillhex), colour = I(contour), shape = I(shp)), size = MK_END, stroke = MK_STROKE) +
        geom_point(data = M[shp != 21L & cond == "Tregs ON"], aes(x, y, fill = I(fillhex), colour = I(contour), shape = I(shp_on)), size = MK_END, stroke = MK_STROKE) +
        geom_point(data = M[shp == 21L & cond == "Tregs OFF"], aes(x, y, colour = I(contour)), shape = 4, size = MK_IC + 0.3, stroke = MK_IC_HALO) +
        geom_point(data = M[shp == 21L & cond == "Tregs OFF"], aes(x, y, colour = I(fillhex)), shape = 4, size = MK_IC, stroke = MK_IC_STROKE)
    } else {
      g = g + geom_point(data = M[shp == 21L], aes(x, y, fill = I(fillhex), colour = I(contour)), shape = 21, size = 4.8, stroke = 0.9) +
        geom_point(data = M[shp != 21L], aes(x, y, fill = I(fillhex), colour = I(contour), shape = I(shp)), size = 6.4, stroke = 0.9)
    }
    if (overlay) {
      ## Both keys on EVERY overlaid panel, as layers with no drawable points. In the
      ## two-row ON_mono figure each row already carries its own colour bar, because
      ## the two branches are on different scales, so giving each row the line-type
      ## and marker keys too is what makes all four overlaid panels in the set
      ## carry an identical legend rather than one of them carrying a partial one.
      ## The line-type key box has to be wide for a dash pattern to read at all at
      ## 30 pt, but the width is set ON THIS GUIDE, not in theme(): guide_colourbar()
      ## sets barheight and leaves barwidth unset, and an unset barwidth falls back
      ## to legend.key.width -- so widening it in theme() widened the colour bar from
      ## 0.61 cm to 2.61 cm. Scoped here, the bar keeps its usual width.
      MKL = c("Initial condition, both arms", "Steady state, Tregs OFF", "Steady state, Tregs ON")
      g = g + geom_path(data = data.table(x = NA_real_, y = NA_real_,
                                          lt = factor(rep(c("Tregs OFF", "Tregs ON"), each = 2), levels = c("Tregs OFF", "Tregs ON"))),
                        aes(x, y, linetype = lt), na.rm = TRUE) +
        scale_linetype_manual(values = c("Tregs OFF" = "22", "Tregs ON" = "solid"), name = NULL,
                              guide = guide_legend(order = 2, theme = theme(legend.key.width = unit(2.6, "cm")),
                                                   override.aes = list(colour = "grey25", linewidth = 2.0))) +
        geom_point(data = data.table(x = NA_real_, y = NA_real_, mk = factor(MKL, levels = MKL)), aes(x, y, shape = mk), na.rm = TRUE) +
        scale_shape_manual(values = setNames(c(4L, 24L, 21L), MKL), name = NULL,
                           guide = guide_legend(order = 3, override.aes = list(fill = "white", colour = "black", size = 8, stroke = 1.2)))
    } else if (!overlay && br == branches[length(branches)])   # marker key once, on the bottom drawn row: a layer with no drawable points
      g = g + geom_point(data = data.table(x = NA_real_, y = NA_real_, mk = factor(c("Initial condition", "Steady state"), levels = c("Initial condition", "Steady state"))), aes(x, y, shape = mk), na.rm = TRUE) +
        scale_shape_manual(values = c("Initial condition" = 21, "Steady state" = 24), name = NULL,
                           guide = guide_legend(order = 3, override.aes = list(fill = "white", colour = "black", size = 7, stroke = 1)))
    g = g + geom_segment(data = R2, aes(x = x0, xend = x1, y = yl, yend = yl), colour = "grey30", linewidth = 1.1, arrow = arrow(length = unit(0.4, "cm"), type = "closed")) +
      geom_text(data = R2, aes(x = xt, y = yl, label = label, hjust = hj, vjust = vj), size = 10, colour = "grey30", lineheight = 0.9) +
      facet_wrap2(~ title, ncol = 2, scales = "free_x", axes = "margins") +
      (if (overlay) scale_x_log10(breaks = XBRK, labels = XBRK_LAB) else scale_x_log10()) + scale_y_continuous(breaks = seq(-10, 10, 5)) + coord_cartesian(ylim = YLp, clip = "off") +
      labs(x = expression((E[a] + E[n]) / K[g]), y = NULL) +
      theme_minimal(base_size = 30) + theme(legend.position = "right", panel.grid.minor = element_blank(), strip.text = element_text(hjust = 0, face = "bold", size = 30, colour = "grey15", lineheight = 0.95, margin = margin(b = 8)), strip.clip = "off",
                                            panel.spacing.x = unit(4, "lines"),   # the two columns have different x axes
                                            panel.grid.major = element_blank(),   # grid is drawn as layers above (GX, GY)
                                            legend.text = element_text(size = 30), legend.key.spacing.y = unit(0.4, "lines"), legend.title = element_text(size = 32), axis.title = element_text(size = 36), axis.text = element_text(size = 30),
                                            plot.margin = margin(t = 20, r = 10, b = 20, l = 10))
    g }
  ## the arrows' y position, read by row_of through aes(y = yl)
  REG[, yl := YL[2] - 0.17 * diff(YL)]
  GR = lapply(branches, function(br) ggplotGrob(row_of(br)))
  if (length(unique(vapply(GR, function(x) length(x$widths), integer(1)))) > 1L)
    stop("[fig_load_balance] the rows do not share a column layout; cannot stack them")
  g = Reduce(function(a, b) rbind(a, b, size = "max"), GR)
  pan = g$layout[grepl("^panel", g$layout$name), ]
  g = gtable::gtable_add_cols(g, grid::unit(3.6, "cm"), pos = 0)   # two lines at 36 pt, either label
  gtable::gtable_add_grob(g, grid::textGrob(if (short_ylab) YLAB_1ROW else YLAB,
                                            rot = 90, gp = grid::gpar(fontsize = YLAB_FS)),
                          t = min(pan$t), b = max(pan$b), l = 1, clip = "off")
}
