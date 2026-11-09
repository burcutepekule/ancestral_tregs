## ============================================================================
## lib/fig_skeleton.R -- the two-forms-of-rescue figure, one page per branch:
## Figs. 3A, 3B and Extended Data Figs. 2A, 2B (sterile and pathogenic, Tregs help),
## Extended Data Fig. 1 (the attractor shifts and rescue gaps behind them)
## and Supplementary Fig. 2 (pathogenic, Tregs harm).
## ----------------------------------------------------------------------------
## Layout (one block per Tregs-ON class; two blocks unless mono_A is non-empty):
##   A  phase portrait of the bistable combination with the largest rescue gap
##   B  attractor shift, Tregs ON - Tregs OFF, bistable class (violins)
##   C  rescue gap of the bistable class, on the bisection lattice
##   D  phase portrait of the ON_mono combination with the largest drop in the
##      low-burden composite
##   E  attractor shift, ON_mono class        F  left blank (no resolved gap in ON_mono)
##   G-I the same for mono_A, drawn only if that class is non-empty (it is empty in
##      all three branches of the stored selection)
## The sterile page draws the portraits in (Eh, M1) -- P_LP is identically 0 there;
## the pathogenic pages in (Eh, P_LP), the plane of the composite ruler.
## Language: LOW-BURDEN / HIGH-BURDEN attractor, as Fig. 2 (figure_scripts/Fig_2.R) words
## it. Tregs ON solid, Tregs OFF dashed + darker.
##
## The drawing code is the earlier per-branch scripts' code, unchanged; the only
## differences between branches are the portrait plane, the P_LP violins and the
## panel-D title, exactly as those scripts had them.
##
## in : ./data/skeleton_<BRANCH>.rds, ./data/rescue_<BRANCH>.rds, ./data/tau_grid_v6_*.rds
## ============================================================================
fig_skeleton = function(BRANCH, out_png, effect = NULL, panels = NULL) {
## lib/skeleton_defs.R reads BRANCH from the global environment (PL0, PLP0)
assign("BRANCH", BRANCH, envir = globalenv())
suppressMessages({library(data.table); library(deSolve)})
MONO_LAB = Sys.getenv("MONO_LAB", "ON_mono")
source("./lib/skeleton_defs.R")   # into the global environment, as the figure scripts always did
HMAX = as.numeric(Sys.getenv("HMAX", "5"))
CK = c("tau_id","Kg","beta_ros","r_L_P")
LOW  = "#2E86FF"; HIGH  = "#FF7070"     # Tregs-ON  (lighter)
LOWd = "#261CC1"; HIGHd = "#B02020"     # Tregs-OFF (saturated)
## initial conditions of the portraits: Eh0 = 0, 10, ..., 90 (the sweep's probe
## without 99, whose trajectory lies on top of the Eh0 = 90 one)
EH0_ICS = as.numeric(strsplit(Sys.getenv("EH0_ICS", "0,10,20,30,40,50,60,70,80,90"), ",")[[1]])
TMAX    = as.numeric(Sys.getenv("TMAX", "12000"))
sterile = identical(BRANCH, "sterile")

S = as.data.table(readRDS(sprintf("./data/skeleton_%s.rds", BRANCH)))
RS = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", BRANCH)))[, c(CK,"Eh0_off","Eh0_on"), with=FALSE]
S  = merge(S, unique(RS), by = CK, all.x = TRUE)
## effect = "help" keeps only the combinations of that Tregs-effect class, as
## the pathogenic branches already do through their own skeleton files; the
## sterile branch carries help and neutral combinations in one file.
if (!is.null(effect)) {
  stopifnot("treg_effect" %in% names(S))
  S = S[treg_effect == effect]
  if (!nrow(S)) stop("no combinations with treg_effect == ", effect, " in skeleton_", BRANCH, ".rds")
}
BI = S[rescue_mode == "bistable" & is.finite(d_hi) & is.finite(rescue_gap)]
MB = S[rescue_mode == "mono_B"   & is.finite(lo_on) & is.finite(hi_off)]
MA = if (sterile) S[0] else S[rescue_mode == "mono_A" & is.finite(lo_on) & is.finite(hi_off)]
hasBI = nrow(BI) > 0; hasMB = nrow(MB) > 0; hasMA = nrow(MA) > 0
exBI = if (hasBI) BI[which.max(rescue_gap)] else NULL
exMB = if (hasMB) MB[which.max(lo_off - lo_on)] else NULL
exMA = if (hasMA) MA[which.max(lo_on - lo_off)] else NULL
cat(sprintf("[skeleton] %s | bistable %d%s | ON_mono %d%s | mono_A %d\n", BRANCH,
    nrow(BI), if (hasBI) sprintf(" (tau_id %d Kg %g beta_ros %g r_L_P %g, gap %.4f)", exBI$tau_id, exBI$Kg, exBI$beta_ros, exBI$r_L_P, exBI$rescue_gap) else "",
    nrow(MB), if (hasMB) sprintf(" (tau_id %d Kg %g beta_ros %g r_L_P %g, dlow %.4f)", exMB$tau_id, exMB$Kg, exMB$beta_ros, exMB$r_L_P, exMB$lo_off - exMB$lo_on) else "",
    nrow(MA)))
if (!hasBI && !hasMB && !hasMA) stop("no bistable, ON_mono or mono_A combinations in skeleton_", BRANCH, ".rds")
mk_ics = function(cc) EH0_ICS
traj = function(cc, theta, ics) {
  ## tg is keyed on (batch, row_ind) -- see lib/skeleton_defs.R. Looking up on
  ## row_ind alone returned whichever batch happened to be loaded, so a batch-2
  ## combination was integrated with batch 1's parameters and the trajectories
  ## did not converge to the fixed points drawn beside them.
  ## nomatch = NULL, deliberately. On a keyed data.table a MISS still returns one
  ## row, with the QUERY values sitting in the key columns and NA everywhere else
  ## -- so testing is.na(trow$row_ind) never fires and the panel would be drawn
  ## from a row of NA parameters. nomatch = NULL returns zero rows instead.
  trow = tg[.(cc$batch, cc$row_ind), nomatch = NULL]
  if (!nrow(trow))
    stop(sprintf("no tau row for batch %s row_ind %s -- is ./tau_grid_v6_%s.rds present?",
                 cc$batch, cc$row_ind, cc$batch))
  trow = trow[1]
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]]
  p$Kg=cc$Kg; p$beta_ros=cc$beta_ros; p$r_L_P=cc$r_L_P; p$r_LP_P=cc$r_L_P; p$theta_treg=theta
  lapply(ics, function(e0) {
    ## ode_ext = lsodar + the extinction root + the sweep's atol/rtol, trimmed
    ## back to the requested time grid. See lib/extinction_floor.R. This matters
    ## more on the pathogenic branch than the sterile one: without the floor a
    ## cleared P_LP can regrow, and the trajectory drawn would not be the one the
    ## sweep classified.
    o = tryCatch(ode_ext(seed(e0), seq(0, TMAX, by = 2), rhs_app, p, maxsteps = 1e5),
                 error = function(e) NULL)
    if (is.null(o)) NULL else list(Eh = o[,"Eh"], M1 = o[,"M1"], PLP = o[,"PLP"]) })
}
## ---- direction markers on the trajectories ---------------------------------
## Ported from figure_scripts/Fig_2.R (panel A), so the two figures mark direction the
## same way: arrowheads sitting ON the curve, no shafts, placed by arc length in
## inches, N_HEAD per trajectory, skipping any head within HEAD_SEP of one already
## drawn in the panel -- converging trajectories would otherwise stack heads on
## the attractor. `phase` shifts the row along the path so a Tregs-ON and a
## Tregs-OFF trajectory from the same initial condition alternate.
## The sizes are figure_scripts/Fig_2.R's scaled for this page, whose base type is larger.
HEAD_LEN = as.numeric(Sys.getenv("HEAD_LEN", "0.27"))   # head length along the path, inches
HEAD_WID = as.numeric(Sys.getenv("HEAD_WID", "0.15"))   # half-width at its base, inches
N_HEAD   = as.integer(Sys.getenv("N_HEAD", "3"))        # per trajectory, before the spacing filter
HEAD_SEP = as.numeric(Sys.getenv("HEAD_SEP", "0.50"))   # minimum head separation, inches
.head_drawn = new.env(parent = emptyenv())
head_reset = function() assign("p", matrix(numeric(0), nrow = 2), envir = .head_drawn)
head_reset()
add_heads = function(x, y, col, n = N_HEAD, len = HEAD_LEN, wid = HEAD_WID, phase = 0) {
  xi = grconvertX(x, "user", "inches"); yi = grconvertY(y, "user", "inches")
  ok = is.finite(xi) & is.finite(yi); xi = xi[ok]; yi = yi[ok]
  if (length(xi) < 2L) return(invisible())
  d = c(0, cumsum(sqrt(diff(xi)^2 + diff(yi)^2))); L = d[length(d)]
  if (!is.finite(L) || L < 1.5 * len) return(invisible())
  at = function(s) {
    i = max(1L, min(length(d) - 1L, findInterval(s, d)))
    f = if (d[i+1] > d[i]) (s - d[i]) / (d[i+1] - d[i]) else 0
    c(xi[i] + f * (xi[i+1] - xi[i]), yi[i] + f * (yi[i+1] - yi[i]))
  }
  for (s in (seq_len(n) + phase) / (n + 1) * L) {
    pt = at(s)
    u = at(min(s + len/2, L)) - at(max(s - len/2, 0))
    m = sqrt(sum(u^2)); if (!is.finite(m) || m == 0) next
    u = u / m; v = c(-u[2], u[1])
    seen = get("p", envir = .head_drawn)
    if (ncol(seen) && min(sqrt(colSums((seen - pt)^2))) < HEAD_SEP) next
    assign("p", cbind(seen, pt), envir = .head_drawn)
    tri = cbind(pt + u * len/2, pt - u * len/2 + v * wid, pt - u * len/2 - v * wid)
    polygon(grconvertX(tri[1, ], "inches", "user"),
            grconvertY(tri[2, ], "inches", "user"), col = col, border = "white", lwd = 0.6)
  }
}
## ---------------------------------------------------------------------------
## left_main(): panel titles are left-aligned on every panel of every figure, the
## way the main figures set them. plot(main=) centres, and mtext() draws "\n"
## literally, so a two-line title is placed one line at a time.
## ---------------------------------------------------------------------------
left_main = function(ttl, cex = 1.25) {
  ln = strsplit(ttl, "\n", fixed = TRUE)[[1]]
  for (k in seq_along(ln))
    mtext(ln[k], side = 3, line = 0.5 + (length(ln) - k) * 1.3 * cex / 1.25,
          adj = 0, font = 2, cex = par("cex") * cex)
}
add_arrow = function(x, y, col) arrows(x[1], y[1], x[2], y[2], length=.06, angle=22, lwd=2, col=col)
callout = function(tx, ty, ex, ey, lab, col, cex=0.95) {
  ## The leader runs from the label's centre, which the label is then drawn over,
  ## to the attractor itself. It is not shortened: the attractor circles are drawn
  ## again after the callouts, so the end of the leader is covered by its circle
  ## and meets the rim whatever the symbol size.
  segments(tx, ty, ex, ey, col = adjustcolor(col,.8), lwd = 2.5)
  text(tx, ty, lab, col = col, font = 2, cex = cex) }

phase_pathogenic = function(cc, ttl, anchors = NULL, heads = FALSE,
                         ic_marks = FALSE, ic_contour = TRUE, ic_brace = FALSE, ic_drop = NULL, arrows = FALSE,
                         lab_cex = 0.95, ttl_cex = 1.25) {
  ## arrows = FALSE by default: the direction markers are off. Every trajectory
  ## starts at a black cross on the E_h axis and ends on a filled attractor circle,
  ## which carries the direction without adding marks over the curves.
  dir_mark = function(x, y, col, ph)
    if (!arrows) invisible() else if (heads) add_heads(x, y, col, phase = ph) else add_arrow(x, y, col)
  ## PATHOGENIC PORTRAIT: y is P_LP, not M1.
  ## The composite ruler is sqrt((K_E - Eh)^2 + P_LP^2), so (Eh, P_LP) IS the plane
  ## the classification lives in -- the two attractors are separated along it by
  ## construction, and the distance between them on screen is the hi_lo_gap the
  ## CUTOFF gates on. M1 is not in the ruler, which is why the M1 version showed
  ## the lower basin sitting ABOVE the higher one in M1 and read as contradictory.
  ## The sterile figure keeps M1: P_LP is identically 0 there, so an (Eh, P_LP)
  ## plane would collapse onto a line.
  mid_off = (cc$lo_off + cc$hi_off)/2
  mid_on  = if (is.finite(cc$hi_on)) (cc$lo_on + cc$hi_on)/2 else Inf
  ics = mk_ics(cc)
  TON = traj(cc,1,ics); TOFF = traj(cc,0,ics)
  cmp = function(z) sqrt((K_E - tail(z$Eh,1))^2 + tail(z$PLP,1)^2)
  ym = max(c(1, sapply(c(TON,TOFF), function(z) if (is.null(z)) 0 else max(z$PLP)),
             cc$lo_PLP_off, cc$hi_PLP_off, cc$lo_PLP_on, cc$hi_PLP_on), na.rm=TRUE) * 1.14
  plot(NA, xlim=c(-2,102), ylim=c(-0.6,ym), xlab=expression(E[h]), ylab=expression(P[LP]),
       main=""); grid(col="grey92", lty=1)
  left_main(ttl, ttl_cex)
  if (heads) head_reset()
  ONL = OFFL = rep(NA, length(ics))
  for (k in seq_along(ics)) {
    z = TOFF[[k]]; if (!is.null(z)) { OFFL[k] = cmp(z) <= mid_off; co = if (OFFL[k]) LOWd else HIGHd
      lines(z$Eh, z$PLP, col=adjustcolor(co,.9), lwd=1.5, lty=2); dir_mark(z$Eh, z$PLP, co, 0) }
    z = TON[[k]];  if (!is.null(z)) { ONL[k] = cmp(z) <= mid_on; cn = if (ONL[k]) LOW else HIGH
      lines(z$Eh, z$PLP, col=adjustcolor(cn,.9), lwd=1.5, lty=1); dir_mark(z$Eh, z$PLP, cn, 0.5) }
    points(ics[k], PLP0, pch=4, col="black", cex=1.5, lwd=2.5)  # initial condition (black cross), at PLP0
  }
  if (ic_marks) ic_squares(ics, ONL, OFFL, contour = ic_contour, brace = ic_brace,
                           drop = if (is.null(ic_drop)) formals(ic_squares)$drop else ic_drop)
  fp = function(x,y,col) if (is.finite(x) && is.finite(y))
    points(x, y, pch=21, bg=col, col="white", cex=2.0, lwd=1.5)
  fp(cc$hi_Eh_off, cc$hi_PLP_off, HIGHd); fp(cc$lo_Eh_off, cc$lo_PLP_off, LOWd)
  fp(cc$hi_Eh_on,  cc$hi_PLP_on,  HIGH);  fp(cc$lo_Eh_on,  cc$lo_PLP_on,  LOW)
  ## Label anchors, (x, fraction of the y range). The defaults suit the example with
  ## the largest rescue gap; a caller reusing this panel overrides any whose callout
  ## would otherwise run across a trajectory.
  AN = list(hi_off = c(28, .97), hi_on = c(28, .76), lo_on = c(74, .33), lo_off = c(74, .56))
  if (!is.null(anchors)) AN[names(anchors)] = anchors
  ## Label anchors: every trajectory starts on the E_h axis and moves right, so the
  ## upper-left corner is empty. Both high-burden labels are stacked there, clear of
  ## the curves and of each other; the low-burden labels stay mid-height right.
  if (is.finite(cc$hi_Eh_off)) callout(AN$hi_off[1], ym*AN$hi_off[2], cc$hi_Eh_off, cc$hi_PLP_off, "High-burden attractor\n(Tregs-OFF)", HIGHd, cex = lab_cex)
  if (is.finite(cc$hi_Eh_on))  callout(AN$hi_on[1],  ym*AN$hi_on[2], cc$hi_Eh_on,  cc$hi_PLP_on,  "High-burden attractor\n(Tregs-ON)",  HIGH, cex = lab_cex)
  if (is.finite(cc$lo_Eh_on))  callout(AN$lo_on[1],  ym*AN$lo_on[2], cc$lo_Eh_on,  cc$lo_PLP_on,  "Low-burden attractor\n(Tregs-ON)",   LOW, cex = lab_cex)
  if (is.finite(cc$lo_Eh_off)) callout(AN$lo_off[1], ym*AN$lo_off[2], cc$lo_Eh_off, cc$lo_PLP_off, "Low-burden attractor\n(Tregs-OFF)",  LOWd, cex = lab_cex)
  ## the circles again, over the leader ends
  fp(cc$hi_Eh_off, cc$hi_PLP_off, HIGHd); fp(cc$lo_Eh_off, cc$lo_PLP_off, LOWd)
  fp(cc$hi_Eh_on,  cc$hi_PLP_on,  HIGH);  fp(cc$lo_Eh_on,  cc$lo_PLP_on,  LOW)
  legend("topright", bty="n", cex=1.1, seg.len=2.4, legend=c("Tregs-ON","Tregs-OFF"),
         col="#555555", lwd=2.2, lty=c(1,2))
}
phase_sterile = function(cc, ttl, anchors = NULL, heads = FALSE,
                         ic_marks = FALSE, ic_contour = TRUE, ic_brace = FALSE, ic_drop = NULL, arrows = FALSE,
                         lab_cex = 0.95, ttl_cex = 1.25) {
  ## arrows = FALSE by default: the direction markers are off. Every trajectory
  ## starts at a black cross on the E_h axis and ends on a filled attractor circle,
  ## which carries the direction without adding marks over the curves.
  dir_mark = function(x, y, col, ph)
    if (!arrows) invisible() else if (heads) add_heads(x, y, col, phase = ph) else add_arrow(x, y, col)
  ## each condition is judged against ITS OWN two attractors. Using OFF's midpoint for the
  ## ON trajectories mislabels an ON higher-basin endpoint as "lower basin" whenever ON's
  ## higher basin has moved below OFF's midpoint -- common on the pathogenic branch.
  mid_off = (cc$lo_off + cc$hi_off)/2
  mid_on  = if (is.finite(cc$hi_on)) (cc$lo_on + cc$hi_on)/2 else Inf
  ics = mk_ics(cc)
  TON = traj(cc,1,ics); TOFF = traj(cc,0,ics)
  cmp = function(z) sqrt((K_E - tail(z$Eh,1))^2 + tail(z$PLP,1)^2)
  ym = max(sapply(c(TON,TOFF), function(z) if (is.null(z)) 0 else max(z$M1))) * 1.14
  plot(NA, xlim=c(-2,102), ylim=c(-0.6,ym), xlab=expression(E[h]), ylab=expression(M[1]),
       main=""); grid(col="grey92", lty=1)
  left_main(ttl, ttl_cex)
  if (heads) head_reset()
  ONL = OFFL = rep(NA, length(ics))
  for (k in seq_along(ics)) {
    z = TOFF[[k]]; if (!is.null(z)) { OFFL[k] = cmp(z) <= mid_off; co = if (OFFL[k]) LOWd else HIGHd
      lines(z$Eh, z$M1, col=adjustcolor(co,.9), lwd=1.5, lty=2); dir_mark(z$Eh, z$M1, co, 0) }
    z = TON[[k]];  if (!is.null(z)) { ONL[k] = cmp(z) <= mid_on; cn = if (ONL[k]) LOW else HIGH
      lines(z$Eh, z$M1, col=adjustcolor(cn,.9), lwd=1.5, lty=1); dir_mark(z$Eh, z$M1, cn, 0.5) }
    points(ics[k], 0, pch=4, col="black", cex=1.5, lwd=2.5)   # initial condition (black cross)
  }
  if (ic_marks) ic_squares(ics, ONL, OFFL, contour = ic_contour, brace = ic_brace,
                           drop = if (is.null(ic_drop)) formals(ic_squares)$drop else ic_drop)
  fp = function(x,y,col) if (is.finite(x) && is.finite(y))
    points(x, y, pch=21, bg=col, col="white", cex=2.0, lwd=1.5)
  fp(cc$hi_Eh_off, cc$hi_M1_off, HIGHd); fp(cc$lo_Eh_off, cc$lo_M1_off, LOWd)
  fp(cc$hi_Eh_on,  cc$hi_M1_on,  HIGH);  fp(cc$lo_Eh_on,  cc$lo_M1_on,  LOW)
  ## Label anchors, (x, fraction of the y range). The defaults suit the example with
  ## the largest rescue gap; a caller reusing this panel overrides any whose callout
  ## would otherwise run across a trajectory.
  AN = list(hi_off = c(28, .97), hi_on = c(28, .76), lo_on = c(74, .33), lo_off = c(74, .56))
  if (!is.null(anchors)) AN[names(anchors)] = anchors
  ## Label anchors. Every trajectory starts on the E_h axis and moves right, so the
  ## upper-left corner of the panel is empty in both the bistable and the monostable
  ## case: both high-burden labels go there, stacked, well clear of the curves and of
  ## each other. The two low-burden labels stay at mid-height on the right,
  ## where the only thing above them is the legend.
  if (is.finite(cc$hi_Eh_off)) callout(AN$hi_off[1], ym*AN$hi_off[2], cc$hi_Eh_off, cc$hi_M1_off, "High-burden attractor\n(Tregs-OFF)", HIGHd, cex = lab_cex)
  if (is.finite(cc$hi_Eh_on))  callout(AN$hi_on[1],  ym*AN$hi_on[2], cc$hi_Eh_on,  cc$hi_M1_on,  "High-burden attractor\n(Tregs-ON)",  HIGH, cex = lab_cex)
  if (is.finite(cc$lo_Eh_on))  callout(AN$lo_on[1],  ym*AN$lo_on[2], cc$lo_Eh_on,  cc$lo_M1_on,  "Low-burden attractor\n(Tregs-ON)",   LOW, cex = lab_cex)
  if (is.finite(cc$lo_Eh_off)) callout(AN$lo_off[1], ym*AN$lo_off[2], cc$lo_Eh_off, cc$lo_M1_off, "Low-burden attractor\n(Tregs-OFF)",  LOWd, cex = lab_cex)
  ## the circles again, over the leader ends
  fp(cc$hi_Eh_off, cc$hi_M1_off, HIGHd); fp(cc$lo_Eh_off, cc$lo_M1_off, LOWd)
  fp(cc$hi_Eh_on,  cc$hi_M1_on,  HIGH);  fp(cc$lo_Eh_on,  cc$lo_M1_on,  LOW)
  legend("topright", bty="n", cex=1.1, seg.len=2.4, legend=c("Tregs-ON","Tregs-OFF"),
         col="#555555", lwd=2.2, lty=c(1,2))
}
## ---------------------------------------------------------------------------
## basinstrip(): which basin each initial condition reaches, one row per arm --
## the layout of figure_scripts/Fig_2.R panel B, with Tregs in place of efferocytosis.
## A third row underneath marks, in black, the initial conditions that reach the
## low-burden basin WITH Tregs but not without: the rescue gap, as the set of
## injuries it actually covers rather than as a number. Basins are read off the
## integrated trajectories, not off the stored thresholds, so the row and the
## portrait above it are the same computation.
## ---------------------------------------------------------------------------
## xlim defaults to the portraits' own E_h range, so a strip placed under one
## lines up with it column for column and can be read downwards.
## WCAG relative luminance -> the higher-contrast ink for text drawn on a fill.
ink = function(fill) {
  v = col2rgb(fill)[, 1] / 255
  v = ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
  L = sum(c(0.2126, 0.7152, 0.0722) * v)
  if ((L + 0.05) / 0.05 > 1.05 / (L + 0.05)) "black" else "white"
}
## ---------------------------------------------------------------------------
## ic_squares(): under each initial-condition cross, a small square split on its
## bottom-left/top-right diagonal. The lower-right half takes the colour of the
## attractor reached with Tregs ON, the upper-left half the one reached with Tregs
## OFF, so a square of two different colours IS an initial condition whose outcome
## Tregs changed. A black outline marks the ones inside the rescue gap (lower basin
## with Tregs, higher without). Drawn in the bottom margin, between the tick labels
## and the axis title, so it reads as part of the x axis rather than a second panel.
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## brace_under(): a horizontal curly brace spanning [x0, x1] in user x, hanging
## from y_top. Built in INCHES out of four quarter circles joined by two straight
## shoulders and only then mapped back to user coordinates, so the curls stay
## round whatever the two axis scales are. xpd = NA: it is drawn in the margin,
## under the squares of ic_squares(). Returns its midpoint in user x.
## ---------------------------------------------------------------------------
brace_under = function(x0, x1, y_top, depth_in = 0.14, lwd = 4, col = "black", n = 24) {
  a = grconvertX(x0, "user", "inches"); b = grconvertX(x1, "user", "inches")
  y = grconvertY(y_top, "user", "inches")
  r = min(depth_in / 2, (b - a) / 4); m = (a + b) / 2
  arc = function(cx, cy, t0, t1) { t = seq(t0, t1, length.out = n); cbind(cx + r * cos(t), cy + r * sin(t)) }
  P = rbind(arc(a + r, y,           pi,         3 * pi / 2),
            cbind(c(a + r, m - r),  c(y - r, y - r)),
            arc(m - r, y - 2 * r,   pi / 2,     0),
            arc(m + r, y - 2 * r,   pi,         pi / 2),
            cbind(c(m + r, b - r),  c(y - r, y - r)),
            arc(b - r, y,           3 * pi / 2, 2 * pi))
  op = par(xpd = NA); on.exit(par(op))
  lines(grconvertX(P[, 1], "inches", "user"), grconvertY(P[, 2], "inches", "user"),
        lwd = lwd, col = col, lend = 1, ljoin = 1)
  invisible(grconvertX(m, "inches", "user"))
}
ic_squares = function(x, on_low, off_low, contour = TRUE, half_in = 0.190, drop = 1.76,
                      contour_lwd = 5.0, gap_in = 0.055, tile_gap_lwd = 3,
                      brace = FALSE, brace_lab = "rescue gap", brace_cex = 1.15,
                      brace_depth = 0.26, brace_pad = 0.025, brace_lab_gap = 0.03, brace_lwd = 4) {
  u = par("usr"); p = par("pin")
  dx = half_in / p[1] * (u[2] - u[1])
  dy = half_in / p[2] * (u[4] - u[3])
  ## the outline stands off the fill by gap_in, leaving a white ring: drawn on the
  ## edge it read as part of the colour and muddied the two halves.
  gx = (half_in + gap_in) / p[1] * (u[2] - u[1])
  gy = (half_in + gap_in) / p[2] * (u[4] - u[3])
  y  = u[3] - drop * par("cxy")[2]
  op = par(xpd = NA); on.exit(par(op))
  for (k in seq_along(x)) {
    if (is.na(on_low[k]) || is.na(off_low[k])) next
    cON  = if (on_low[k])  LOW  else HIGH
    cOFF = if (off_low[k]) LOWd else HIGHd
    polygon(c(x[k]-dx, x[k]+dx, x[k]+dx), c(y-dy, y-dy, y+dy), col = cON,  border = NA)
    polygon(c(x[k]-dx, x[k]-dx, x[k]+dx), c(y-dy, y+dy, y+dy), col = cOFF, border = NA)
    ## a white rule along the diagonal, so the two fills never touch
    segments(x[k]-dx, y-dy, x[k]+dx, y+dy, col = "white", lwd = tile_gap_lwd)
    if (contour && on_low[k] && !off_low[k])
      rect(x[k]-gx, y-gy, x[k]+gx, y+gy, border = "black", lwd = contour_lwd)
  }
  ## THE RESCUE GAP, named under the squares it covers. The contoured squares are
  ## exactly the initial conditions from which only the Tregs-ON system reaches the
  ## low-burden attractor, so the brace spans the outer edges of the first and last
  ## of them. The word hangs off the TONGUE of the brace, brace_lab_gap inches under
  ## it, rather than off a margin line: it then stays tight to the brace whatever
  ## mgp and drop do with the axis label and the squares above it.
  if (brace && contour) {
    sel = which(!is.na(on_low) & !is.na(off_low) & on_low & !off_low)
    if (length(sel) >= 2L) {
      bot_in = grconvertY(y - gy, "user", "inches") - brace_pad      # just under the contour rects
      mid = brace_under(min(x[sel]) - gx, max(x[sel]) + gx,
                        grconvertY(bot_in, "inches", "user"),
                        depth_in = brace_depth, lwd = brace_lwd)
      ## adj = c(0.5, 1) sets the TOP of the word, so the distance to the tongue is
      ## brace_lab_gap whatever the type size.
      text(mid, grconvertY(bot_in - brace_depth - brace_lab_gap, "inches", "user"),
           brace_lab, adj = c(0.5, 1), cex = brace_cex, xpd = NA)
    }
  }
}
set.seed(1)
## violin + inner box + explicit outlier points. gstep = fixed grid spacing (NULL = dense auto).
## ---------------------------------------------------------------------------
## gaphist(): the rescue gap for ONE rescue_mode, as discrete bars.
##
## Replaces the stability scatters that used to occupy this column. 06c's own
## artefact test returned rho = -0.40 / -0.35 between percentage change and
## baseline |Re(lambda)|, i.e. the largest percentage stiffening is the smallest
## baseline -- so the direction survives (a sign test is unaffected) but the
## magnitude is not quotable and did not earn four panels. The gap belongs here
## instead: it is the quantity the whole rescue argument rests on, and its shape
## differs sharply BETWEEN the rows this figure already splits by.
##
## Drawn as bars on the EH0_SCAN lattice, never a density -- the values are
## multiples of the bisection step and a smooth curve would invent structure between
## them. "0" means below the detection limit, not "no effect".
##
## The censored class (Tregs ON reaches the low basin from every probed Eh0, so
## Eh0_on* does not exist and the gap is a lower bound) gets its OWN bar at the
## right, separated by a
## break, because it is the largest effect and averaging it in would understate
## the rescue.
## ---------------------------------------------------------------------------
## THE SCALE CHANGED IN v6, AND SO DID THE REASON FOR IT.
## v5 read the threshold off the probe spacing: Eh0* was the midpoint of the two
## bracketing values, so every gap landed on a multiple of 2.5 and the detection
## floor was 5. The bars were 5 wide because nothing finer existed, and 98% of the
## sterile mass sat in the single bar at zero.
## v6 BISECTS the threshold, so the floor is now EH0_TOL, and the smallest gap the
## run can express is one bisection step:
##       STEP = PROBE_BY / 2^ceiling(log2(PROBE_BY / EH0_TOL))
## which at the defaults (10, 0.1) is 0.078125. Measured on the first v6 batch:
## 138 sterile gaps take 19 distinct values, every one a multiple of that step,
## with 36% at exactly 0 and 28% at one step. So the distribution is real but its
## bottom is still pressed against a floor -- just a floor 64x finer than v5's.
## FLOOR is what the legend calls "at floor" and it is the same quantity as STEP:
## a gap below one bisection step is not a small measurement, it is an unresolved
## one. Do not set them independently.
##
## Values are BINNED by rounding onto the display lattice, never matched exactly.
## v5 tabulated with factor(v, levels = seq(..., by = ray)), which turned every
## off-lattice value into NA and dropped it -- the bars summed to well under 1 with
## no warning and half the mono_B mass disappeared. The hard check below refuses to
## draw a panel whose bars do not account for every observation.
STEP_DEFAULT = local({
  pb = as.numeric(Sys.getenv("PROBE_BY", "10")); tol = as.numeric(Sys.getenv("EH0_TOL", "0.1"))
  pb / 2^ceiling(log2(pb / tol)) })
## ---------------------------------------------------------------------------
## LOG Y, ONE BIN WIDTH. The bars decay so steeply (on every branch the second
## lattice bar is already below a quarter of the peak, and ~95% of the mass lies
## below 1) that on a linear y axis the panel read as a spike at 0 and a flat line
## to the maximum, whatever was done to x. An earlier version cut the x axis with a
## // and rebinned the tail about 1 unit wide, which made the tail no more visible
## and left two bin widths in one panel, so a bar's height meant different things on
## either side of the break. The fraction is now on a log y axis and the bisection
## lattice runs unbroken from 0 to the maximum: every bar is the same width and the
## four decades of the decay are all legible. Minor lattice ticks are not drawn --
## there are hundreds of them at this step and they merge into a black band.
## The in-panel text reports min, mean and max of the gap, nothing else.
## ---------------------------------------------------------------------------
gaphist = function(X, col, ttl, step = STEP_DEFAULT, floor_at = STEP_DEFAULT) {
  ## isTRUE() takes a length-1 logical -- on a COLUMN it returns a single FALSE, so
  ## !isTRUE(gap_censored) was always TRUE and the censored combos were included in
  ## the bars AND counted again in the "all" bar, doubling the denominator. Compare
  ## elementwise.
  cens = if ("gap_censored" %in% names(X)) X$gap_censored %in% TRUE else rep(FALSE, nrow(X))
  v  = X[is.finite(rescue_gap) & !cens, rescue_gap]
  nc = sum(cens)
  if (!length(v)) { blank("no gap data"); return(invisible()) }
  tot = length(v) + nc
  ## HEAD: lattice bars, binned by rounding (never matched exactly, see above)
  rng = range(v)
  br  = seq(min(0, floor(rng[1]/step)*step), max(step, ceiling(rng[2]/step)*step), by = step)
  vb  = round(v / step) * step
  tb  = table(factor(vb, levels = br))
  if (sum(tb) != length(v))
    stop(sprintf("gaphist: %d of %d values fell outside the lattice -- the bars would understate the data",
                 length(v) - sum(tb), length(v)))
  ys  = as.numeric(tb) / tot
  xl0 = br[1] - step/2; xr = br[length(br)] + step/2
  ## one lattice bin width across the whole range, so every bar means the same
  ## thing; the fraction carries the dynamic range instead (see above).
  ytop = max(ys); ybot = 10^floor(log10(min(ys[ys > 0])))
  plot(NA, xlim = c(xl0, xr), ylim = c(ybot, ytop * 2.2), xaxs = "i", log = "y",
       xlab = expression(Delta * E[h] * "(0)*   (OFF - ON)"), ylab = "Fraction of combinations",
       main = "", yaxt = "n")
  left_main(ttl)
  ## powers of ten: 0.00001 and the like ran into the axis title
  dec = 10^(floor(log10(ybot)):0)
  axis(2, at = dec, labels = parse(text = sprintf("10^%d", log10(dec))), las = 1)
  abline(h = dec, col = "grey92")
  abline(v = 0, lty = 3, col = "grey45")
  ## no minor lattice ticks: at this step there are hundreds of them and they merge
  ## into a solid black band under the axis.
  ## The bars carried their percentage in text; the labels crowded the head of the
  ## distribution and the y axis already gives the fraction, so they are not drawn.
  rect(br - step/2, ybot, br + step/2, pmax(ys, ybot), col = col, border = NA)
  box()
  legend("topright", bty = "n", cex = 1.05,
         legend = c(sprintf("min  %.3g", min(v)), sprintf("mean %.3g", mean(v)), sprintf("max  %.3g", max(v)),
                    if (nc) sprintf("ON reaches low basin throughout: %.0f%%", 100*nc/tot)))
  cat(sprintf("[skeleton gap] %-34s n=%6d | min %.3g mean %.3g max %.3g | unresolved %3.0f%% | %d bins of %.3g | y from %.3g | censored %5.1f%%\n",
              ttl, tot, min(v), mean(v), max(v), 100*mean(abs(v) < floor_at),
              length(br), step, ybot, 100*nc/tot))
}
## ---------------------------------------------------------------------------
## gapdens(): the same lattice distribution as gaphist(), in a SHORT panel, for a
## page that carries two branches. Identical binning and log y axis -- a kernel
## density is not used here, because about a third of the gaps are exactly 0 and
## smoothing turns that atom into a hump spread over neighbouring values that no
## combination takes. The three numbers that the long panel puts in a legend block
## are printed on one line instead.
## ---------------------------------------------------------------------------
gapdens = function(X, col, ttl, step = STEP_DEFAULT, floor_at = STEP_DEFAULT,
                   xlab = expression(Delta * E[h] * "(0)*   (OFF - ON)"),
                   ttl_cex = 1.25, xlab_cex = par("cex.lab"), xaxis_cex = par("cex.axis")) {
  cens = if ("gap_censored" %in% names(X)) X$gap_censored %in% TRUE else rep(FALSE, nrow(X))
  v = X[is.finite(rescue_gap) & !cens, rescue_gap]
  if (!length(v)) { blank("no gap data"); return(invisible()) }
  tot = length(v) + sum(cens)
  rng = range(v)
  br = seq(min(0, floor(rng[1]/step)*step), max(step, ceiling(rng[2]/step)*step), by = step)
  tb = table(factor(round(v/step) * step, levels = br))
  if (sum(tb) != length(v))
    stop(sprintf("gapdens: %d of %d values fell outside the lattice", length(v) - sum(tb), length(v)))
  ys = as.numeric(tb) / tot
  ytop = max(ys); ybot = 10^floor(log10(min(ys[ys > 0])))
  plot(NA, xlim = c(br[1] - step/2, br[length(br)] + step/2), ylim = c(ybot, ytop * 4),
       xaxs = "i", log = "y", xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "fraction")
  axis(1, cex.axis = xaxis_cex)
  title(xlab = xlab, cex.lab = xlab_cex)
  left_main(ttl, ttl_cex)
  dec = 10^(floor(log10(ybot)):0)
  axis(2, at = dec, labels = parse(text = sprintf("10^%d", log10(dec))), las = 1)
  abline(h = dec, col = "grey92"); abline(v = 0, lty = 3, col = "grey45")
  rect(br - step/2, ybot, br + step/2, pmax(ys, ybot), col = col, border = NA)
  box()
  legend("topright", bty = "n", cex = 0.95,
         legend = sprintf("%.0f%% at 0    %.0f%% > 0    max %.1f",
                          100*mean(abs(v) < floor_at), 100*mean(v >= floor_at), max(v)))
  cat(sprintf("[gapdens] %-28s n=%6d | at 0 %3.0f%% | >0 %3.0f%% | max %.3g\n",
              ttl, tot, 100*mean(abs(v) < floor_at), 100*mean(v >= floor_at), max(v)))
}
vio = function(L, cols, ttl, yl, gstep = NULL, labs = names(L), ttl_cex = 1.25, lab_cex = 1.1) {
  L  = lapply(L, function(x) x[is.finite(x)])
  av = unlist(L, use.names = FALSE); n = length(L)
  yr = range(av); pad = diff(yr) * .07; ylim = c(yr[1] - pad, yr[2] + pad)
  plot(NA, xlim = c(.4, n + .6), ylim = ylim, xaxt = "n", xlab = "", ylab = yl, main = "")
  left_main(ttl, ttl_cex)
  ticks = if (is.null(gstep)) pretty(ylim, n = 12) else
          seq(floor(ylim[1]/gstep)*gstep, ceiling(ylim[2]/gstep)*gstep, by = gstep)
  abline(h = ticks, col = "grey92"); abline(h = 0, lty = 2, col = "grey40")
  ## tick labels drawn by hand at 45 degrees: axis() silently drops labels that
  ## would overlap, and at the larger font the six pathogenic labels did.
  axis(1, at = seq_len(n), labels = FALSE)
  text(seq_len(n), par("usr")[3] - 0.035 * diff(par("usr")[3:4]), labels = labs,
       srt = 60, adj = 1, xpd = TRUE, cex = lab_cex)
  axis(2, at = ticks, labels = FALSE, tcl = -.25)
  for (j in seq_len(n)) {
    x = L[[j]]; if (!length(x)) next
    if (length(x) > 2 && diff(range(x)) > 0) {
      d = density(x, from = min(x), to = max(x)); w = .38 * d$y / max(d$y)
      polygon(c(j - w, rev(j + w)), c(d$x, rev(d$x)),
              col = adjustcolor(cols[j], .45), border = cols[j], lwd = 1.3)
    }
    q = as.numeric(quantile(x, c(.25, .5, .75))); iqr = q[3] - q[1]
    lo = max(min(x), q[1] - 1.5*iqr); hi = min(max(x), q[3] + 1.5*iqr)
    segments(j, lo, j, hi, col = "grey20", lwd = 1.1)
    rect(j - .055, q[1], j + .055, q[3], col = "white", border = "grey20", lwd = 1.1)
    points(j, q[2], pch = 19, cex = 1.1, col = "grey10")
    out = x[x < q[1] - 1.5*iqr | x > q[3] + 1.5*iqr]
    if (length(out)) points(j + runif(length(out), -.07, .07), out,
                            pch = 19, cex = .5, col = adjustcolor(cols[j], .9))
  }
}
## stability change as a PERCENTAGE of the OFF return rate.
## positive = |Re(lambda)| larger under ON = faster return = stiffer well.
pct = function(on, off) { z = 100 * (abs(on) - abs(off)) / abs(off); z[is.finite(z)] }
## a class can be empty on a branch (pathogenic may have no mono_B at all), so the
## corresponding row is drawn as a labelled blank rather than erroring on an empty table.
blank = function(msg) { plot.new(); box(col="grey85")
  text(0.5, 0.5, msg, cex = 1.0, col = "grey40") }
phase = if (sterile) phase_sterile else phase_pathogenic
BLOCKS = list(1:3, 4:6)                    # bistable, ON_mono
if (hasMA) BLOCKS[[3]] = 7:9               # mono_A, only when the class is non-empty
NBLK = length(BLOCKS)
.top = unlist(lapply(BLOCKS, function(b) rep(b[1], 2)))   # phase portrait, 2 cells wide
.bot = unlist(lapply(BLOCKS, function(b) b[2:3]))         # the two summary panels
## panels = a function: the CALLER owns the device and the layout, and draws its
## own arrangement from the helpers below. It can be used to put panels
## from two branches on one page. Default (NULL) keeps the one-branch page.
if (is.null(panels)) {
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  png(out_png, width = 3400 * NBLK, height = 5200, res = 250, pointsize = 38)
  layout(matrix(c(.top, .top, .bot), nrow = 3, byrow = TRUE),
         widths  = rep(1, 2 * NBLK), heights = c(1, 1, 1.85))
}
MBIG = c(5.2, 5.6, 4.0, 1.6); MSML = c(8.6, 5.4, 4.4, 1.2); MGAP = c(8.6, 5.4, 3.2, 1.2)
par(cex.main=1.25, cex.lab=1.3, cex.axis=1.15, mgp=c(3.4,0.9,0))
## attractor shifts: composite, Eh and (pathogenic only) P_LP, low (L) and high (H)
shifts = function(X, bistable) {
  L = if (bistable) list(`comp L`=X$d_lo, `comp H`=X$d_hi, `dEh L`=X$d_lo_Eh, `dEh H`=X$d_hi_Eh,
                         `dPLP L`=X$d_lo_PLP, `dPLP H`=X$d_hi_PLP)
      else list(`comp L`=X$lo_on-X$lo_off, `comp H`=X$lo_on-X$hi_off,
                `dEh L`=X$lo_Eh_on-X$lo_Eh_off, `dEh H`=X$lo_Eh_on-X$hi_Eh_off,
                `dPLP L`=X$lo_PLP_on-X$lo_PLP_off, `dPLP H`=X$lo_PLP_on-X$hi_PLP_off)
  if (sterile) L[1:4] else L }
vcols = rep(c(LOW, HIGH), if (sterile) 2 else 3)
## tick labels name the subtraction: bistable compares each basin with itself
## (L^ON - L^OFF, H^ON - H^OFF); ON_mono has one basin, compared with both of OFF's
labs_bi = expression(b ~ (L^ON - L^OFF), b ~ (H^ON - H^OFF),
                     E[h] ~ (L^ON - L^OFF), E[h] ~ (H^ON - H^OFF),
                     P[LP] ~ (L^ON - L^OFF), P[LP] ~ (H^ON - H^OFF))
labs_mb = expression(b ~ (L^ON - L^OFF), b ~ (L^ON - H^OFF),
                     E[h] ~ (L^ON - L^OFF), E[h] ~ (L^ON - H^OFF),
                     P[LP] ~ (L^ON - L^OFF), P[LP] ~ (L^ON - H^OFF))
if (sterile) { labs_bi = labs_bi[1:4]; labs_mb = labs_mb[1:4] }
if (!is.null(panels)) { panels(environment()); return(invisible(environment())) }
## ---- block 1: Tregs ON bistable (A, B, C) ----
par(mar=MBIG); if (hasBI) phase(exBI, "A   Tregs-ON bistable") else blank("no bistable combos")
par(mar=MSML); if (hasBI) vio(shifts(BI, TRUE), vcols, sprintf("B   Fixed point shift\n(n = %d)", nrow(BI)),
                              "ON - OFF", gstep = 10, labs = labs_bi) else blank("")
par(mar=MSML); if (hasBI) gaphist(BI, LOW, "C   Rescue gap\n(bistable)") else blank("")
## ---- block 2: Tregs ON monostable, lower basin (D, E, F) ----
par(mar=MBIG); if (hasMB) phase(exMB, "D   Tregs-ON monostable") else blank(paste0("no ", MONO_LAB, " combos"))
par(mar=MSML); if (hasMB) vio(shifts(MB, FALSE), vcols, sprintf("E   Fixed point shift\n(n = %d)", nrow(MB)),
                              "ON - OFF", gstep = 10, labs = labs_mb) else blank("")
## F is intentionally empty: the ON_mono class has no resolved rescue gap (every
## combination is censored), so the bar panel carried no information. The cell
## stays blank to keep the two blocks aligned.
par(mar=MSML); plot.new()
## ---- block 3: Tregs ON monostable, higher basin (G, H, I), only if non-empty ----
if (hasMA) {
  par(mar=MBIG); phase(exMA, "G   Tregs-ON monostable")
  par(mar=MSML); vio(shifts(MA, FALSE), vcols, sprintf("H   Fixed point shift\n(n = %d)", nrow(MA)), "ON - OFF", gstep = 10, labs = labs_mb)
  par(mar=MSML); gaphist(MA, HIGHd, "I   Rescue gap - mono_A")
}
dev.off()
cat(sprintf("[skeleton] wrote %s\n", out_png))
invisible(list(BI = BI, MB = MB, MA = MA, exBI = exBI, exMB = exMB, exMA = exMA))
}

## ===========================================================================
## Page builders. Both drive fig_skeleton() through its `panels` callback, so a
## page can mix branches and every panel is one the function already draws.
##   fig_portraits() -- the two phase portraits of a branch, each initial condition
##     carrying a split square under the axis (see ic_squares).
##   fig_summaries() -- per branch, the two fixed-point-shift panels and the
##     rescue-gap distribution of the bistable class.
## layout() picks cex from the grid shape (0.83 at two rows/columns, 0.66 above),
## so it is pinned here and every page in the set renders text at the same size.
## ===========================================================================
## shift and rescue-gap panels of Extended Data Fig. 1 and Supplementary Fig. 2
## (fig_summaries(), figure_scripts/Supplementary_Fig_2.R): panel titles, x-axis text, margins and the gap axis title
SUM_TTL = 1.6; SUM_XLAB = 1.5; SUM_XAX = 1.45
SUM_MAR = c(10.4, 5.4, 5.0, 1.2)
SUM_GAP_XLAB = expression(Delta * E[h] * "(0)*" == E[h] * "(0)*"[OFF] - E[h] * "(0)*"[ON])
fig_portraits = function(branch, out_png, effect = NULL, lab = "", anchors = NULL,
                         width = 6400, height = 2900) {
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  png(out_png, width = width, height = height, res = 250, pointsize = 38)
  layout(matrix(1:2, nrow = 1), widths = c(1, 1))
  par(cex = 0.66)
  fig_skeleton(branch, NULL, effect = effect, panels = function(e) {
    MSQ = c(6.2, e$MBIG[2], e$MBIG[3], e$MBIG[4])   # room under the axis for the squares
    par(mar = MSQ, mgp = c(4.0, 0.9, 0))
    e$phase(e$exBI, sprintf("A   %s, Tregs-ON bistable", lab),
            anchors = anchors, ic_marks = TRUE, ic_contour = TRUE)
    par(mar = MSQ, mgp = c(4.0, 0.9, 0))
    e$phase(e$exMB, sprintf("B   %s, Tregs-ON monostable", lab),
            ic_marks = TRUE, ic_contour = FALSE)
  })
  dev.off()
  cat(sprintf("[fig] wrote %s\n", out_png))
}

fig_summaries = function(specs, out_png, width = 9600, height = 6000) {
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  png(out_png, width = width, height = height, res = 250, pointsize = 38)
  layout(matrix(seq_len(3 * length(specs)), nrow = length(specs), byrow = TRUE))
  par(cex = 0.66)
  for (sp in specs) {
    fig_skeleton(sp$branch, NULL, effect = sp$effect, panels = function(e) {
      par(mar = SUM_MAR)
      e$vio(e$shifts(e$BI, TRUE), e$vcols,
            sprintf("%s   %s, Tregs-ON bistable\n(n = %d)", sp$tags[1], sp$lab, nrow(e$BI)),
            "ON - OFF", gstep = 10, labs = e$labs_bi, ttl_cex = SUM_TTL, lab_cex = SUM_XLAB)
      par(mar = SUM_MAR)
      e$vio(e$shifts(e$MB, FALSE), e$vcols,
            sprintf("%s   %s, Tregs-ON monostable\n(n = %d)", sp$tags[2], sp$lab, nrow(e$MB)),
            "ON - OFF", gstep = 10, labs = e$labs_mb, ttl_cex = SUM_TTL, lab_cex = SUM_XLAB)
      par(mar = SUM_MAR)
      e$gapdens(e$BI, e$LOW, sprintf("%s   %s, rescue gap\n(bistable)", sp$tags[3], sp$lab),
                xlab = SUM_GAP_XLAB, ttl_cex = SUM_TTL, xlab_cex = SUM_XLAB, xaxis_cex = SUM_XAX)
    })
  }
  dev.off()
  cat(sprintf("[fig] wrote %s\n", out_png))
}

## ===========================================================================
## fig_portrait_row() -- the A, B phase-portrait half of the merged v2 figures:
## figure_scripts/Fig_3.R (sterile, Fig. 3) and figure_scripts/Extended_Data_Fig_2.R
## (pathogenic, Extended Data Fig. 2).
## It is one row of the former 2 x 2 page (sterile over pathogenic) and its device,
## defined ONCE so the two figures cannot drift apart: the same MSQ, MGP_LAB and
## IC_DROP, the same 6400 x 2900 at res 250 with cex pinned at 0.66, and the
## same hand-tuned callout anchors, which are a fixed property of each branch's
## portrait and so live here rather than at the call sites.
## use: fig_portrait_row(out_png, "sterile")  /  ..., "pathogenic_help")
## ===========================================================================
FIG_PORTRAIT_SPEC = list(
  ## the sterile call: the lower-basin Tregs-OFF label sat to the
  ## RIGHT of its attractor and the higher-basin Tregs-ON callout ran across both
  ## red trajectories
  sterile = list(
    effect = "help", lab = "Sterile", lab_cex = 0.95,
    bi = list(lo_off = c(30, .50), hi_on = c(48, .76)), mb = NULL),
  ## the pathogenic call: every trajectory climbs in P_LP between E_h = 20 and 45
  ## and then sweeps right and down, so the far left above the crosses and the
  ## upper right below the legend are the empty corners
  pathogenic_help = list(
    effect = NULL, lab = "Pathogenic", lab_cex = 0.85,
    bi = list(hi_off = c(17, .99), hi_on = c(52, .84),
              lo_off = c(82, .66), lo_on = c(86, .44)),
    mb = list(hi_off = c(17, .99), lo_off = c(54, .99), lo_on = c(84, .74))),
  ## Supplementary Fig. 2, rendered at res 375 on its 9600-wide page (see figure_scripts/Supplementary_Fig_2.R)
  ## every trajectory peaks near P_LP = 90 and falls away to the right, so the strip
  ## above the peaks and the corner under the legend are the empty places
  pathogenic_harm = list(
    effect = NULL, lab = "Tregs harm", lab_cex = 0.85,
    bi = list(hi_off = c(12, .97), hi_on = c(49, .97),
              lo_off = c(85, .74), lo_on = c(88, .60)),
    mb = list(hi_off = c(17, .97), lo_off = c(86, .76), lo_on = c(87, .63))))

fig_portrait_row = function(out_png, branch, tags = c("A", "B"),
                            width = 6400, height = 2900, res = 250) {
  sp = FIG_PORTRAIT_SPEC[[branch]]
  if (is.null(sp)) stop("[fig_portrait_row] no spec for branch ", branch)
  dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
  png(out_png, width = width, height = height, res = res, pointsize = 38)
  layout(matrix(1:2, nrow = 1), widths = c(1, 1))
  par(cex = 0.66)   # pinned, so text matches the former 2 x 2 page
  fig_skeleton(branch, NULL, effect = sp$effect, panels = function(e) {
    ## The strip under the axis runs: tick labels, the E_h label (MGP_LAB), the
    ## squares (IC_DROP), the brace and the words "rescue gap". MGP_LAB lifts the
    ## label above the squares and IC_DROP pushes the squares clear of it; MSQ[1] is
    ## then the height that whole stack needs, with no empty band left under it.
    MSQ = c(6.6, e$MBIG[2], e$MBIG[3], e$MBIG[4])
    MGP_LAB = c(2.1, 0.9, 0); IC_DROP = 2.64
    par(mar = MSQ, mgp = MGP_LAB)
    e$phase(e$exBI, sprintf("%s   %s, Tregs-ON bistable", tags[1], sp$lab),
            anchors = sp$bi, ic_marks = TRUE, ic_contour = TRUE, ic_brace = TRUE,
            ic_drop = IC_DROP, lab_cex = sp$lab_cex)
    par(mar = MSQ, mgp = MGP_LAB)
    e$phase(e$exMB, sprintf("%s   %s, Tregs-ON monostable", tags[2], sp$lab),
            anchors = sp$mb, ic_marks = TRUE, ic_contour = FALSE,
            ic_drop = IC_DROP, lab_cex = sp$lab_cex)
  })
  invisible(dev.off()); cat(sprintf("[fig] wrote %s\n", out_png))
}
