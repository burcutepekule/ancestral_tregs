## ============================================================================
## Fig_2.R -- Fig. 2: the efferocytosis knockout.
##
## Panels
##   A  phase portrait (Eh, M1) of one sterile example system, efferocytosis ON
##      (solid) and OFF (dashed, the knockout), 8 initial conditions, t = 0..4000
##   B  basin map of the same system over 80 initial conditions, ON and OFF
##   C  sterile:    final states in (Eh, M1)  -- OFF, ON high-burden, ON low-burden
##   D  pathogenic: final states in (Eh, M1)  -- same three groups
##   E  pathogenic: final states in (Eh, P_LP) -- same three groups
##   F  pathogenic: median Delta P_LP (ON - OFF) against beta_ros, 95% bootstrap CI
##   G  pathogenic: median Delta P_LP (ON - OFF) against r_L_P,    95% bootstrap CI
##
## "Low-burden basin" in this figure is NOT a fixed-point computation: a line
## (one unit under one arm, 21 final states at t = 12000) has one when its finite
## final composites span more than GAP_CUT = 20 composite units. Panels A/B
## re-integrate the example to t = 4000 only.
##
## in : ./data/effero_knockout.rds   (pipeline/06_effero_knockout.R)
##      ./data/tau_grid_v6_<n>.rds    (via lib/skeleton_defs.R)
## out: ./figures/main/Fig_2.png
##      ./source_data/Fig_2_numbers.csv   every count / endpoint number behind the panels
##      ./source_data/Fig_2_panelFG.csv   the plotted medians and bootstrap CIs of F and G
##      ./source_data/Fig_2_panelCDE.csv  per-group summaries and density peaks of C, D, E
## RUN: Rscript figure_scripts/Fig_2.R            (about 30 s)
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
rm(list = setdiff(ls(), "SRC"))

suppressMessages({library(data.table); library(deSolve); library(MASS)})
## skeleton_defs.R gives rhs_app (with the eff_ko multiplier on g), STRUCTURAL,
## EFFERO, KIN, STATE_ORDER, cmpf, ode_ext and tg keyed on (batch, row_ind).
## BRANCH only sets PL0/PLP0 for its own seed(), which this script does not use.
BRANCH = "sterile"
source("./lib/skeleton_defs.R")
FP_FIXED = list(K_E = STRUCTURAL$K_E, K_M = STRUCTURAL$K_M, K_T = STRUCTURAL$K_T)
fp_make_rhs = function() rhs_app
fp_rhs_params = function(par0, Kg, beta_ros, r_L_P, theta_treg = 0L) {
  p = c(STRUCTURAL, EFFERO)
  for (nm in KIN) p[[nm]] = par0[[nm]]
  p$Kg = Kg; p$beta_ros = beta_ros; p$r_L_P = r_L_P; p$r_LP_P = r_L_P
  p$theta_treg = theta_treg; p }
## hmax caps the solver's internal step (lsoda otherwise strides the stiff tail
## and can return NaN or a wrong basin). 5, as in the sweep.
HMAX = as.numeric(Sys.getenv("HMAX", "5"))


## ############################################################################
## knobs
## ############################################################################
## Relative heights of the right column's two rows:
##   RH_SC  final-state density panels (C/D/E)   RH_DRV  driver curves (F/G)
RH_DENS = 1.15; RH_SC = 1.215; RH_DRV = 1.135
CD_MAR_R = 0.5
## FS multiplies every cex (text only). FIG_POINTSIZE is the device base type
## size; FIG_SCALE the device size in pixels as a multiple of 3300 x HEIGHT.
FS = as.numeric(Sys.getenv("FONT_SIZE", "1.4"))
FIG_POINTSIZE = as.numeric(Sys.getenv("FIG_POINTSIZE", "15"))
FIG_SCALE = as.numeric(Sys.getenv("FIG_SCALE", "1.15"))
## PANEL_A_HEADS: direction on panel A's trajectories as arrowheads sitting on the
## curve (no shafts). See add_heads().
PANEL_A_HEADS = nzchar(Sys.getenv("PANEL_A_HEADS", "1"))
## ALT_LABELS: the larger typography (default). ALT_LABELS= (empty) gives the
## smaller sizes and prints n on C/D/E.
ALT_LABELS = nzchar(Sys.getenv("ALT_LABELS", "1"))
altv = function(base, alt) if (ALT_LABELS) alt else base
## C/D/E draw highest-density regions. HDR_P: mass fraction each band encloses,
## outermost first. HDR_A: its fill alpha. Groups smaller than N_CONTOUR_MIN are
## drawn as points.
HDR_P = c(0.90, 0.70, 0.50, 0.25)
HDR_A = c(0.10, 0.15, 0.22, 0.32)
N_CONTOUR_MIN = 25L

## ---- GAP_CUT: the low-burden-basin criterion, used by every panel ----------
## A line (the 21 final states of one unit under one arm) "has a low-burden
## basin" when its finite final composites span more than GAP_CUT composite units.
## comp = sqrt((K_E - Eh)^2 + P_LP^2); in the sterile branch P_LP = 0, so a gap of
## G is G units of epithelium. No fixed points are computed here.
## Sensitivity:  GAP_CUT=15 Rscript figure_scripts/Fig_2.R
GAP_CUT   = as.numeric(Sys.getenv("GAP_CUT", "20"))
has_lower = function(v) { v = v[is.finite(v)]; length(v) >= 2L && (max(v) - min(v)) > GAP_CUT }

## ---- panel A/B example --------------------------------------------------------
## Chosen from the data below (sterile, low-burden basin intact, none under
## knockout, widest spread). EX_ROW (= batch*1e5 + row_ind) / EX_KG / EX_BETA
## override the automatic pick.
ROW_IND = NA_integer_; Kg = NA_real_; beta_ros = NA_real_; r_L_P = 0
EH0_FAN = c(2, 16, 30, 44, 58, 72, 86, 97)      # panel A initial conditions
EH0_MAP = seq(1, 99, length.out = 80)           # panel B initial conditions
TMAX = 4000                                     # panels A/B integration horizon
## THR: Eh above which a panel-A/B run is coloured as the low-burden basin. Set
## from the data below: the midpoint between the smallest and largest final Eh of
## the 16 panel-A runs (8 ON + 8 OFF).
THR = NA_real_
LOW  = "#2E86FF"; HIGH  = "#FF7070"     # efferocytosis ON: lower / higher burden
LOWd = "#261CC1"; HIGHd = "#B02020"     # efferocytosis OFF
STATE_ORDER = c("PL","PLP","CLP","Eh","Ea","En","M0","M1","R","T0","Tr","D","S")
rhs = fp_make_rhs()

## initial condition: the line Eh = e, En = K_E - e, M0 = K_M, T0 = K_T, rest 0;
## pathogenic adds PL = 500, PLP = 10 (unused here: the example is sterile)
build_state = function(Eh0, patho) c(
  PL = if (patho) 500 else 0, PLP = if (patho) 10 else 0, CLP = 0, Eh = Eh0, Ea = 0,
  En = FP_FIXED$K_E - Eh0, M0 = FP_FIXED$K_M, M1 = 0, R = 0, T0 = FP_FIXED$K_T,
  Tr = 0, D = 0, S = 0)
## ode_ext: lsodar + extinction root, atol 1e-14, rtol 1e-8, hmax 5
## (lib/extinction_floor.R); output every 2 time units to TMAX
run = function(Eh0, eff_ko, theta = 0) {
  rp = c(fp_rhs_params(par0, Kg, beta_ros, r_L_P = r_L_P, theta_treg = theta), list(eff_ko = eff_ko))
  out = ode_ext(build_state(Eh0, r_L_P > 0)[STATE_ORDER], seq(0, TMAX, by = 2),
                function(t, y, p) list(rhs(t, y, p)[[1]]), rp, maxsteps = 1e5)
  list(Eh = out[, "Eh"], M1 = out[, "M1"])
}
endpt = function(z) c(tail(z$Eh, 1), tail(z$M1, 1))
bcol  = function(e) if (e > THR) LOW  else HIGH
bcold = function(e) if (e > THR) LOWd else HIGHd
add_arrow = function(x, y, col)
  arrows(x[1], y[1], x[2], y[2], length = 0.06, angle = 22, lwd = 2, col = col)
## ---- direction markers on the trajectories (PANEL_A_HEADS) ------------------
## Arrowheads only, placed by arc length in inches on the page, N_HEAD per
## trajectory, skipping any head within HEAD_SEP of one already drawn in the panel
## (converging trajectories would otherwise stack heads at the attractor).
HEAD_LEN = 0.16   # arrowhead length along the path, inches
HEAD_WID = 0.090  # arrowhead half-width at its base, inches
N_HEAD   = 3L     # heads per trajectory, before the spacing filter
HEAD_SEP = 0.30   # minimum distance between any two heads in the panel, inches
.head_drawn = new.env(parent = emptyenv())
head_reset = function() assign("p", matrix(numeric(0), nrow = 2), envir = .head_drawn)
head_reset()
## phase shifts the row of heads along the path by a fraction of the spacing, so
## ON and OFF heads of one initial condition alternate instead of stacking
add_heads = function(x, y, col, n = N_HEAD, len = HEAD_LEN, wid = HEAD_WID, phase = 0) {
  xi = grconvertX(x, "user", "inches"); yi = grconvertY(y, "user", "inches")
  ok = is.finite(xi) & is.finite(yi); xi = xi[ok]; yi = yi[ok]
  if (length(xi) < 2L) return(invisible())
  d = c(0, cumsum(sqrt(diff(xi)^2 + diff(yi)^2))); L = d[length(d)]
  if (!is.finite(L) || L < 1.5 * len) return(invisible())
  at = function(s) {                                   # point at arc length s
    i = max(1L, min(length(d) - 1L, findInterval(s, d)))
    f = if (d[i+1] > d[i]) (s - d[i]) / (d[i+1] - d[i]) else 0
    c(xi[i] + f * (xi[i+1] - xi[i]), yi[i] + f * (yi[i+1] - yi[i]))
  }
  ss = (seq_len(n) + phase) / (n + 1) * L              # interior, evenly spaced
  for (s in ss) {
    p = at(s)
    ## tangent from a chord half a head long either side
    u = at(min(s + len/2, L)) - at(max(s - len/2, 0))
    m = sqrt(sum(u^2)); if (!is.finite(m) || m == 0) next
    u = u / m; v = c(-u[2], u[1])                       # unit tangent, unit normal
    seen = get("p", envir = .head_drawn)
    if (ncol(seen) && min(sqrt(colSums((seen - p)^2))) < HEAD_SEP) next
    assign("p", cbind(seen, p), envir = .head_drawn)
    tri = cbind(p + u * len/2, p - u * len/2 + v * wid, p - u * len/2 - v * wid)
    polygon(grconvertX(tri[1, ], "inches", "user"),
            grconvertY(tri[2, ], "inches", "user"), col = col, border = "white", lwd = 0.6)
  }
}
callout = function(tx, ty, ex, ey, label, col, cex = 0.9 * FS) {
  segments(tx, ty, tx + 0.9 * (ex - tx), ty + 0.9 * (ey - ty), col = adjustcolor(col, .7), lwd = 1)
  text(tx, ty, label, col = col, font = 2, cex = cex)
}

## ---- data ------------------------------------------------------------------
## KO_RDS: one file, a comma-separated list or a glob. Chunk checkpoints
## (_partNNN.rds) are excluded, since their union is the assembled file. Pooled
## files are de-duplicated on (parameter set, Kg, eff_ko, Eh0). The figure as
## published reads ./data/effero_knockout.rds only.
IN_RDS = Sys.getenv("KO_RDS", "./data/effero_knockout.rds")
.kof = unlist(strsplit(IN_RDS, ",", fixed = TRUE))
.kof = trimws(.kof)
.kof = unlist(lapply(.kof, function(f) if (grepl("[*?]", f)) Sys.glob(f) else f))
.kof = .kof[!grepl("_part[0-9]+\\.rds$", .kof)]
if (!length(.kof)) stop("KO_RDS matched no files: ", IN_RDS)
.miss = .kof[!file.exists(.kof)]
if (length(.miss)) stop("no ", paste(.miss, collapse = ", "),
                        " -- run pipeline/06_effero_knockout.R first")
KD = rbindlist(lapply(.kof, function(f) {
  z = as.data.table(readRDS(f))
  cat(sprintf("[Fig. 2] %-42s %9s rows\n", f, format(nrow(z), big.mark = ",")))
  z }), fill = TRUE)
N_ROWS_READ = nrow(KD)
KD = KD[is.finite(Eh) & is.finite(PLP) & is.finite(comp)]

if (length(.kof) > 1L) {
  .UK = intersect(c("batch", "row_ind", "beta_ros", "r_L_P", "Kg", "eff_ko", "Eh0"),
                  names(KD))
  .n0 = nrow(KD); KD = unique(KD, by = .UK)
  cat(sprintf("[Fig. 2] pooled %d files: %s rows, %s after dropping %s duplicate arms\n",
              length(.kof), format(.n0, big.mark = ","),
              format(nrow(KD), big.mark = ","), format(.n0 - nrow(KD), big.mark = ",")))
}
## row_ind is per tau grid; replace it by batch*1e5 + row_ind, unique across
## grids, so that every expression keyed on row_ind means one kinetic vector.
KD[, row_ind0 := row_ind]
KD[, row_ind := batch * 100000L + row_ind]

## ---- the panel A/B example, chosen from the data ---------------------------
## Candidates: sterile systems (unit, Kg) with a low-burden basin intact and
## none under knockout; the one with the largest spread of final composites.
.cand = merge(
  KD[eff_ko == 1L & r_L_P == 0, .(on = has_lower(comp)),
     by = .(row_ind, beta_ros, Kg)][on == TRUE],
  KD[eff_ko == 0L & r_L_P == 0, .(off = has_lower(comp)),
     by = .(row_ind, beta_ros)][off == FALSE],
  by = c("row_ind","beta_ros"))
if (!nrow(.cand)) stop("no sterile set is bistable intact and monostable under knockout")
if (nzchar(Sys.getenv("EX_ROW"))) {
  ROW_IND = as.integer(Sys.getenv("EX_ROW"))
  Kg      = as.numeric(Sys.getenv("EX_KG", "20"))
  beta_ros = as.numeric(Sys.getenv("EX_BETA", "0.005"))
} else {
  .sep = KD[eff_ko == 1L & r_L_P == 0, .(sep = diff(range(comp))), by = .(row_ind, beta_ros, Kg)]
  .cand = merge(.cand, .sep, by = c("row_ind","beta_ros","Kg"))[order(-sep)]
  ROW_IND = .cand$row_ind[1]; Kg = .cand$Kg[1]; beta_ros = .cand$beta_ros[1]
}
.exb = ROW_IND %/% 100000L; .exr = ROW_IND %% 100000L
par0 = as.list(tg[.(.exb, .exr), nomatch = NULL][1])
if (!length(par0)) stop("no tau grid row for batch ", .exb, " row_ind ", .exr)
for (.nm in c("row_ind","batch","treg_eff","phi_instability")) par0[[.nm]] = NULL
cat(sprintf("[Fig. 2] panel A/B example: batch %d row_ind %d | Kg %g | beta_ros %g | r_L_P %g\n",
            .exb, .exr, Kg, beta_ros, r_L_P))
cat(sprintf("        %s sterile candidates were bistable intact and monostable under knockout\n",
            format(nrow(.cand), big.mark = ",")))
K_E = 100; FLOOR = 0.1
K_M_CAP = FP_FIXED$K_M            # macrophage carrying capacity, for panel axis caps
key = c("row_ind","beta_ros","r_L_P","pathogenic","Eh0")

## ---- the unit gate for panels F/G --------------------------------------------
## A unit qualifies if its intact line has a low-burden basin at one or more Kg.
## All 7 Kg arms and all 21 initial conditions of a qualifying unit are kept.
GATE_SETS = !identical(tolower(Sys.getenv("GATE_SETS", "TRUE")), "false")
QSET = unique(KD[eff_ko == 1L, .(b = has_lower(comp)), by = .(row_ind, beta_ros, r_L_P, Kg)
                 ][b == TRUE, .(row_ind, beta_ros, r_L_P)])
KDq  = if (GATE_SETS) merge(KD, QSET, by = c("row_ind","beta_ros","r_L_P")) else copy(KD)
cat(sprintf("[gate] GAP_CUT=%g GATE_SETS=%s | %d of %d units qualify | %s of %s rows kept\n",
            GAP_CUT, GATE_SETS, nrow(QSET),
            uniqueN(KD[, .(row_ind, beta_ros, r_L_P)]),
            format(nrow(KDq), big.mark = ","), format(nrow(KD), big.mark = ",")))

koff = KDq[arm == "KO",  .(Eh = Eh[1], PLP = PLP[1]),           by = key]

## ---- panels F/G: paired Delta P_LP --------------------------------------------
## dPLP per (unit, Eh0, Kg) = max(P_LP intact at Kg, 0.1) - max(P_LP knockout, 0.1)
## (the knockout arm is Kg-independent).
konK = KDq[eff_ko == 1L, .(PLP_on = median(PLP)), by = c(key, "Kg")]
PK   = merge(konK, koff[, .(row_ind, beta_ros, r_L_P, pathogenic, Eh0, PLP_off = PLP)], by = key)
PK[, dPLP := pmax(PLP_on, FLOOR) - pmax(PLP_off, FLOOR)]
pPK = PK[pathogenic == TRUE]
pPK[, src := "production"]

## ---- extra pathogenic units for F/G only (pipeline/06b_effero_fg_extra_units.R)
## Drawn to fill the thin (beta_ros, r_L_P) cells; gated and paired exactly as
## above. They enter panels F and G and nothing else. FG_EXTRA_RDS is a
## comma-separated list (missing files are skipped); "" skips them all.
FG_EXTRA_RDS = Sys.getenv("FG_EXTRA_RDS", "./data/effero_fg_extra.rds,./data/effero_fg_extra2.rds")
.fgx = Filter(file.exists, trimws(strsplit(FG_EXTRA_RDS, ",")[[1]]))
N_FG_EXTRA_UNITS = 0L; N_FG_EXTRA_DRAWN = 0L
if (length(.fgx)) {
  KX = rbindlist(lapply(.fgx, function(f) as.data.table(readRDS(f))), fill = TRUE)
  KX = KX[is.finite(Eh) & is.finite(PLP) & is.finite(comp)]
  stopifnot(!anyDuplicated(KX, by = c("batch", "row_ind", "beta_ros", "r_L_P", "Kg", "eff_ko", "Eh0")))
  KX[, row_ind := batch * 100000L + row_ind]
  QX = unique(KX[eff_ko == 1L, .(b = has_lower(comp)), by = .(row_ind, beta_ros, r_L_P, Kg)
                 ][b == TRUE, .(row_ind, beta_ros, r_L_P)])
  KXq = merge(KX, QX, by = c("row_ind", "beta_ros", "r_L_P"))
  kxoff = KXq[arm == "KO", .(PLP_off = PLP[1]), by = key]
  kxon  = KXq[eff_ko == 1L, .(PLP_on = median(PLP)), by = c(key, "Kg")]
  PX = merge(kxon, kxoff, by = key)[pathogenic == TRUE]
  PX[, dPLP := pmax(PLP_on, FLOOR) - pmax(PLP_off, FLOOR)]
  PX[, src := "fg_extra"]
  stopifnot(!nrow(merge(unique(PX[, .(row_ind, beta_ros, r_L_P)]),
                        unique(pPK[, .(row_ind, beta_ros, r_L_P)]), by = c("row_ind", "beta_ros", "r_L_P"))))
  N_FG_EXTRA_UNITS = uniqueN(PX[, .(row_ind, beta_ros, r_L_P)])
  N_FG_EXTRA_DRAWN = uniqueN(KX[, .(row_ind, beta_ros, r_L_P)])
  cat(sprintf("[F/G] %s: %d of %d extra units have a low-burden basin; %s pairs added\n",
              paste(.fgx, collapse = ", "), N_FG_EXTRA_UNITS, N_FG_EXTRA_DRAWN,
              format(nrow(PX), big.mark = ",")))
  pPK = rbind(pPK, PX, fill = TRUE)
}
## ---- the plotted statistic in F/G ---------------------------------------------
## FG_STAT = "balanced": at each level of the driver on the x axis
##   (beta_ros in F, r_L_P in G), the weighted median of dPLP in which every
##   stratum of the other driver (r_L_P in F, beta_ros in G) has equal total
##   weight: a pair has weight 1 / (number of pairs in its stratum). Every unit
##   carries all 7 Kg x 21 Eh0, so Kg and Eh0 need no weight. This is the
##   statistic that equal subsampling of the strata estimates, using all units.
## FG_STAT = "kg": median over the 7 Kg of the per-Kg median of dPLP (strata
##   weighted by their number of units; the statistic before 2026-10-01).
## FG_STAT = "subsample" (default): every (beta_ros, r_L_P) cell is subsampled to the same
##   number of units, FG_K (default: the smallest cell over both panels), and the
##   plain median of dPLP is taken over the pooled subsample. The point is the
##   median of NSUB such medians (subsamples without replacement); the band is the
##   2.5-97.5% range of NBOOT medians of subsamples drawn with replacement.
FG_STAT = Sys.getenv("FG_STAT", "subsample")
STRAT   = c(beta_ros = "r_L_P", r_L_P = "beta_ros")   # stratum variable per driver
wmed = function(x, w) { o = order(x); x = x[o]; w = w[o]; x[which(cumsum(w) >= sum(w) / 2)[1L]] }
fg_stat = function(dd, xcol) {
  if (identical(FG_STAT, "kg")) return(median(dd[, .(m = median(dPLP)), by = Kg]$m))
  st = dd[[STRAT[[xcol]]]]
  wmed(dd$dPLP, 1 / as.numeric(table(st)[as.character(st)]))
}
## per-Kg value of the statistic at each driver value (BAND = "kg" and the table)
mB = pPK[, .(m = fg_stat(.SD, "beta_ros")), by = .(Kg, beta_ros)]
mR = pPK[, .(m = fg_stat(.SD, "r_L_P")),    by = .(Kg, r_L_P)]

## ---- the band in F/G ----------------------------------------------------------
## BAND = "boot" (default): 95% percentile bootstrap interval on the plotted
##   statistic, NBOOT resamples of the UNITS within each driver level (the 21 x 7
##   rows of one unit share a parameter vector and are not independent). With
##   FG_STAT = "balanced" the units are resampled within each stratum, so every
##   replicate keeps the stratum sizes.
## BAND = "kg": min-max envelope of the per-Kg values (a spread, not an
##   uncertainty).
BAND  = Sys.getenv("BAND", "boot")
NBOOT = as.integer(Sys.getenv("NBOOT", "400"))
pPK[, .set := paste(row_ind, beta_ros, r_L_P)]

FG_K = NA_integer_
if (identical(FG_STAT, "subsample")) {
  FG_K = as.integer(Sys.getenv("FG_K", min(unique(pPK[, .(.set, beta_ros, r_L_P)])[, .N, by = .(beta_ros, r_L_P)]$N)))
  NSUB = as.integer(Sys.getenv("NSUB", "1000"))
}
sub_draw = function(S, U, xcol, repl) median(S[J(U[, .(.set = .set[sample.int(.N, FG_K, replace = repl)]),
                                                    by = c(STRAT[[xcol]])]$.set), dPLP, allow.cartesian = TRUE])
if (identical(FG_STAT, "subsample")) {       # per-Kg values of the statistic (table only), 200 subsamples each
  set.seed(as.integer(Sys.getenv("BOOT_SEED", "1")) + 1L)
  perkg = function(xcol) pPK[, {
    S = copy(.SD); setkey(S, .set); U = unique(S[, c(".set", STRAT[[xcol]]), with = FALSE])
    .(m = median(vapply(1:200, function(i) sub_draw(S, U, xcol, FALSE), numeric(1)))) }, by = c("Kg", xcol)]
  mB = perkg("beta_ros"); mR = perkg("r_L_P")
}
sub_band = function(D, xcol) {
  xs = sort(unique(D[[xcol]]))
  rbindlist(lapply(xs, function(xv) {
    S = D[get(xcol) == xv]; setkey(S, .set)
    U = unique(S[, c(".set", STRAT[[xcol]]), with = FALSE])
    pts  = vapply(seq_len(NSUB),  function(i) sub_draw(S, U, xcol, FALSE), numeric(1))
    reps = vapply(seq_len(NBOOT), function(i) sub_draw(S, U, xcol, TRUE),  numeric(1))
    data.table(x = xv, med = median(pts), n_sets = nrow(U),
               n_min_stratum = min(U[, .N, by = c(STRAT[[xcol]])]$N),
               lo = unname(quantile(reps, 0.025)), hi = unname(quantile(reps, 0.975)))
  }))[, setnames(.SD, "x", xcol)] }

boot_band = function(D, xcol) {
  if (identical(FG_STAT, "subsample")) return(sub_band(D, xcol))
  xs = sort(unique(D[[xcol]]))
  out = rbindlist(lapply(xs, function(xv) {
    S = D[get(xcol) == xv]
    U = unique(S[, c(".set", STRAT[[xcol]]), with = FALSE])
    reps = vapply(seq_len(NBOOT), function(i) {
      pick = if (identical(FG_STAT, "kg")) data.table(.set = sample(U$.set, nrow(U), replace = TRUE))
             else U[, .(.set = .set[sample.int(.N, .N, replace = TRUE)]), by = c(STRAT[[xcol]])][, .(.set)]
      fg_stat(merge(pick, S, by = ".set", allow.cartesian = TRUE), xcol) }, numeric(1))
    data.table(x = xv, med = fg_stat(S, xcol), n_sets = nrow(U),
               n_min_stratum = min(U[, .N, by = c(STRAT[[xcol]])]$N),
               lo = unname(quantile(reps, 0.025, na.rm = TRUE)),
               hi = unname(quantile(reps, 0.975, na.rm = TRUE))) }), fill = TRUE)
  setnames(out, "x", xcol)[] }

if (identical(BAND, "boot")) {
  set.seed(as.integer(Sys.getenv("BOOT_SEED", "1")))
  cB = boot_band(pPK, "beta_ros"); cR = boot_band(pPK, "r_L_P")
  cat(sprintf("[F/G] statistic = %s | band = 95%% bootstrap CI, %d resamples of the units\n", FG_STAT, NBOOT))
  if (identical(FG_STAT, "subsample")) cat(sprintf("      subsample: %d units per (beta_ros, r_L_P) cell, %d subsamples for the point\n", FG_K, NSUB))
  cat("      smallest stratum per level -- beta_ros: ", paste(cB$n_min_stratum, collapse = " "), " | r_L_P: ", paste(cR$n_min_stratum, collapse = " "), "\n", sep = "")
  cat("      units per level -- beta_ros: ", paste(cB$n_sets, collapse = " "), "\n", sep = "")
  cat("      units per level -- r_L_P   : ", paste(cR$n_sets, collapse = " "), "\n", sep = "")
} else {
  cB = mB[, .(med = median(m), lo = min(m), hi = max(m)), by = beta_ros][order(beta_ros)]
  cR = mR[, .(med = median(m), lo = min(m), hi = max(m)), by = r_L_P][order(r_L_P)]
  cat("[F/G] band = min-max across the Kg arms (a spread, not an uncertainty)\n")
}

## panel letters sit LETTER_DROP lines below the top of their figure region
LETTER_DROP = 1.6
## fit_cex: the given cex if the string fits `avail` inches, else the largest that
## does. fit_avail: the room a label has (to the edge of the figure region).
fit_cex = function(s, cex, avail) {
  w = strwidth(s, cex = cex, units = "inches")
  if (!is.finite(w) || w <= 0 || w <= avail) cex else cex * avail / w
}
FIT_PAD = 0.06   # inches
fit_avail = function(centred = TRUE) {
  if (!centred) return(par("pin")[1] + par("mar")[4] * par("csi") - FIT_PAD)
  c0 = par("mar")[2] * par("csi") + par("pin")[1] / 2   # plot centre, from the figure's left edge
  2 * min(c0, par("fin")[1] - c0) - 2 * FIT_PAD
}
tagit = function(lab, drop = LETTER_DROP) {
  mtext(
    lab,
    side = 3,
    line = par("mar")[3] - drop,
    adj = 0,
    at = par("usr")[1],
    font = 2,
    cex = 1.15 * FS
  )
}

## ---- panels A/B: trajectories and basin map (layout-independent) --------------
trajON  = lapply(EH0_FAN, function(e) run(e, 1))
trajOFF = lapply(EH0_FAN, function(e) run(e, 0))
.ends = c(sapply(trajON, function(z) tail(z$Eh, 1)), sapply(trajOFF, function(z) tail(z$Eh, 1)))
THR = (min(.ends, na.rm = TRUE) + max(.ends, na.rm = TRUE)) / 2
cat(sprintf("[Fig. 2] example row %d, Kg=%g, beta_ros=%g | final Eh %.1f..%.1f -> THR=%.1f\n",
            ROW_IND, Kg, beta_ros, min(.ends, na.rm=TRUE), max(.ends, na.rm=TRUE), THR))
if (!is.finite(THR) || diff(range(.ends, na.rm = TRUE)) < 2)
  cat("[Fig. 2] *** this example is NOT bistable on the panel-A initial conditions; A/B will be flat.\n")
ymax = max(sapply(c(trajON, trajOFF), function(z) max(z$M1))) * 1.06
loE = endpt(run(max(EH0_FAN), 1)); hiE = endpt(run(min(EH0_FAN), 1)); ofE = endpt(run(min(EH0_FAN), 0))
onV  = sapply(EH0_MAP, function(e) if (tail(run(e, 1)$Eh, 1) > THR) 1L else 2L)
offV = sapply(EH0_MAP, function(e) if (tail(run(e, 0)$Eh, 1) > THR) 3L else 4L)

COL_B = "#4E1F6E"; COL_R = "#0C513B"      # F: beta_ros = purple, G: r_L_P = green

## ---- the figure ---------------------------------------------------------------
## Left column: A over B. Right column: C/D/E over F/G.
## FIG_TALL stretches the whole figure vertically; B_FRAC is B's share of it.
FIG_TALL = as.numeric(Sys.getenv("FIG_TALL", "1.305"))
B_FRAC   = 0.21
render_fig = function(OUTFILE) {
  hs     = RH_SC + RH_DRV
  HEIGHT = round(596 * hs * FIG_TALL)                    # 596 px per height unit
  eBot   = 1 - RH_SC / hs                                # bottom of C/D/E = top of F/G
  png(OUTFILE, width = round(3300 * FIG_SCALE), height = round(HEIGHT * FIG_SCALE),
      res = 200, pointsize = FIG_POINTSIZE)
  par(xaxs = "i", cex.axis = 0.95 * FS, cex.lab = 1.05 * FS)

## ---- A: phase portrait ----
par(fig = c(0, 0.5, B_FRAC, 1.0), mar = c(4.4, 4.8, 2.2, 1.5),
    cex.lab = altv(1.05 * FS, 1.50 * FS), new = FALSE)
plot(NA, xlim = c(-2, 102), ylim = c(-0.6, ymax),
     xlab = expression("" * E[h]),
     ylab = expression("" * M[1]), main = "")
grid(col = "grey90", lty = 1)
.dir = function(x, y, col, phase)
  if (PANEL_A_HEADS) add_heads(x, y, col, phase = phase) else add_arrow(x, y, col)
head_reset()                       # the spacing filter is per panel
## the OFF high-burden state is grey45, the OFF colour in C/D/E
.acold = function(e) if (e > THR) LOWd else "grey45"
for (k in seq_along(EH0_FAN)) {
  zoff = trajOFF[[k]]; co = .acold(tail(zoff$Eh,1))
  lines(zoff$Eh, zoff$M1, col = adjustcolor(co, .9), lwd = 1.5, lty = 2); .dir(zoff$Eh, zoff$M1, co, 0.5)
  zon  = trajON[[k]];  cn = bcol(tail(zon$Eh,1))
  lines(zon$Eh,  zon$M1,  col = adjustcolor(cn, .9), lwd = 1.5, lty = 1); .dir(zon$Eh, zon$M1, cn, 0.0)
  points(EH0_FAN[k], 0, pch = 4, col = "black", cex = 1.5, lwd = 2.5)   # initial condition (black cross)
}
points(loE[1], loE[2], pch = 21, bg = bcol(loE[1]),  col = "white", cex = 2.1, lwd = 1.6)
points(hiE[1], hiE[2], pch = 21, bg = bcol(hiE[1]),  col = "white", cex = 2.1, lwd = 1.6)
points(ofE[1], ofE[2], pch = 21, bg = .acold(ofE[1]), col = "white", cex = 2.1, lwd = 1.6)
callout(48, ymax * 0.72, hiE[1], hiE[2], "High-burden attractor\n(efferocytosis ON)",  HIGH)
callout(83, ymax * 0.50, loE[1], loE[2], "Low-burden attractor\n(efferocytosis ON)",   LOW)
callout(38, ymax * 0.97, ofE[1], ofE[2], "High-burden attractor\n(efferocytosis OFF)", "grey45")
legend("topright", bty = "n", cex = 0.92 * FS, seg.len = 2.6,
       legend = c("Efferocytosis ON", "Efferocytosis OFF"),
       col = c("#555555", "#555555"), lwd = c(2.2, 2.2), lty = c(1, 2))
tagit("A")
par(cex.lab = 1.05 * FS)

## ---- B: basin map (same plot width as A) ----
par(fig = c(0, 0.5, 0.0, B_FRAC), mar = c(4.0, 4.8, 1.8, 1.5), new = TRUE)
ZB = cbind(offV, onV)
image(EH0_MAP, 1:2, ZB, col = c(LOW, HIGH, LOWd, "grey45"), breaks = c(.5,1.5,2.5,3.5,4.5),
      xlim = c(-2, 102), ylim = c(0.5, 2.5),
      xlab = "", ylab = "", yaxt = "n", main = "")
.bxl = expression("Initial epithelial health  " * E[h] * "(0)")
title(xlab = .bxl, cex.lab = fit_cex(.bxl, 1.05 * FS, fit_avail(TRUE)))
axis(2, at = 1:2, labels = c("OFF", "ON"), las = 1, tick = FALSE, cex.axis = altv(0.78 * FS, 1.02 * FS))
mtext("Efferocytosis", side = 2, line = 3.7, cex = altv(0.85 * FS, 0.92 * FS))
abline(h = 1.5, col = "white", lwd = 3)
onlow = onV == 1; i0 = which(onlow)[1]
if (any(onlow) && !all(onlow) && !is.na(i0) && i0 > 1) {
  s = (EH0_MAP[i0-1] + EH0_MAP[i0]) / 2
  segments(s, 2-0.47, s, 2+0.47, col = "white", lwd = 2.4) }
## a label inside a coloured block is fitted to the block's width and height
.bcex = function(lab, xs, want = 1.20 * FS, air = 0.86, row_h = 0.94, air_h = 0.80) {
  cx = fit_cex(lab, want, air * diff(range(xs)) / diff(par("usr")[1:2]) * par("pin")[1])
  h_in = strheight(lab, cex = cx, units = "inches"); h_ok = air_h * row_h / diff(par("usr")[3:4]) * par("pin")[2]
  if (is.finite(h_in) && h_in > h_ok) cx * h_ok / h_in else cx
}
if (any(onV == 2)) text(mean(EH0_MAP[onV == 2]), 2, "High-burden basin", col = "white", font = 2,
                        cex = altv(0.70 * FS, .bcex("High-burden basin", EH0_MAP[onV == 2])))
if (any(onV == 1)) text(mean(EH0_MAP[onV == 1]), 2, "Low-burden basin",  col = "white", font = 2,
                        cex = altv(0.70 * FS, .bcex("Low-burden basin", EH0_MAP[onV == 1])))
if (!any(offV == 3)) text(50, 1, "No low-burden basin (monostable)", col = "white", font = 2,
                          cex = altv(0.74 * FS, .bcex("No low-burden basin (monostable)", EH0_MAP)))
tagit("B")

## ---- C / D / E: OFF -> ON final-state densities -------------------------------
## Sterile contributes one panel (M1); pathogenic two (M1 and P_LP).
.scTop = 0.985
.gut = 0.008
.pw = (1.0 - 0.5 - 2 * .gut) / 3
.px = 0.5 + (0:2) * (.pw + .gut)
draw_shift(FALSE, "M1",  expression(M[1]),  .px[1], .px[1] + .pw, eBot, .scTop, "C", "Sterile",    FALSE)
draw_shift(TRUE,  "M1",  expression(M[1]),  .px[2], .px[2] + .pw, eBot, .scTop, "D", "Pathogenic", FALSE)
draw_shift(TRUE,  "PLP", expression(P[LP]), .px[3], .px[3] + .pw, eBot, .scTop, "E", "Pathogenic", FALSE)
## one shared legend, horizontal, in the top margin above C/D/E
par(fig = c(0.5, 1.0, eBot, .scTop), mar = c(0, 0, 0, 0), new = TRUE)
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
.fig_h = par("fin")[2]                                                 # this band's height, inches
.yleg  = 1 - (SHIFT_MAR_TOP - SHIFT_LEG_LINE) * par("csi") / .fig_h    # NDC within the band
.leg   = c("OFF attractor", "ON high-burden attractor", "ON low-burden attractor")
.lw    = function(cx) sum(strwidth(.leg, cex = cx, units = "inches")) + 3 * 2.5 * strwidth("M", cex = cx, units = "inches")
.lcx   = 1.00 * FS; if (.lw(.lcx) > 0.94 * par("pin")[1]) .lcx = .lcx * 0.94 * par("pin")[1] / .lw(.lcx)
legend(0.5, .yleg, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n", cex = .lcx,
       pch = 23, col = "white", pt.lwd = 1.2, pt.cex = 1.5, pt.bg = c("grey45", HIGH, LOW),
       legend = .leg, x.intersp = 0.8, text.width = NA)
cat(sprintf("[C/D/E legend] cex = %.2f | panel plot box %.2f x %.2f in\n", .lcx, .pw_in[1], .pw_in[2]))

## ---- F / G: median Delta P_LP against each driver -------------------------
## Separate x axes (beta_ros and r_L_P differ by two orders of magnitude), shared y.
COL_B = "#4E1F6E"; COL_R = "#0C513B"      # beta_ros = purple, r_L_P = green
yr = range(0, cB$lo, cB$hi, cR$lo, cR$hi, na.rm = TRUE)
draw_driver = function(d, xcol, col, pch) {
  x = d[[xcol]]
  polygon(c(x, rev(x)), c(d$lo, rev(d$hi)), col = adjustcolor(col, .12), border = NA)
  lines(x, d$lo, col = adjustcolor(col, .75), lwd = 1, lty = 2)
  lines(x, d$hi, col = adjustcolor(col, .75), lwd = 1, lty = 2)
  lines(x, d$med, col = col, lwd = 1.9); points(x, d$med, col = col, pch = pch, cex = 1.05, lwd = 1.4)
}
leg_driver = function(where = "bottomright")
  legend(where, bty = "o", bg = "white", box.col = "grey80", box.lwd = 0.8,
         cex = altv(0.72 * FS, 0.98 * FS), seg.len = altv(2.0, 2.4),
         legend = if (identical(BAND, "boot")) "95% CI" else expression(K[g] ~ "range"),
         col = "grey45", lty = 2, lwd = altv(1.2, 2.6))

par(fig = c(0.500, 0.755, 0.0, eBot), mar = c(4.0, 5, 2.4, 0.6), mgp = c(2.2, 0.7, 0),
    cex.axis = altv(0.95 * FS, 1.18 * FS), new = TRUE)
.xlF = expression(beta[r])
plot(NA, xlim = range(cB$beta_ros), ylim = yr, xaxs = "r", xaxt = "n",
     xlab = if (ALT_LABELS) "" else .xlF,
     ylab = expression(paste("median  ", Delta, P[LP], "  (ON - OFF)")), main = "")
if (ALT_LABELS) title(xlab = .xlF, cex.lab = 1.45 * FS)
## beta_ros ticks in units of 1e-3
.bt = pretty(range(cB$beta_ros), 4)
axis(1, at = .bt, labels = formatC(.bt * 1e3, format = "g", digits = 2))
mtext(expression("(" %*% 10^-3 * ")"), side = 1, line = altv(1.9, 3.1), adj = 1,
      cex = altv(0.72 * FS, 0.98 * FS))
grid(col = "grey90", lty = 1, nx = NA, ny = NULL); abline(v = .bt, col = "grey90")
abline(h = 0, lty = 3, col = "grey55")
draw_driver(cB, "beta_ros", COL_B, 19)
leg_driver("topright")
box(); tagit("F", drop = 1.2)

.kL = 0.75 + (4.4 - 0.6) * par("csi") / par("din")[1]
par(fig = c(.kL, 1.000, 0.0, eBot), mar = c(4.0, 0.6, 2.4, 1.0), mgp = c(2.2, 0.7, 0),
    cex.axis = altv(0.95 * FS, 1.18 * FS), new = TRUE)
.xlG = expression(r[P])
plot(NA, xlim = range(cR$r_L_P), ylim = yr, xaxs = "r", yaxt = "n",
     xlab = if (ALT_LABELS) "" else .xlG, ylab = "", main = "")
if (ALT_LABELS) title(xlab = .xlG, cex.lab = 1.45 * FS)
grid(col = "grey90", lty = 1); abline(h = 0, lty = 3, col = "grey55")
draw_driver(cR, "r_L_P", COL_R, 8)
leg_driver("topleft")
box(); tagit("G", drop = 1.2)

  dev.off()
  cat(sprintf("wrote %s\n", OUTFILE))
}

## ============================================================================
## C/D/E data: per system, the low- and high-burden final states and the
## knockout outcome of the same unit
## ----------------------------------------------------------------------------
## system = (unit, Kg) for the intact arm. lb = its minimum-composite final state,
## hb = its maximum-composite final state (over the 21 initial conditions).
## offH = the unit's knockout outcome, the MEAN over its knockout final states.
## Only systems whose spread hb - lb exceeds GAP_CUT are drawn.
## ============================================================================
onH = KD[eff_ko == 1L]
onH[, sys := paste(row_ind, beta_ros, r_L_P, Kg, sep = "|")]
setorder(onH, sys, comp)
lb = onH[onH[, .I[1L], by = sys]$V1]     # min composite per system = low-burden final state
hb = onH[onH[, .I[.N], by = sys]$V1]     # max composite per system = high-burden final state
offH = KD[arm == "KO", .(Eh = mean(Eh), M1 = mean(M1), PLP = mean(PLP)),
          by = .(row_ind, beta_ros, r_L_P, pathogenic)]
offH[, sk := paste(row_ind, beta_ros, r_L_P, sep = "|")]

.gap = merge(lb[, .(sys, pathogenic, c_l = comp, Eh_l = Eh, PLP_l = PLP)],
             hb[, .(sys, c_h = comp, Eh_h = Eh, PLP_h = PLP)], by = "sys")
.gap[, hi_lo_gap := c_h - c_l]
KEEP_SYS = .gap[hi_lo_gap > GAP_CUT, sys]
N_KEEP   = c("FALSE" = 0L, "TRUE" = 0L)
for (.p in c(FALSE, TRUE)) {
  .x = .gap[pathogenic == .p]; .k = .x[hi_lo_gap > GAP_CUT]
  N_KEEP[[as.character(.p)]] = nrow(.k)
  cat(sprintf("[C/D/E] %-11s: %d of %d systems have spread > %g (%.1f%%) | low-burden state has higher Eh AND no higher P_LP: %.1f%%\n",
              if (.p) "pathogenic" else "sterile", nrow(.k), nrow(.x), GAP_CUT,
              100 * nrow(.k) / max(nrow(.x), 1L),
              100 * mean(.k$Eh_l > .k$Eh_h & .k$PLP_l <= .k$PLP_h)))
}
## axis limits: 0.5% and 99.5% quantiles of the plotted points, padded by 6% of
## that span, clamped to [0, K_E] for Eh and [0, K_M] or [0, K_E] for M1 / P_LP
lim_pad = function(v, lo = -Inf, hi = Inf, q = 0.995) {
  v = v[is.finite(v)]
  if (!length(v)) return(c(0, 1))
  ab = as.numeric(stats::quantile(v, c(1 - q, q), na.rm = TRUE))
  a = ab[1]; b = ab[2]
  if (!is.finite(a) || !is.finite(b) || b <= a) { a = min(v); b = a + max(abs(a), 1) * 0.05 }
  d = 0.06 * (b - a)
  c(max(lo, a - d), min(hi, b + d))
}
best_corner = function(x, y, xl, yl) {
  ok = is.finite(x) & is.finite(y); x = x[ok]; y = y[ok]
  if (!length(x)) return("topright")
  fx = (x - xl[1]) / diff(xl); fy = (y - yl[1]) / diff(yl)
  cn = c(topleft     = sum(fx < 0.45 & fy > 0.55),
         topright    = sum(fx > 0.55 & fy > 0.55),
         bottomleft  = sum(fx < 0.45 & fy < 0.45),
         bottomright = sum(fx > 0.55 & fy < 0.45))
  names(cn)[which.min(cn)]
}
## top margin of C/D/E in lines; subtitle at 0.35, legend at SHIFT_LEG_LINE,
## letter at SHIFT_LET_LINE
SHIFT_MAR_TOP  = 4.0
SHIFT_LEG_LINE = 2.15
SHIFT_LET_LINE = 3.3
PEAKS = list()       # density peaks per panel and group, for the numbers export
GROUPS = list()      # the plotted points per panel and group, for the numbers export
draw_shift = function(pat, yv, ylab_expr, xL, xR, yb, yt, letter, ttl, showleg) {
  blb = lb[sys %in% KEEP_SYS & pathogenic == pat]
  bhb = hb[sys %in% KEEP_SYS & pathogenic == pat]
  blb[, sk := paste(row_ind, beta_ros, r_L_P, sep = "|")]
  offp = offH[pathogenic == pat]; ob = offp[match(blb$sk, offp$sk)]
  XA = c(ob$Eh, bhb$Eh, blb$Eh); YA = c(ob[[yv]], bhb[[yv]], blb[[yv]])
  xlim = lim_pad(XA, lo = 0, hi = K_E)
  ylim = lim_pad(YA, lo = 0, hi = if (yv == "PLP") K_E else K_M_CAP)
  par(fig = c(xL, xR, yb, yt), mar = c(3.4, 3.4, SHIFT_MAR_TOP, 0.4), mgp = c(2.1, 0.6, 0), new = TRUE)
  plot(NA, xlim = xlim, ylim = ylim, xlab = expression(E[h]), ylab = ylab_expr, main = "", xaxt = "n")
  axis(1, at = seq(0, 100, by = 40))   # regular ticks; R would drop overlapping labels at by = 20
  .pw_in <<- par("pin")
  grid(col = "grey92")
  ## Highest-density regions of each group's own mass: MASS::kde2d (default
  ## bandwidth) on a 120 x 120 grid over the window extended by 10% per side; the
  ## band enclosing fraction p is bounded by the contour at the (1 - p) quantile of
  ## the density evaluated at the group's points.
  hdr = function(x, y, col) {
    ok = is.finite(x) & is.finite(y); x = x[ok]; y = y[ok]
    if (length(x) < N_CONTOUR_MIN || length(unique(x)) < 3L || length(unique(y)) < 3L) {
      points(x, y, pch = 16, col = adjustcolor(col, 0.45), cex = 0.5); return(invisible())
    }
    ex = 0.10 * diff(xlim); ey = 0.10 * diff(ylim)
    k = tryCatch(MASS::kde2d(x, y, n = 120,
                             lims = c(xlim[1] - ex, xlim[2] + ex, ylim[1] - ey, ylim[2] + ey)),
                 error = function(e) NULL)
    if (is.null(k)) { points(x, y, pch = 16, col = adjustcolor(col, 0.35), cex = 0.45); return(invisible()) }
    ix = pmin(pmax(findInterval(x, k$x), 1L), length(k$x))
    iy = pmin(pmax(findInterval(y, k$y), 1L), length(k$y))
    dz = k$z[cbind(ix, iy)]
    for (i in seq_along(HDR_P)) {
      lv = as.numeric(quantile(dz, 1 - HDR_P[i], na.rm = TRUE))
      if (!is.finite(lv) || lv <= 0) next
      for (p in contourLines(k$x, k$y, k$z, levels = lv))
        polygon(p$x, p$y, col = adjustcolor(col, HDR_A[i]), border = NA)
    }
    lv = as.numeric(quantile(dz, 1 - HDR_P[1], na.rm = TRUE))
    if (is.finite(lv) && lv > 0)
      for (p in contourLines(k$x, k$y, k$z, levels = lv))
        lines(p$x, p$y, col = adjustcolor(col, 0.65), lwd = 1.1)
  }
  hdr(ob$Eh, ob[[yv]], "grey45"); hdr(bhb$Eh, bhb[[yv]], HIGH); hdr(blb$Eh, blb[[yv]], LOW)

  ## ---- group markers: density modes ------------------------------------------
  ## Local maxima of the same kernel surface (8-neighbourhood), at least PEAK_FRAC
  ## of the highest, more than PEAK_SEP apart in axis-normalised units. Every point
  ## is assigned to its nearest peak (axis-normalised distance); a peak is kept if
  ## it holds >= MASS_FRAC of the group (else the largest one is kept); peaks are
  ## ranked by assigned mass; the largest is the filled marker, at most MAX_LOBES
  ## are drawn. CENTRE = "median" uses componentwise medians instead.
  CENTRE    = Sys.getenv("CENTRE", "mode")
  MASS_FRAC = as.numeric(Sys.getenv("MASS_FRAC", "0.25"))
  MAX_LOBES = as.integer(Sys.getenv("MAX_LOBES", "2"))
  PEAK_FRAC = as.numeric(Sys.getenv("PEAK_FRAC", "0.05"))
  PEAK_SEP  = as.numeric(Sys.getenv("PEAK_SEP", "0.12"))
  kde_peaks = function(xx, yy) {
    ok = is.finite(xx) & is.finite(yy); xx = xx[ok]; yy = yy[ok]
    fallback = list(main = c(median(xx), median(yy)), sec = NULL)
    if (identical(CENTRE, "median")) return(fallback)
    if (length(xx) < N_CONTOUR_MIN || length(unique(xx)) < 3L || length(unique(yy)) < 3L)
      return(fallback)
    ex = 0.10 * diff(xlim); ey = 0.10 * diff(ylim)
    k = tryCatch(MASS::kde2d(xx, yy, n = 120,
                             lims = c(xlim[1] - ex, xlim[2] + ex, ylim[1] - ey, ylim[2] + ey)),
                 error = function(e) NULL)
    if (is.null(k)) return(fallback)
    z = k$z; nx = nrow(z); ny = ncol(z)
    core = z[2:(nx-1), 2:(ny-1)]
    hi = matrix(TRUE, nx-2, ny-2)
    for (dx in -1:1) for (dy in -1:1) {
      if (dx == 0 && dy == 0) next
      hi = hi & (core >= z[(2:(nx-1)) + dx, (2:(ny-1)) + dy]) }
    idx = which(hi, arr.ind = TRUE)
    if (!nrow(idx)) return(fallback)
    vals = core[cbind(idx[,1], idx[,2])]
    o = order(-vals); idx = idx[o, , drop = FALSE]; vals = vals[o]
    px = k$x[idx[,1] + 1L]; py = k$y[idx[,2] + 1L]
    keep = vals >= PEAK_FRAC * vals[1]
    px = px[keep]; py = py[keep]
    sel = 1L
    if (length(px) > 1L) for (i in 2:length(px)) {
      d = sqrt(((px[i] - px[sel])/diff(xlim))^2 + ((py[i] - py[sel])/diff(ylim))^2)
      if (all(d > PEAK_SEP)) sel = c(sel, i) }
    px = px[sel]; py = py[sel]
    dm = outer(xx, px, function(a, b) ((a - b)/diff(xlim))^2) +
         outer(yy, py, function(a, b) ((a - b)/diff(ylim))^2)
    who = max.col(-dm, ties.method = "first")
    mass = tabulate(who, nbins = length(px)) / length(who)
    ok = mass >= MASS_FRAC
    if (!any(ok)) ok[which.max(mass)] = TRUE
    px = px[ok]; py = py[ok]; mass = mass[ok]
    o = order(-mass); px = px[o]; py = py[o]; mass = mass[o]
    keep_n = min(length(px), max(1L, MAX_LOBES))
    list(main = c(px[1], py[1]), mass = mass, drawn = keep_n, all = cbind(px, py),
         n_cand = length(sel),
         sec = if (keep_n > 1L) cbind(px[2:keep_n], py[2:keep_n]) else NULL) }

  PK_ob = kde_peaks(ob$Eh, ob[[yv]]); PK_hb = kde_peaks(bhb$Eh, bhb[[yv]])
  PK_lb = kde_peaks(blb$Eh, blb[[yv]])
  PEAKS[[letter]] <<- list(OFF = PK_ob, ON_high = PK_hb, ON_low = PK_lb, yv = yv,
                           xlim = xlim, ylim = ylim)
  GROUPS[[letter]] <<- list(OFF = ob, ON_high = bhb, ON_low = blb, yv = yv)
  mx = function(pk) pk$main[1]; my = function(pk) pk$main[2]
  ## arrows from the main OFF peak to each main ON peak, trimmed (in inches) to
  ## the diamond's edge; omitted when shorter than 2.4 x the trim
  MDI_CEX = 2.2
  .mr = 0.375 * sqrt(2) * MDI_CEX * par("cex") * par("cin")[1]
  arrow_to = function(x0, y0, x1, y1) {
    if (!all(is.finite(c(x0, y0, x1, y1)))) return(invisible())
    xr = diff(par("usr")[1:2]) / par("pin")[1]
    yr = diff(par("usr")[3:4]) / par("pin")[2]
    dx = (x1 - x0) / xr; dy = (y1 - y0) / yr
    L  = sqrt(dx^2 + dy^2)                         # arrow length, inches
    if (!is.finite(L) || L <= 0) return(invisible())
    back = .mr / (abs(dx) + abs(dy)) * L + 0.008   # diamond edge along this direction
    if (L <= 2.4 * back) return(invisible())
    f0 = back / L; f1 = 1 - back / L
    arrows(x0 + dx * f0 * xr, y0 + dy * f0 * yr,
           x0 + dx * f1 * xr, y0 + dy * f1 * yr,
           col = "grey15", lwd = 2.4, length = 0.13, angle = 22)
  }
  ## drawing order: secondary peaks, arrows, main peaks
  SEC_ALPHA = as.numeric(Sys.getenv("SEC_ALPHA", "0.65"))
  sdi = function(pk, col) if (!is.null(pk$sec))
    points(pk$sec[,1], pk$sec[,2], pch = 23, bg = adjustcolor("white", SEC_ALPHA),
           col = col, cex = MDI_CEX * 0.85, lwd = 2.2)
  sdi(PK_ob, "grey45"); sdi(PK_hb, HIGH); sdi(PK_lb, LOW)
  arrow_to(mx(PK_ob), my(PK_ob), mx(PK_hb), my(PK_hb))   # OFF -> high-burden
  arrow_to(mx(PK_ob), my(PK_ob), mx(PK_lb), my(PK_lb))   # OFF -> low-burden
  mdi = function(pk, col) points(mx(pk), my(pk), pch = 23, bg = col, col = "white",
                                 cex = MDI_CEX, lwd = 1.8)
  mdi(PK_ob, "grey45"); mdi(PK_hb, HIGH); mdi(PK_lb, LOW)
  .pm = function(pk) {
    if (is.null(pk$mass)) return("median")
    lab = sprintf("%.0f%%", 100 * pk$mass)
    if (length(lab) > pk$drawn) lab[(pk$drawn + 1L):length(lab)] =
      paste0("(", lab[(pk$drawn + 1L):length(lab)], ")")
    paste(lab, collapse = "/") }
  cat(sprintf("[%s] %s | centre = %s | peak mass  OFF %s | ON higher %s | ON lower %s\n",
              letter, ttl, CENTRE, .pm(PK_ob), .pm(PK_hb), .pm(PK_lb)))
  box(); title(main = ttl, adj = 0, font.main = 1, line = 0.35,
               cex.main = altv(0.88 * FS, fit_cex(ttl, 1.18 * FS, fit_avail(TRUE))))
  if (!ALT_LABELS)
    mtext(sprintf("n = %s", format(nrow(blb), big.mark = ",")), side = 3, line = 0.35,
          adj = 1, cex = 0.72 * FS, col = "grey25")
  if (showleg) {
    .leg = c("OFF attractor", "ON high-burden attractor", "ON low-burden attractor")
    .lcx = altv(0.5 * FS,
                fit_cex(.leg[2], 0.78 * FS,
                        0.92 * par("pin")[1] - 2.4 * strwidth("M", cex = 0.78 * FS, units = "inches")))
    legend(best_corner(XA, YA, xlim, ylim), bty = "o", bg = "white", box.col = "grey80", box.lwd = 0.8,
           cex = .lcx, pch = 23, col = "white", pt.lwd = 1.2, pt.cex = altv(1.0, 1.25),
           pt.bg = c("grey45", HIGH, LOW), legend = .leg)
  }
  if (nchar(letter)) mtext(letter, side = 3, line = SHIFT_LET_LINE, adj = 0, at = par("usr")[1], font = 2, cex = 1.15 * FS)
}

## ============================================================================
## render
## ============================================================================
cat(sprintf("[Fig. 2] GAP_CUT = %g\n", GAP_CUT))
dir.create("./figures/main", showWarnings = FALSE, recursive = TRUE); dir.create("./source_data", showWarnings = FALSE)
render_fig("./figures/main/Fig_2.png")

## ============================================================================
## numbers behind the figure -> ./source_data/Fig_2_*.csv
## (after rendering, and using no random numbers, so it cannot alter the figure)
## ============================================================================
NUM = list()
put = function(k, v) NUM[[k]] <<- v
q3 = function(x) as.numeric(quantile(x, c(.25, .5, .75), names = FALSE, na.rm = TRUE))

## -- integration bookkeeping
put("rows_read", N_ROWS_READ)
put("rows_expected", 10500L * 8L * 21L)
put("rows_missing", 10500L * 8L * 21L - N_ROWS_READ)
put("rows_nonfinite", N_ROWS_READ - nrow(KD))
put("rows_used", nrow(KD))
put("units", uniqueN(KD[, .(row_ind, beta_ros, r_L_P)]))
put("units_sterile", uniqueN(KD[pathogenic == FALSE, .(row_ind, beta_ros, r_L_P)]))
put("units_pathogenic", uniqueN(KD[pathogenic == TRUE, .(row_ind, beta_ros, r_L_P)]))
put("kinetic_vectors_distinct", uniqueN(KD$row_ind))

## -- lines and the gap criterion
LN_on = KD[eff_ko == 1L, .(n = .N, gap = max(comp) - min(comp)), by = .(row_ind, beta_ros, r_L_P, pathogenic, Kg)]
LN_ko = KD[eff_ko == 0L, .(n = .N, gap = max(comp) - min(comp),
                           rEh = max(Eh) - min(Eh)), by = .(row_ind, beta_ros, r_L_P, pathogenic)]
for (p in c(FALSE, TRUE)) {
  b = if (p) "pathogenic" else "sterile"
  put(sprintf("systems_%s", b), nrow(LN_on[pathogenic == p]))
  put(sprintf("systems_%s_lines_lt21", b), nrow(LN_on[pathogenic == p & n < 21]))
  put(sprintf("systems_%s_gap_gt_cut", b), nrow(LN_on[pathogenic == p & gap > GAP_CUT]))
  put(sprintf("units_%s_with_gap_any_Kg", b),
      uniqueN(LN_on[pathogenic == p & gap > GAP_CUT, .(row_ind, beta_ros, r_L_P)]))
  put(sprintf("ko_lines_%s", b), nrow(LN_ko[pathogenic == p]))
  put(sprintf("ko_lines_%s_lt21", b), nrow(LN_ko[pathogenic == p & n < 21]))
  put(sprintf("ko_lines_%s_gap_gt_cut", b), nrow(LN_ko[pathogenic == p & gap > GAP_CUT]))
  put(sprintf("ko_lines_%s_gap_max", b), max(LN_ko[pathogenic == p, gap]))
  for (g in sort(unique(LN_on$Kg)))
    put(sprintf("frac_%s_gap_gt_cut_Kg%02d", b, g), mean(LN_on[pathogenic == p & Kg == g, gap] > GAP_CUT))
}
put("units_qualifying_total", nrow(QSET))
put("example_candidates_sterile", nrow(.cand))
put("rows_KDq", nrow(KDq))

## -- C/D/E populations; the knockout outcome is a mean over a unit's knockout
##    final states -- how far apart those states are, for the units drawn
for (L in names(GROUPS)) {
  G = GROUPS[[L]]
  sk_used = unique(G$ON_low$sk)
  put(sprintf("panel%s_n_systems", L), nrow(G$ON_low))
  put(sprintf("panel%s_n_units", L), length(sk_used))
  put(sprintf("panel%s_off_missing", L), sum(is.na(G$OFF$Eh)))
  ko_r = LN_ko[paste(row_ind, beta_ros, r_L_P, sep = "|") %in% sk_used]
  put(sprintf("panel%s_ko_range_Eh_max", L), max(ko_r$rEh))
  put(sprintf("panel%s_ko_gap_comp_max", L), max(ko_r$gap))
}
.ks = .gap[hi_lo_gap > GAP_CUT & pathogenic == TRUE]
put("pathogenic_kept_lower_higherEh_and_noHigherPLP", mean(.ks$Eh_l > .ks$Eh_h & .ks$PLP_l <= .ks$PLP_h))
put("pathogenic_kept_lower_higherEh", mean(.ks$Eh_l > .ks$Eh_h))
put("pathogenic_kept_lower_lowerPLP", mean(.ks$PLP_l < .ks$PLP_h))
.ks0 = .gap[hi_lo_gap > GAP_CUT & pathogenic == FALSE]
put("sterile_kept_lower_higherEh", mean(.ks0$Eh_l > .ks0$Eh_h))
put("gap_median_sterile_kept", median(.ks0$hi_lo_gap))
put("gap_median_pathogenic_kept", median(.ks$hi_lo_gap))

CDE = rbindlist(lapply(names(GROUPS), function(L) {
  G = GROUPS[[L]]; P = PEAKS[[L]]; yv = G$yv
  rbindlist(lapply(c("OFF", "ON_high", "ON_low"), function(gn) {
    d = G[[gn]]; pk = P[[gn]]
    qe = q3(d$Eh); qy = q3(d[[yv]])
    allpk = if (!is.null(pk$all)) pk$all else matrix(pk$main, 1)
    ms = if (!is.null(pk$mass)) pk$mass else NA_real_
    data.table(panel = L, y = yv, group = gn, n = nrow(d),
               Eh_q25 = qe[1], Eh_med = qe[2], Eh_q75 = qe[3],
               y_q25 = qy[1], y_med = qy[2], y_q75 = qy[3],
               peak1_Eh = allpk[1, 1], peak1_y = allpk[1, 2], peak1_mass = ms[1],
               peak2_Eh = if (nrow(allpk) > 1) allpk[2, 1] else NA_real_,
               peak2_y  = if (nrow(allpk) > 1) allpk[2, 2] else NA_real_,
               peak2_mass = if (length(ms) > 1) ms[2] else NA_real_,
               n_peak_candidates = if (!is.null(pk$n_cand)) pk$n_cand else NA_integer_,
               n_peaks_kept = nrow(allpk), n_peaks_drawn = if (!is.null(pk$drawn)) pk$drawn else 1L,
               xlim_lo = P$xlim[1], xlim_hi = P$xlim[2], ylim_lo = P$ylim[1], ylim_hi = P$ylim[2])
  }))
}))
fwrite(CDE, "./source_data/Fig_2_panelCDE.csv")

## -- F/G
put("FG_units", uniqueN(pPK$.set))
put("FG_units_production", uniqueN(pPK[src == "production"]$.set))
put("FG_units_extra", N_FG_EXTRA_UNITS)
put("FG_units_extra_drawn", N_FG_EXTRA_DRAWN)
put("FG_pairs", nrow(pPK))
put("FG_pairs_per_Kg_min", min(pPK[, .N, by = Kg]$N))
put("FG_pairs_per_Kg_max", max(pPK[, .N, by = Kg]$N))
put("FG_frac_ties", mean(pPK$dPLP == 0))
put("FG_frac_positive", mean(pPK$dPLP > 0))
put("FG_frac_negative", mean(pPK$dPLP < 0))
put("FG_median_all", median(pPK$dPLP))
put("FG_intact_rows_unpaired", nrow(KDq[eff_ko == 1L & pathogenic == TRUE]) - nrow(pPK[src == "production"]))
put("FG_stat_balanced", as.integer(identical(FG_STAT, "balanced")))
put("FG_stat_subsample", as.integer(identical(FG_STAT, "subsample")))
put("FG_subsample_k", FG_K)
FG = rbind(cB[, .(panel = "F", driver = "beta_ros", level = beta_ros, n_units = n_sets, n_min_stratum, med, lo, hi)],
           cR[, .(panel = "G", driver = "r_L_P",    level = r_L_P,    n_units = n_sets, n_min_stratum, med, lo, hi)])
FG[, n_pairs := c(pPK[, .N, by = beta_ros][order(beta_ros)]$N, pPK[, .N, by = r_L_P][order(r_L_P)]$N)]
## the per-Kg medians the statistic is the median of, as a range
FG[, perKg_min := c(mB[, min(m), by = beta_ros][order(beta_ros)]$V1, mR[, min(m), by = r_L_P][order(r_L_P)]$V1)]
FG[, perKg_max := c(mB[, max(m), by = beta_ros][order(beta_ros)]$V1, mR[, max(m), by = r_L_P][order(r_L_P)]$V1)]
fwrite(FG, "./source_data/Fig_2_panelFG.csv")

## -- A/B example: the stored t = 12000 line of the same system, and the
##    t = 4000 re-integration drawn in A/B
.exbr = beta_ros; .exkg = Kg
.exd = KD[row_ind == ROW_IND & abs(beta_ros - .exbr) < 1e-12 & r_L_P == 0]
put("ex_batch", .exb); put("ex_row_ind", .exr); put("ex_uid", .exd$uid[1])
put("ex_Kg", .exkg); put("ex_beta_ros", .exbr)
put("ex_spread_comp", .cand$sep[1])
.on12 = .exd[arm == sprintf("Kg%02d", .exkg)]; .ko12 = .exd[arm == "KO"]
put("ex_t12000_on_n", nrow(.on12)); put("ex_t12000_ko_n", nrow(.ko12))
put("ex_t12000_on_Eh_min", min(.on12$Eh)); put("ex_t12000_on_Eh_max", max(.on12$Eh))
put("ex_t12000_on_M1_at_Ehmin", .on12[which.min(Eh), M1]); put("ex_t12000_on_M1_at_Ehmax", .on12[which.max(Eh), M1])
put("ex_t12000_ko_Eh_min", min(.ko12$Eh)); put("ex_t12000_ko_Eh_max", max(.ko12$Eh))
put("ex_t12000_ko_M1_mean", mean(.ko12$M1))
put("ex_t12000_on_first_lower_Eh0", min(.on12[Eh > THR, Eh0]))
put("ex_t12000_on_last_higher_Eh0", max(.on12[Eh <= THR, Eh0]))
put("ex_t4000_THR", THR)
put("ex_t4000_ends_Eh_min", min(.ends)); put("ex_t4000_ends_Eh_max", max(.ends))
put("ex_t4000_loE_Eh", loE[1]); put("ex_t4000_loE_M1", loE[2])
put("ex_t4000_hiE_Eh", hiE[1]); put("ex_t4000_hiE_M1", hiE[2])
put("ex_t4000_ofE_Eh", ofE[1]); put("ex_t4000_ofE_M1", ofE[2])
put("ex_t4000_A_on_n_lower", sum(sapply(trajON, function(z) tail(z$Eh, 1)) > THR))
put("ex_t4000_A_off_n_lower", sum(sapply(trajOFF, function(z) tail(z$Eh, 1)) > THR))
put("ex_B_on_n_lower", sum(onV == 1)); put("ex_B_on_n_higher", sum(onV == 2))
put("ex_B_off_n_lower", sum(offV == 3)); put("ex_B_off_n_higher", sum(offV == 4))
put("ex_B_on_boundary", if (any(onV == 1) && !all(onV == 1)) {
  .i0 = which(onV == 1)[1]; (EH0_MAP[.i0 - 1] + EH0_MAP[.i0]) / 2 } else NA_real_)
put("ex_B_on_monotone", all(diff(onV) <= 0))
## how settled the t = 4000 endpoints of A are: |Eh(4000) - Eh(3900)| / 100
.slope = function(z) { n = length(z$Eh); abs(z$Eh[n] - z$Eh[n - 50]) / 100 }
put("ex_t4000_A_max_dEh_dt_last100", max(sapply(c(trajON, trajOFF), .slope)))
for (nm in KIN) put(sprintf("ex_par_%s", nm), par0[[nm]])

fwrite(data.table(quantity = names(NUM), value = vapply(NUM, function(v) format(v, digits = 10), "")),
       "./source_data/Fig_2_numbers.csv")
cat("[Fig. 2] wrote ./source_data/Fig_2_numbers.csv, _panelCDE.csv, _panelFG.csv\n")
