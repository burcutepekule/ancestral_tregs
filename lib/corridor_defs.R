## ============================================================================
## lib/corridor_defs.R -- the machinery behind the corridor panels.
##
## Sourced by 08_corridor_composite.R (one coordinate, the main figure) and by
## 09_corridor_supp.R (all five, one per panel). It exists so those two cannot
## drift: the boundary sweep, the basin rule and the panel geometry are defined
## once. BRANCH must be set before sourcing.
##
## THE BOUNDARY IS RE-BISECTED AT EACH SAMPLED TIME, from the state the
## trajectory has reached -- "given everything else as it is now, how far along
## this coordinate is the flip". It is not the t = 0 wall carried forward by the
## flow; the walls move as the rest of the state moves, and that movement is what
## the panel shows.
## ============================================================================
suppressMessages({library(data.table); library(deSolve); library(parallel)
                  library(ggplot2); library(patchwork); library(scales); library(grid)})
if (!exists("BRANCH")) BRANCH = Sys.getenv("BRANCH", "sterile")
source("./lib/skeleton_defs.R")

CK        = c("tau_id","Kg","beta_ros","r_L_P")
TMAX      = as.numeric(Sys.getenv("TMAX", "12000"))
HMAX      = as.numeric(Sys.getenv("HMAX", "5"))
## SAME ATTRACTOR RULE AS THE SWEEP AND THE t = 0 CORRIDOR: Chebyshev on (Eh, PLP)
## at RAY_CLUST_TOL = 1 (lib/attractor_cluster.R). This file used to break clusters
## on a gap > 2 in the sorted composite, the rule pipeline/01_sweep.R and
## pipeline/21_corridor_tzero.R abandoned because the composite merges states far apart in
## P_LP. The one-combination panels must label basins the way the population does.
RAY_CLUST_TOL = as.numeric(Sys.getenv("RAY_CLUST_TOL", "1.0"))
source("./lib/attractor_cluster.R")
PROBE_BY  = as.numeric(Sys.getenv("PROBE_BY", "10"))
EH0_LO    = as.numeric(Sys.getenv("EH0_LO", "0")); EH0_HI = as.numeric(Sys.getenv("EH0_HI", "99"))
PROBE     = unique(c(seq(EH0_LO, EH0_HI, PROBE_BY), EH0_HI))
N_CORES   = as.integer(Sys.getenv("N_CORES", "1"))
T_END     = as.numeric(Sys.getenv("T_END", "350"))
NT_D      = as.integer(Sys.getenv("NT_D", "1200"))
NSAMP     = as.integer(Sys.getenv("NSAMP", "40"))
NSCAN_D   = as.integer(Sys.getenv("NSCAN_D", "40"))
BIS_D     = as.integer(Sys.getenv("BIS_D", "20"))
EXPAND_MAX= as.integer(Sys.getenv("EXPAND_MAX", "8"))
MARGIN    = as.numeric(Sys.getenv("MARGIN", "0.5"))
## sigma sets where the x axis stops being linear. Small enough that the OFF wall
## of a zero-based coordinate (x_off ~ 0.1 to 1) lands in the log part.
XSIGMA    = as.numeric(Sys.getenv("XSIGMA", "0.1"))
BASE_SIZE = as.numeric(Sys.getenv("BASE_SIZE", "13"))
HARD_D = list(CLP=c(0,300), Eh=c(0,100), Ea=c(0,600), En=c(0,800), M1=c(0,M_TOT),
              R=c(0,50), D=c(0,50), PLP=c(0,300))
COL_SEED  = c(below = "#2a78d6", corridor = "#c2189a", above = "#eb6834")
COL_BASIN = c(resolved = "#2a78d6", chronic = "#eb6834")
ALPHA_FILL = 0.10
SEP_LWD    = 1.1      # separatrix lines in panel A
TRD_LWD    = 1.3      # trajectories in panel A, as in panel B
SEP_PT     = 2.6      # markers on S^ON in panel A
SEP_N_MARK = 10L      # markers on S^ON, at equal distances along the drawn curve
SEP_ASPECT = 1.4      # panel A plot width / height, for distances along the curve
LT  = c(`Tregs-OFF` = "dashed", `Tregs-ON` = "solid")
SYM = list(En = "E[n]", Ea = "E[a]", M1 = "M[1]", R = "R", D = "D", CLP = "C[LP]", Eh = "E[h]")
symof = function(w) if (w %in% names(SYM)) SYM[[w]] else w

## Legend row 3 of panel A (Separatrix: S^OFF dotted, S^ON dotted with markers).
## Panel B draws the same row invisibly, so the Trajectory and Tregs OFF / ON rows
## of the two panels sit at the same heights.
sep_scales = function(visible = TRUE) {
  th = if (visible) NULL else theme(legend.text = element_text(colour = "transparent"),
                                    legend.title = element_text(colour = "transparent"))
  ov_l = list(linewidth = SEP_LWD, alpha = if (visible) 1 else 0)
  ov_p = list(size = SEP_PT, alpha = if (visible) 1 else 0)
  list(scale_linetype_manual(values = c(S_OFF = "dotted", S_ON = "dotted"), name = "Separatrix",
         labels = parse(text = c('S^{OFF}', 'S^{ON}')), drop = FALSE,
         guide = guide_legend(order = 3, override.aes = ov_l, theme = th)),
       scale_shape_manual(values = c(S_OFF = NA, S_ON = 16), name = "Separatrix",
         labels = parse(text = c('S^{OFF}', 'S^{ON}')), drop = FALSE,
         guide = guide_legend(order = 3, override.aes = ov_p, theme = th)))
}

TH_C = function(bs = BASE_SIZE) theme_minimal(base_size = bs) +
  theme(legend.position = "top", legend.box = "vertical",
        legend.key.height = unit(10, "pt"), legend.key.width = unit(bs * 1.4, "pt"), legend.key.spacing.x = unit(3, "pt"),
        legend.margin = margin(0,0,0,0), legend.spacing.y = unit(0, "pt"), legend.box.spacing = unit(4, "pt"),
        legend.text = element_text(size = bs * 0.85, margin = margin(l = 3, r = 8)),
        legend.title = element_text(size = bs * 0.85),
        axis.title = element_text(size = bs * 1.0),
        axis.text  = element_text(size = bs * 0.85, colour = "grey20"),
        plot.tag   = element_text(size = bs * 1.4, face = "plain"))

load_corridor = function(branch = BRANCH) {
  f = sprintf("./data/corridor_q_%s%s.rds", branch, Sys.getenv("CORR_TAG", ""))  # see pipeline/22_corridor_attractors.R
  if (!file.exists(f)) stop("no ", f, " -- run pipeline/23_corridor_position.R first")
  as.data.table(readRDS(f)) }

## ---- choosing the worked example -------------------------------------------
## "widest" gives the largest corridor in WIN, which is what makes the panel
## legible and is by construction an outlier -- the percentile is printed and
## belongs in the caption. "gap" takes the largest rescue gap on the Eh0 line,
## "median" the median gap within the top decile.
##
## "gap0" IS "gap" WITH THE PANEL'S OWN REQUIREMENT ADDED. The corridor band is
## the interval between the two walls, so it can only be drawn at times where
## BOTH exist. A combination whose OFF wall is undefined at t = 0 therefore gives
## a panel whose band begins partway down -- measured on the pathogenic En panel,
## the OFF wall was NA for the first 14 sampled times and the band started at
## t = 2.9. That is a property of the example, not of the branch, so the fix is
## to choose an example that has both walls at t = 0 rather than to explain the
## missing stretch in the caption.
##
## COST. Each candidate needs its own engine (the fixed-point sets for both
## conditions) plus two wall searches, i.e. of order a hundred integrations to
## TMAX. Candidates are tried in descending rescue gap and the first that passes
## wins, so the usual cost is one or two. PICK_TRY caps it.
PICK_TRY = as.integer(Sys.getenv("PICK_TRY", "12"))
pick_combo = function(L, win, rule = "widest", quiet = FALSE) {
  Cd = L[coord == win & q_kind == "measured" & rescue_mode == "bistable" &
         !gap_censored & is.finite(rescue_gap)]
  if (!nrow(Cd)) stop("no usable combination for ", win)
  Cd[, corr_w := x_on - x_off]
  ## Both walls at t = 0, from the same seed the panel sweeps from, using the
  ## same find_wall the panel uses.
  walls_at_zero = function(cand) {
    eng = tryCatch(combo_engine(cand), error = function(e) NULL)
    if (is.null(eng)) return(c(off = NA_real_, on = NA_real_))
    y0 = eng$SETC(eng$ybase, win, make_seeds(cand)[[BND_SEED]])
    c(off = unname(find_wall(eng, win, y0, 0L)[["b"]]),
      on  = unname(find_wall(eng, win, y0, 1L)[["b"]])) }
  pk = switch(rule,
    widest = Cd[which.max(corr_w)],
    gap    = Cd[which.max(rescue_gap)],
    gap0   = {
      ord = Cd[order(-rescue_gap)]
      ntry = min(nrow(ord), PICK_TRY); hit = NULL
      for (i in seq_len(ntry)) {
        cand = ord[i]; w0 = walls_at_zero(cand)
        if (!quiet)
          cat(sprintf("[corridor] %s: candidate %d/%d gap %.3f -> wall at t=0 OFF %s ON %s%s\n",
              win, i, ntry, cand$rescue_gap,
              if (is.finite(w0[["off"]])) sprintf("%.4g", w0[["off"]]) else "NA",
              if (is.finite(w0[["on"]]))  sprintf("%.4g", w0[["on"]])  else "NA",
              if (all(is.finite(w0))) "  <- taken" else ""))
        if (all(is.finite(w0))) { hit = cand; break } }
      if (is.null(hit))
        stop("no combination in the top ", ntry, " by rescue gap has both walls at t = 0 in ",
             win, " -- raise PICK_TRY, or use EX_PICK=gap and say in the caption where the band starts")
      hit },
    median = { hi = Cd[rescue_gap >= quantile(rescue_gap, 0.90, na.rm = TRUE)]
               hi[order(rescue_gap)][ceiling(.N/2)] },
    stop("rule must be widest, gap, gap0 or median"))
  if (!quiet) {
    cat(sprintf("[corridor] %s | tau_id %d Kg %g beta_ros %g r_L_P %g | gap %.3f\n",
                win, pk$tau_id, pk$Kg, pk$beta_ros, pk$r_L_P, pk$rescue_gap))
    cat(sprintf("[corridor] %s: base %.5g x_off %.5g x_on %.5g | width %.4g = %.2f percentile of %s\n",
                win, pk$base, pk$x_off, pk$x_on, pk$corr_w,
                100 * mean(Cd$corr_w <= pk$corr_w, na.rm = TRUE),
                format(nrow(Cd), big.mark = ",")))
  }
  pk }

## ---- one combination's dynamics --------------------------------------------
combo_engine = function(pick) {
  trow = tg[.(pick$tau_id %/% 10000L, pick$tau_id %% 10000L), nomatch = NULL]
  if (!nrow(trow)) stop("no tau grid row for tau_id ", pick$tau_id)
  trow = trow[1]
  par = c(STRUCTURAL, EFFERO); for (nm in KIN) par[[nm]] = trow[[nm]]
  par$Kg = pick$Kg; par$beta_ros = pick$beta_ros
  par$r_L_P = pick$r_L_P; par$r_LP_P = pick$r_L_P
  ybase = seed(pick$e0_base)
  ## M1 is the one coordinate with a conservation constraint: M0 + M1 is fixed,
  ## so setting M1 alone would inflate the pool instead of displacing along it.
  SETC = function(y, v, x) {
    if (v == "M1") { y["M1"] = x; y["M0"] = max(M_TOT - x, 0) } else y[v] = x
    y }
  ## ode_ext: lsodar + the extinction root + the sweep's atol/rtol, trimmed to
  ## `times` (lib/extinction_floor.R), so tail_state's final-time check below is
  ## unaffected by the root row lsodar would otherwise append.
  run = function(y, th, times) tryCatch(
    ode_ext(y, times, rhs_app, c(par, list(theta_treg = th)), maxsteps = 1e5),
    error = function(e) NULL)
  tail_state = function(y, th) {
    o = run(y, th, c(0, TMAX)); if (is.null(o) || nrow(o) < 2) return(NULL)
    if (o[nrow(o), 1] < TMAX * (1 - 1e-9)) return(NULL)
    setNames(as.numeric(o[nrow(o), STATE_ORDER]), STATE_ORDER) }
  FPV = c("Eh","Ea","En","PLP","CLP","M1","R","D","Tr","S","PL")
  snap_fp = function(e, fps) {
    if (is.null(fps) || !nrow(fps) || is.na(e["Eh"])) return(NA_integer_)
    vs = intersect(FPV, names(fps)); ev = e[vs]
    dd = vapply(seq_len(nrow(fps)), function(i) {
      fv = unlist(fps[i, vs]); sum(((ev - fv)/(abs(fv) + 1))^2, na.rm = TRUE) }, numeric(1))
    if (!length(dd) || all(!is.finite(dd))) return(NA_integer_)
    which.min(dd) }
  mkseeds = function(Z) { Z = Z[!vapply(Z, is.null, logical(1))]
    if (!length(Z)) return(NULL)
    m = do.call(rbind, lapply(Z, function(z) c(M1 = unname(z["M1"]), PLP = unname(z["PLP"]),
                                               L = unname(z["Ea"] + z["En"]))))
    m = m[is.finite(rowSums(m)), , drop = FALSE]; if (nrow(m)) unique(m) else NULL }
  fpset = function(th) {
    Z = lapply(PROBE, function(e0) tail_state(seed(e0), th))
    z = tryCatch(as.data.frame(finder(c(par, list(theta_treg = th)), pick$Kg, pick$beta_ros,
                                      pick$r_L_P, th, extra_seeds = mkseeds(Z))),
                 error = function(e) NULL)
    if (is.null(z) || !nrow(z)) stop("no fixed points at theta_treg=", th)
    z$cmp = cmpf(z$Eh, z$PLP); z = z[order(z$cmp), , drop = FALSE]; rownames(z) = NULL
    st = which(z$type == "stable"); if (!length(st)) stop("no stable fixed point at ", th)
    .cl = fp_cluster(z$Eh[st], z$PLP[st], z$cmp[st])
    list(fp = z, lo_ids = st[.cl == 1L]) }
  FS = list(`0` = fpset(0L), `1` = fpset(1L))
  ## Basin by ATTRACTOR IDENTITY -- match the final state to a fixed point and
  ## test membership of the lowest-burden cluster. The same rule 03 and the sweep
  ## use, so a trajectory cannot be labelled one way here and another there.
  lab = function(y, th) {
    z = tail_state(y, th); if (is.null(z)) return(NA_integer_)
    f = FS[[as.character(th)]]; j = snap_fp(z, f$fp)
    if (is.na(j) || !identical(f$fp$type[j], "stable")) return(NA_integer_)
    as.integer(j %in% f$lo_ids) }
  list(par = par, ybase = ybase, SETC = SETC, run = run, lab = lab, FS = FS) }

make_seeds = function(pick) {
  w = pick$x_on - pick$x_off
  s = c(below = pick$x_off - MARGIN * max(w, 1e-9),
        corridor = 0.5 * (pick$x_off + pick$x_on),
        above = pick$x_on + MARGIN * max(w, 1e-9))
  s["below"] = max(s[["below"]], pick$base); s }

run_traj = function(eng, win, seeds, t_end = T_END, nt = NT_D) {
  tt = sort(unique(c(0, exp(seq(log(1e-3), log(t_end), length.out = nt)))))
  TR = rbindlist(lapply(names(seeds), function(nm) rbindlist(lapply(c(0L, 1L), function(th) {
    o = eng$run(eng$SETC(eng$ybase, win, seeds[[nm]]), th, tt); if (is.null(o)) return(NULL)
    d = as.data.table(o); d[, `:=`(seed_key = nm, theta = th)] }), fill = TRUE)), fill = TRUE)
  TR[, seed_key := factor(seed_key, levels = c("below","corridor","above"))]
  TR[, system := factor(ifelse(theta == 0, "Tregs-OFF", "Tregs-ON"),
                        levels = c("Tregs-OFF","Tregs-ON"))]
  TR }

## THE WALL AT ONE STATE, one condition. Lifted out of sweep_bnd so that the
## example picker can ask the same question the panel will ask, with the same
## scan resolution and the same refusal rules, instead of a cheaper proxy that
## could disagree with the figure it is choosing for.
##
## Returns NA when the scan finds no basin change in range, when it finds more
## than one (there is then no single wall, and returning one would invent it), or
## when a label is NA during bisection. The three are not distinguished.
find_wall = function(eng, win, yb, th) {
  rng = HARD_D[[win]]; if (is.null(rng)) rng = c(0, 1e3)
  lb = function(x) eng$lab(eng$SETC(yb, win, x), th)
  hi = min(rng[2], max(as.numeric(yb[[win]]) * 4, 50))
  for (a in seq_len(EXPAND_MAX)) {
    g = seq(rng[1], hi, length.out = NSCAN_D); l = vapply(g, lb, integer(1))
    ok = !is.na(l); g = g[ok]; l = l[ok]
    if (length(l) >= 2 && !all(l == l[1])) {
      tr2 = which(l[-1] != l[-length(l)])
      if (length(tr2) != 1) return(c(b = NA_real_, lo = NA_real_))
      k = tr2[1]; lo2 = g[k]; hi2 = g[k+1]; blo = l[k]
      for (j in seq_len(BIS_D)) { m = (lo2 + hi2)/2; b = lb(m)
        if (is.na(b)) return(c(b = NA_real_, lo = NA_real_))
        if (b == blo) lo2 = m else hi2 = m }
      return(c(b = (lo2 + hi2)/2, lo = as.numeric(blo)))
    }
    if (hi >= rng[2]) break
    hi = min(hi * 2, rng[2])
  }
  c(b = NA_real_, lo = NA_real_) }

sweep_bnd = function(eng, win, pick, seeds, TR, nsamp = NSAMP) {
  bstar = function(yb, th) find_wall(eng, win, yb, th)
  tsamp = sort(unique(c(seq(0, min(5, T_END), length.out = max(10L, nsamp %/% 2L)),
                        exp(seq(log(min(2, T_END)), log(T_END), length.out = nsamp)))))
  jobs = CJ(seed_key = names(seeds), theta = c(0L, 1L))
  one = function(j) {
    nm = as.character(jobs$seed_key[j]); th = jobs$theta[j]
    X = TR[as.character(seed_key) == nm & theta == th][order(time)]
    rbindlist(lapply(tsamp, function(ts) {
      k = which.min(abs(X$time - ts))
      yb = setNames(as.numeric(X[k, ..STATE_ORDER]), STATE_ORDER)
      r = bstar(yb, th)
      data.table(t = X$time[k], b = unname(r[["b"]]), lo_basin = unname(r[["lo"]]),
                 seed_key = nm, theta = th) })) }
  key = paste("v6", pick$tau_id, pick$Kg, pick$beta_ros, pick$r_L_P, win,
              round(seeds[["below"]],5), round(seeds[["corridor"]],5), round(seeds[["above"]],5),
              length(tsamp), T_END, NSCAN_D, BIS_D, TMAX, BRANCH,
              paste0("cheb", RAY_CLUST_TOL), sep = "_")
  cf = file.path(getwd(), "data", paste0("cache_bnd_v6_",
       substr(gsub("[^0-9A-Za-z_.-]", "", key), 1, 150), ".rds"))
  if (file.exists(cf) && !nzchar(Sys.getenv("PANELD_NOCACHE", ""))) {
    B = as.data.table(readRDS(cf))
    cat(sprintf("[corridor] %s: loaded sweep from cache\n", win))
  } else {
    cat(sprintf("[corridor] %s: sweep of %d times x %d jobs on %d cores\n",
                win, length(tsamp), nrow(jobs), N_CORES))
    B = rbindlist(if (N_CORES > 1)
      mclapply(seq_len(nrow(jobs)), function(j) tryCatch(one(j), error = function(e) NULL),
               mc.cores = min(N_CORES, nrow(jobs))) else
      lapply(seq_len(nrow(jobs)), function(j) tryCatch(one(j), error = function(e) NULL)),
      fill = TRUE)
    saveRDS(B, cf)
  }
  if (!nrow(B)) stop("the boundary sweep in ", win, " returned nothing")
  ## AGREEMENT AT t = 0. bstar overwrites the swept coordinate, so at t = 0 every
  ## seed gives the same boundary and it must equal the wall 03 bisected. A
  ## mismatch means the stored corridor and this sweep were produced by different
  ## code -- which is exactly what the M0 = M_TOT - x fix causes for M1 until the
  ## re-run is merged and corridor_q is rebuilt. Reported, loudly, not silently
  ## drawn as though the seeds and the separatrices described the same system.
  t0 = B[is.finite(b), min(t)]
  b0_off = B[t == t0 & theta == 0 & is.finite(b), b][1]
  b0_on  = B[t == t0 & theta == 1 & is.finite(b), b][1]
  chk = function(got, want, nm) {
    if (!is.finite(got) || !is.finite(want)) return(invisible(NULL))
    rel = abs(got - want) / max(abs(want), 1e-9)
    if (rel > 0.02)
      cat(sprintf("[corridor] %s WARNING: %s at t=%.3g is %.5g but 03 stored %.5g (%.0f%% apart)\n",
                  win, nm, t0, got, want, 100 * rel)) }
  chk(b0_off, pick$x_off, "x_off"); chk(b0_on, pick$x_on, "x_on")
  if ((is.finite(b0_off) && abs(b0_off - pick$x_off)/max(abs(pick$x_off),1e-9) > 0.02) ||
      (is.finite(b0_on)  && abs(b0_on  - pick$x_on) /max(abs(pick$x_on), 1e-9) > 0.02))
    cat(sprintf("[corridor] %s: the seeds come from the stored walls and the separatrices from this sweep,\n            so they describe different systems -- rebuild corridor_q before using this panel\n", win))
  B[, seed_key := factor(seed_key, levels = c("below","corridor","above"))]
  B[, system := factor(ifelse(theta == 0, "Tregs-OFF", "Tregs-ON"),
                       levels = c("Tregs-OFF","Tregs-ON"))]
  attr(B, "tsamp") = tsamp
  B }

## ---- the corridor panel -----------------------------------------------------
## WHICH SEED THE SEPARATRICES ARE SWEPT FROM.
##
## The corridor seed loses under Tregs OFF by construction -- that is what makes
## it the corridor seed -- and once it leaves t = 0 there is generally no
## displacement along the coordinate that rescues it, so its OFF boundary ceases
## to exist. Measured on Ea: OFF is 1.327 at t = 0, matching 03 exactly, and NA
## at every later sampled time. A single point cannot be drawn as a path, so the
## OFF separatrix and the corridor band vanish from the panel.
##
## En is the exception. Reducing necrotic burden undoes the thing driving the
## failure, so its OFF wall survives to about t = 30 and the panel has a curve.
##
## The "below" seed resolves under BOTH conditions, so both boundaries stay
## defined for as long as the trajectory does and every coordinate gets two
## curves. That is the default here; BND_SEED=corridor reproduces v5's choice,
## which is the right one when the point is the corridor seed's own fate.
BND_SEED = Sys.getenv("BND_SEED", "below")
panel_corridor = function(win, TR, BND, bs = BASE_SIZE, show_legend = TRUE,
                          bnd_seed = BND_SEED) {
  tsamp = attr(BND, "tsamp")
  TRD = TR[time >= min(tsamp) & time <= T_END]
  BR  = BND[is.finite(b)][order(seed_key, theta, t)]
  BRc = BR[seed_key == bnd_seed]
  if (!nrow(BRc)) stop("no boundary located for the ", bnd_seed, " seed in ", win)
  nth = BRc[, uniqueN(theta)]
  cat(sprintf("[corridor] %s: separatrices from the '%s' seed | OFF %d pts, ON %d pts\n",
              win, bnd_seed, BRc[theta == 0, .N], BRc[theta == 1, .N]))
  if (nth < 2)
    cat(sprintf("[corridor] %s: only one condition has a boundary -- no corridor band\n", win))
  XLIM = c(0, max(c(BR$b, TRD[[win]]), na.rm = TRUE) * 1.02)
  CORR = BRc[, {
    off = b[theta == 0]; on = b[theta == 1]
    if (!length(off) || !length(on) || is.na(off[1]) || is.na(on[1])) {
      .(lo = NA_real_, hi = NA_real_)
    } else {
      .(lo = min(off[1], on[1]), hi = max(off[1], on[1]))
    } }, by = t]
  CORR = CORR[is.finite(lo) & is.finite(hi) & hi > lo][order(t)]
  TSc = sort(unique(BRc$t)); CORR[, idx := match(t, TSc)]
  CORR[, run := cumsum(c(1L, as.integer(diff(idx) != 1L)))]
  res_low = isTRUE(mean(BRc$lo_basin == 1, na.rm = TRUE) >= 0.5)
  if (res_low) CORR[, `:=`(RES_LO = XLIM[1], RES_HI = lo, CHR_LO = hi, CHR_HI = XLIM[2])] else
               CORR[, `:=`(RES_LO = hi, RES_HI = XLIM[2], CHR_LO = XLIM[1], CHR_HI = lo)]
  wmed = BRc[is.finite(b), .(mb = median(b)), by = theta]
  chrth = if (res_low) wmed[which.max(mb), theta] else wmed[which.min(mb), theta]
  CHR = BRc[theta == chrth & is.finite(b), .(t, cw = b)][order(t)]
  CHR[, idx := match(t, TSc)]; CHR[, run := cumsum(c(1L, as.integer(diff(idx) != 1L)))]
  if (res_low) CHR[, `:=`(CHR_LO = cw, CHR_HI = XLIM[2])] else
               CHR[, `:=`(CHR_LO = XLIM[1], CHR_HI = cw)]
  sy = symof(win)
  ## THE X SCALE IS CHOSEN PER PANEL, from the walls themselves.
  ##
  ## Neither choice works everywhere. For Ea the walls sit at 1.33 and 46, so a
  ## linear axis collapses the OFF separatrix and the resolving band onto the
  ## left edge and they look absent. For En they sit at 48 and 88 on a panel that
  ## runs to 155, and a pseudo-log axis with a small sigma puts 48 at 84% of the
  ## width -- the band is drawn correctly and reads as though it reached 150.
  ##
  ## So: log only when the walls actually span decades, and place the linear-to-log
  ## knee below the smallest wall rather than at a fixed value, so no wall lands in
  ## the compressed part. XLOG=1 or 0 overrides.
  bb = BRc$b[is.finite(BRc$b) & BRc$b > 0]
  span_ratio = if (length(bb) > 1)
    as.numeric(quantile(bb, 0.95) / max(quantile(bb, 0.05), 1e-12)) else 1
  .xl = Sys.getenv("XLOG", "auto")
  xlog = if (identical(.xl, "1")) TRUE else
         if (identical(.xl, "0")) FALSE else span_ratio > 20
  sig = if (length(bb)) max(min(bb) / 3, 1e-3) else XSIGMA
  cat(sprintf("[corridor] %s: wall span ratio %.1f -> %s x axis%s\n", win, span_ratio,
              if (xlog) "pseudo-log" else "linear",
              if (xlog) sprintf(" (sigma %.3g)", sig) else ""))
  XSCALE = if (xlog) {
    scale_x_continuous(trans = pseudo_log_trans(sigma = sig, base = 10))
  } else {
    scale_x_continuous()
  }
  BRs = copy(BRc)[is.finite(b)]
  BRs[, sep := factor(fifelse(theta == 0L, "S_OFF", "S_ON"), levels = c("S_OFF", "S_ON"))]
  setorder(BRs, theta, t)
  ## markers on S^ON at equal distances ALONG THE DRAWN CURVE: arc length in plot
  ## coordinates (each axis transformed as drawn and scaled to the panel's range),
  ## markers interpolated linearly between samples, i.e. on the drawn segments
  ytr = function(t) asinh(t / 2) / log(10)
  xtr = if (xlog) function(x) asinh(x / (2 * sig)) / log(10) else identity   # pseudo_log_trans
  BRs[, mark := FALSE]
  MK = BRs[sep == "S_ON"][order(t)]
  if (nrow(MK) > 1) {
    xr = diff(range(xtr(XLIM))); yr = diff(range(ytr(c(0, T_END))))
    u = SEP_ASPECT * xtr(MK$b) / xr; v = ytr(MK$t) / yr
    s_arc = c(0, cumsum(sqrt(diff(u)^2 + diff(v)^2)))
    at = seq(0, max(s_arc), length.out = SEP_N_MARK + 2L)[-c(1L, SEP_N_MARK + 2L)]
    yi = approx(s_arc, ytr(MK$t), at)$y
    MKP = data.table(b = approx(s_arc, MK$b, at)$y, t = 2 * sinh(yi * log(10)),
                     sep = factor("S_ON", levels = c("S_OFF", "S_ON")), mark = TRUE)
  } else MKP = BRs[0]
  g = ggplot() +
    geom_ribbon(data = CORR, aes(y = t, xmin = RES_LO, xmax = RES_HI, group = run),
                orientation = "y", fill = COL_BASIN[["resolved"]], alpha = ALPHA_FILL) +
    geom_ribbon(data = CHR, aes(y = t, xmin = CHR_LO, xmax = CHR_HI, group = run),
                orientation = "y", fill = COL_BASIN[["chronic"]], alpha = ALPHA_FILL) +
    geom_ribbon(data = CORR, aes(y = t, xmin = lo, xmax = hi, group = run),
                orientation = "y", fill = "white", alpha = 1) +
    geom_path(data = TRD, aes(x = get(win), y = time, group = interaction(seed_key, theta),
                              colour = seed_key, linetype = system), linewidth = TRD_LWD) +
    scale_colour_manual(values = COL_SEED, name = "Trajectory",
                        labels = c(below = "Below", corridor = "Corridor", above = "Above"),
                        guide = guide_legend(order = 1, override.aes = list(linewidth = 1.3))) +
    ## Rows 1 and 2 of the legend are those of panel B (Trajectory; Tregs OFF / ON).
    scale_linetype_manual(values = LT, name = NULL,
                          guide = guide_legend(order = 2, override.aes = list(linewidth = 1.3, colour = "black"))) +
    ## Row 3: the separatrices, on their own scales so their line types (dotted;
    ## S^ON with markers) cannot be confused with the dashed Tregs-OFF trajectories.
    ## The linetype and shape scales share name and labels, so ggplot merges them
    ## into one legend row. The coordinate is NOT in the key: a supplement draws one
    ## legend for several panels, and the coordinate is on each panel's x axis.
    ggnewscale::new_scale("linetype") +
    geom_path(data = BRs, aes(x = b, y = t, group = interaction(seed_key, theta), linetype = sep),
              colour = "grey15", linewidth = SEP_LWD) +
    geom_point(data = MKP, aes(x = b, y = t, shape = sep),
               colour = "grey15", size = SEP_PT) +
    sep_scales(TRUE) +
    scale_y_continuous(trans = compose_trans(pseudo_log_trans(sigma = 1, base = 10), "reverse"),
                       breaks = { b = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000); b[b <= T_END] }) +
    XSCALE +
    coord_cartesian(xlim = XLIM, ylim = c(T_END, 0), expand = FALSE) +
    labs(x = parse(text = sy), y = "Time") + TH_C(bs)
  if (!show_legend) g = g + theme(legend.position = "none")
  g }

## ---- the phase portrait -----------------------------------------------------
panel_phase = function(eng, win, seeds, TR, bs = BASE_SIZE) {
  ENDS = TR[, .SD[.N], by = .(seed_key, theta, system)]
  STARTS = TR[, .SD[1], by = .(seed_key, theta)]
  ENDS[, cmpv := cmpf(Eh, PLP)]
  ENDS[, basin := {
    b = vapply(seq_len(.N), function(i)
          eng$lab(eng$SETC(eng$ybase, win, seeds[[as.character(seed_key[i])]]), theta[i]),
          integer(1))
    fifelse(is.na(b), "Unassigned", fifelse(b == 1L, "Low-burden attractor", "High-burden attractor")) }]
  ENDS[, txt := sprintf("%s\n(Tregs-%s)", basin, ifelse(system == "Tregs-OFF", "OFF", "ON"))]
  ENDS[, col := COL_SEED[as.character(seed_key)]]
  ## One label per ATTRACTOR, in black: two initial conditions that reach the same
  ## attractor under the same condition end at the same point, and a coloured label
  ## on each read as two attractors. The label sits at the mean of the endpoints it
  ## stands for; the endpoints and paths keep the trajectory colours. The two
  ## low-burden attractors sit at the far right next to each other, so their
  ## labels go to the RIGHT of the points, stacked; the rest are placed by the repel.
  ATT = ENDS[basin != "Unassigned", .(Eh = mean(Eh), M1 = mean(M1), col = "black"), by = .(system, basin, txt)]
  ATT[, right := basin == "Low-burden attractor"]
  ## the high-burden Tregs-ON label is pinned to the right of its point, off the path
  REPH = ATT[right == FALSE & system == "Tregs-ON", .(Eh, M1, txt, col)]
  REP  = rbind(ATT[right == FALSE & system != "Tregs-ON", .(Eh, M1, txt, col)],
               ENDS[basin == "Unassigned", .(Eh, M1, txt, col)])
  REPR = ATT[right == TRUE, .(Eh, M1, txt, col)]
  ## the initial-condition label sits just below its point, clear of the paths
  REPI = data.table(Eh = STARTS$Eh[1], M1 = STARTS$M1[1], txt = "Initial conditions", col = "black")
  lab_geom = if (requireNamespace("ggrepel", quietly = TRUE)) list(
    ## the remaining labels (high-burden, Tregs OFF) sit above-right of their point,
    ## in the empty corner, clear of the path that runs down from it
    geom_segment(data = REP, aes(x = Eh, y = M1, xend = Eh + 3, yend = M1 + 0.6), colour = REP$col,
        linewidth = 0.7, inherit.aes = FALSE, show.legend = FALSE),
    geom_text(data = REP, aes(Eh + 3.5, M1 + 0.6, label = txt), colour = REP$col, hjust = 0,
        vjust = 0, size = bs * 0.23, lineheight = 0.9, inherit.aes = FALSE, show.legend = FALSE),
    ## the right-hand labels are kept to the right of both points (xlim), stacked
    ggrepel::geom_text_repel(data = REPR, aes(Eh, M1, label = txt), colour = REPR$col,
        segment.colour = REPR$col, inherit.aes = FALSE, nudge_x = 9, direction = "y",
        xlim = c(max(REPR$Eh) + 8, NA), hjust = 0, force = 10, box.padding = 0.6,
        point.padding = 0.3, segment.size = 0.7, min.segment.length = 0, max.overlaps = Inf,
        size = bs * 0.23, seed = 1, show.legend = FALSE, lineheight = 0.9),
    ggrepel::geom_text_repel(data = REPH, aes(Eh, M1, label = txt), colour = "black",
        segment.colour = "black", inherit.aes = FALSE, nudge_x = 12, nudge_y = 1, direction = "y",
        xlim = c(max(REPH$Eh) + 8, NA), hjust = 0, segment.size = 0.7, min.segment.length = 0,
        size = bs * 0.23, seed = 1, show.legend = FALSE, lineheight = 0.9),
    ggrepel::geom_text_repel(data = REPI, aes(Eh, M1, label = txt), colour = "black",
        segment.colour = "black", inherit.aes = FALSE, nudge_y = -2.2, direction = "x",
        hjust = 0, segment.size = 0.7, min.segment.length = 0, size = bs * 0.23, seed = 1,
        show.legend = FALSE)) else
    geom_text(data = REP, aes(Eh, M1, label = txt), colour = REP$col, inherit.aes = FALSE,
        size = bs * 0.23, hjust = -0.1, vjust = -0.6, show.legend = FALSE)
  ## Direction: a few arrowheads ON each path, placed by arc length in the plot's
  ## own units (each axis scaled to its range, so a head every so often on the
  ## page rather than in Eh units). No shaft: a segment from the previous sample
  ## to the head's sample carries the head and lies on the path. The six paths
  ## share their first stretch, so Tregs OFF and Tregs ON take different fractions
  ## and the heads do not stack.
  ## One head per path, at HEAD_OFF of normalised arc length (the plot diagonal is
  ## ~1.4), the offset differing between the two conditions. The six paths share
  ## their first stretch, so this draws two heads on it, one per condition; more
  ## heads (HEAD_STEP finite) crowd the turn near the low-burden attractors. The
  ## direction is taken over the last few samples so a wiggle does not turn a head.
  HEAD_STEP = Inf; HEAD_OFF = c(`0` = 0.25, `1` = 0.50)
  rx = diff(range(TR$Eh)); ry = diff(range(TR$M1))
  HEADS = TR[, {
    sl = cumsum(c(0, sqrt((diff(Eh) / rx)^2 + (diff(M1) / ry)^2)))
    at = if (is.finite(HEAD_STEP)) seq(HEAD_OFF[[as.character(theta[1])]], sl[.N] * 0.97, by = HEAD_STEP) else HEAD_OFF[[as.character(theta[1])]]
    ii = vapply(at, function(f) max(4L, which(sl >= f)[1]), integer(1))
    .(x0 = Eh[ii - 3L], y0 = M1[ii - 3L], x1 = Eh[ii], y1 = M1[ii]) }, by = .(seed_key, theta)]
  HEADS[, system := factor(fifelse(theta == 0L, "Tregs-OFF", "Tregs-ON"), levels = names(LT))]
  SEP_DUMMY = data.table(Eh = rep(TR$Eh[1], 4), M1 = rep(TR$M1[1], 4),
                         sep = factor(rep(c("S_OFF", "S_ON"), each = 2), levels = c("S_OFF", "S_ON")))
  list(plot = ggplot(TR, aes(Eh, M1, colour = seed_key, linetype = system,
                             alpha = seed_key, group = interaction(seed_key, theta))) +
    geom_path(linewidth = 1.3) +
    geom_segment(data = HEADS, aes(x = x0, y = y0, xend = x1, yend = y1, colour = seed_key, alpha = seed_key),
                 linewidth = 1.3, linetype = "solid", inherit.aes = FALSE, show.legend = FALSE,
                 arrow = arrow(length = unit(0.34, "cm"), angle = 25, type = "closed", ends = "last")) +
    geom_point(data = STARTS, aes(Eh, M1), shape = 21, fill = "black", colour = "black",
               size = 3.2, stroke = 0.3, inherit.aes = FALSE, show.legend = FALSE) +
    geom_point(data = ENDS, size = 3.2, show.legend = FALSE) +
    lab_geom +
    scale_colour_manual(values = COL_SEED, name = "Trajectory",
                        labels = c(below = "Below", corridor = "Corridor", above = "Above"),
                        guide = guide_legend(order = 1, override.aes = list(linewidth = 1.3, alpha = 1))) +
    scale_alpha_manual(values = c(below = 0.5, corridor = 1, above = 0.5), guide = "none") +
    scale_linetype_manual(values = LT, name = NULL,
                          guide = guide_legend(order = 2, override.aes = list(linewidth = 1.3, colour = "black"))) +
    ## invisible placeholder for panel A's separatrix row (see sep_scales)
    ggnewscale::new_scale("linetype") +
    geom_path(data = SEP_DUMMY, aes(Eh, M1, linetype = sep, group = sep), inherit.aes = FALSE, alpha = 0) +
    geom_point(data = SEP_DUMMY, aes(Eh, M1, shape = sep), inherit.aes = FALSE, alpha = 0) +
    sep_scales(FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.95))) +   # room for the two labels right of the points
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.26))) +
    labs(x = expression(E[h]), y = expression(M[1])) + TH_C(bs),
    ends = ENDS) }
