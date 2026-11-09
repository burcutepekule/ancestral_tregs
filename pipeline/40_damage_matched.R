## ============================================================================
## pipeline/40_damage_matched.R -- the damage-matched experiment behind Fig. 4
## and Extended Data Fig. 5. Everything
## fixed except beta_ros: at the SAME integrated epithelial damage, does the Tregs ON
## system run more aggressively, and what does that do to the pathogen?
## ----------------------------------------------------------------------------
## ONE UNIT = (combination, initial condition E_h(0)): a gated pathogenic combination
## (kinetic row x K_g x r_L_P x beta_r) in which the Tregs OFF system does not
## eliminate the pathogen, drawn by pipeline/43_damage_matched_gated_units.R and read
## from UNIT_FILE (one row per unit: batch, row_ind, tau_id, Kg, r_L_P, b_off, Eh0;
## unit i is row i). The parameter vector, the initial condition and the solver are
## those of the sweep (lib/skeleton_defs.R, lib/extinction_floor.R).
## Per unit, three runs:
##   OFF(beta_r)   Tregs OFF at the unit's own level -- the reference cost
##   ON (beta_r)   Tregs ON at the SAME level -- the brake alone; also the check that the
##                 arms are identical before Tregs act, and the divergence time
##   ON (beta*)    Tregs ON with beta_ros searched so that its integrated damage equals
##                 OFF(beta_r)'s, to a relative TOL
## SKIP_OFF_ELIM=1 (the cluster default): a unit whose Tregs OFF run eliminates the
## pathogen after all is recorded (status off_eliminates) and not searched.
## THE COST THAT IS MATCHED, one search per MATCH mode (the figures use "conv"):
##   conv      integrated damage over [0, max(T_conv OFF, T_conv ON)], where T_conv
##             of a run is the first time after which every state stays within
##             CONV_REL of its final value (or CONV_ABS for states near zero).
##             The window moves with the searched beta, so the match is verified
##             on the final pair's own window.
##   total     integrated damage over the whole run, t = 0 .. TMAX.
##   acute01/05/10  integrated damage over one common window for both arms: from
##             t = 0 until the Tregs OFF reference run's pathogen first falls to
##             1% / 5% / 10% of its peak after the peak (whole run if never).
## THE SEARCH. The cost is not monotone in beta (raising beta kills epithelium through
## ROS but also kills pathogen), so beta* is the LARGEST beta at which the cost crosses
## the target upwards -- above it every beta tested costs Tregs ON more damage than
## Tregs OFF.
##   1. scan: f(beta) = relative damage gap on the log grid beta_r * GRID_STEP^k
##      over the whole [BETA_LO, BETA_HI] (beta_r is a grid point)
##   2. probe: every grid local minimum of f in (0, NEAR) gets N_PROBE extra points
##      between its neighbours, so a dip below the target narrower than a grid
##      step is not missed there
##   3. bracket: the last point with f <= 0 and the next one (f > 0); bisection on
##      log beta towards the crossing, |f| <= TOL_BIS, keeping the closest point;
##      accepted if it meets TOL (else no_match_jump).
##      No point with f <= 0 at all: no_match_lo; f <= 0 at BETA_HI: no_match_hi
##   4. jump test: f at beta* x (0.97, 0.99, 1.01, 1.03); jump_pts = the change of
##      f in percentage points across +-3% (recorded, not filtered here)
##   Recorded too: n_cross_up / n_cross_down on the scanned profile, b_first_up
##   (the first upward crossing at or above beta_r, log-interpolated), and the whole
##   profile in <OUT>_profile.rds (one row per unit x mode x beta). At large beta Tregs
##   ON loses its epithelium and its cost tends to the maximum, K_E x window, which
##   Tregs OFF cannot exceed, so no upward crossing exists past saturation; BETA_HI = 10
##   is past it (E_h at the end ~0.1).
## ARRAY_N / ARRAY_ID: this process runs the units with (i - 1) %% ARRAY_N == ARRAY_ID
## (a Slurm array task, pipeline/cluster/damage_matched_gated.sh); RESUME=1 skips the
## units already in OUT. Every outcome is recorded; no unit is selected or dropped on
## its result.
##
## THE SPLITS: t_div (pathogen) -- the first output time the pathogen of
## ON(beta_r) departs from OFF(beta_r)'s by more than DIV_FRAC; t_div_M1 --
## the same on M1 (the time the figures use). t_arr -- first time Tr exceeds
## ARR_FRAC of its maximum.
##
## RECORDED PER RUN: integrated pathogen and damage over the whole run, before /
## after t_div and before t_div_M1; the pathogen peak and its time; the time to
## extinction; residual pathogen, final epithelium and whether the efferocytosis
## switch is open at the end (Ea + En < Kg); for ON(beta*) the window, the target
## damage, the achieved relative error, and the window pathogen of OFF and ON(beta*).
##
## in : UNIT_FILE (./data/fig04_gated_units.rds), ./data/tau_grid_v6_*.rds, ./lib/*
## out: OUT (one row per unit x MATCH mode), written every FLUSH units, + <OUT>_profile.rds
## RUN: sbatch pipeline/cluster/damage_matched_gated.sh         (the cluster array)
##      UNITS=1,2,3 OUT=./data/test.rds Rscript pipeline/40_damage_matched.R   # a few units
## ============================================================================
suppressMessages({library(data.table); library(parallel); library(deSolve)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
assign("BRANCH", Sys.getenv("BRANCH", "pathogenic_harm"), envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R"))
suppressMessages(source("./lib/extinction_floor.R"))
if (!isTRUE(PL0 > 0)) stop("skeleton_defs.R did not load a pathogenic branch")

BETA_LO  = as.numeric(Sys.getenv("BETA_LO", "0.0005"))
BETA_HI  = as.numeric(Sys.getenv("BETA_HI", "10"))
GRID_STEP = as.numeric(Sys.getenv("GRID_STEP", "1.03")) # scan grid, x on beta
NEAR     = as.numeric(Sys.getenv("NEAR", "0.02"))       # probe local minima of f below this
N_PROBE  = as.integer(Sys.getenv("N_PROBE", "6"))       # extra points per probed minimum
JUMP_X   = c(0.97, 0.99, 1.01, 1.03)                    # the jump test around beta*
ARRAY_N  = as.integer(Sys.getenv("ARRAY_N", "1"))
ARRAY_ID = as.integer(Sys.getenv("ARRAY_ID", "0"))
RESUME   = Sys.getenv("RESUME", "0") == "1"
TOL      = as.numeric(Sys.getenv("TOL", "1e-3"))        # relative damage match
TOL_BIS  = as.numeric(Sys.getenv("TOL_BIS", "1e-4"))    # bisection target, |f| <= TOL_BIS (accept |f| <= TOL)
MAX_BIS  = as.integer(Sys.getenv("MAX_BIS", "30"))
ARR_FRAC = as.numeric(Sys.getenv("ARR_FRAC", "0.01"))   # Treg accumulation: Tr > 1% of its max
DIV_FRAC = as.numeric(Sys.getenv("DIV_FRAC", "0.01"))   # the arms diverge: pathogen differs by > 1%
MATCH    = strsplit(Sys.getenv("MATCH", "conv"), ",")[[1]]
ACUTE    = c(acute01 = 0.01, acute05 = 0.05, acute10 = 0.10)
stopifnot(all(MATCH %in% c("total", "conv", names(ACUTE))))
## conv: the episode ends when BOTH runs of the pair have converged -- every state
## within CONV_REL of its final value (or CONV_ABS for states near zero) from then
## on; both runs are charged over [0, max(T_conv OFF, T_conv ON)]. The window moves
## with the searched beta, so the match is verified on the final pair's own window.
CONV_REL = as.numeric(Sys.getenv("CONV_REL", "0.01"))
CONV_ABS = as.numeric(Sys.getenv("CONV_ABS", "0.01"))
TMAX     = as.numeric(Sys.getenv("TMAX", "12000"))
TT_N     = as.integer(Sys.getenv("TT_N", "400"))
KE       = as.numeric(Sys.getenv("K_E", "100"))
N_CORES  = as.integer(Sys.getenv("N_CORES", "1"))
SEED     = as.integer(Sys.getenv("SEED", "1"))
FLUSH    = as.integer(Sys.getenv("FLUSH", "400"))
OUT      = Sys.getenv("OUT", "./data/fig04_damage_matched_gated.rds")
UNITS    = Sys.getenv("UNITS", "")                      # "" = every unit; else comma list of unit indices i
UNIT_FILE = Sys.getenv("UNIT_FILE", "./data/fig04_gated_units.rds")   # one row per unit (pipeline/43_damage_matched_gated_units.R)
SKIP_OFF_ELIM = Sys.getenv("SKIP_OFF_ELIM", "0") == "1" # 1: a unit whose Tregs-OFF run eliminates the pathogen is recorded
                                                        # (status off_eliminates) and not searched
ELIM_TOL = as.numeric(Sys.getenv("ELIM_TOL", "1e-3"))   # SKIP_OFF_ELIM: eliminated = extinction floor reached, or final P_LP below this
TT = c(0, exp(seq(log(0.05), log(TMAX), length.out = TT_N)))

## ---- the units: the rows of UNIT_FILE (unit i is row i) ----------------------------------
U = as.data.table(readRDS(UNIT_FILE))
need = c("batch", "row_ind", "tau_id", "Kg", "r_L_P", "b_off", "Eh0")
if (!all(need %in% names(U))) stop(UNIT_FILE, " lacks: ", paste(setdiff(need, names(U)), collapse = ", "))
U = U[, ..need]
cat(sprintf("[fig04 exp] %d units from %s (%d cells) | r_L_P mix: %s\n", nrow(U), UNIT_FILE, uniqueN(U[, .(tau_id, Kg, r_L_P)]),
            paste(U[, .N, by = r_L_P][order(r_L_P), sprintf("%.2f:%d", r_L_P, N)], collapse = " ")))

## ---- one integration, everything kept ----------------------------------------------
run = function(p, e0, beta, th) {
  pp = p; pp$beta_ros = beta
  o = tryCatch(ode_ext(seed(e0), TT, rhs_app, c(pp, list(theta_treg = th)), maxsteps = 1e5),
               error = function(e) NULL)
  if (is.null(o) || nrow(o) < length(TT)) return(NULL)
  o = as.data.table(o)
  if (any(!is.finite(o$PLP)) || any(!is.finite(o$Eh))) return(NULL)
  o
}
trap = function(t, y) sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
summ = function(o, t_split, Kg) {
  i  = max(1L, findInterval(t_split, o$time))
  pre = seq_len(i); post = i:nrow(o)
  ip = which.max(o$PLP); iz = which(o$PLP <= EXT_TOL)[1]; n = nrow(o)
  list(auc_P = trap(o$time, o$PLP), auc_D = trap(o$time, KE - o$Eh),
       auc_P_pre = trap(o$time[pre], o$PLP[pre]), auc_P_post = trap(o$time[post], o$PLP[post]),
       auc_D_pre = trap(o$time[pre], KE - o$Eh[pre]), auc_D_post = trap(o$time[post], KE - o$Eh[post]),
       peak_P = o$PLP[ip], t_peak = o$time[ip], t_clear = if (is.na(iz)) NA_real_ else o$time[iz],
       P_end = o$PLP[n], Eh_end = o$Eh[n], open_end = (o$Ea[n] + o$En[n]) < Kg)
}
arrival = function(o) { m = max(o$Tr); if (!is.finite(m) || m <= 0) return(NA_real_)
  o$time[which(o$Tr > ARR_FRAC * m)[1]] }
## the end of the acute phase: first time after the peak the pathogen is at or
## below frac x peak; NA if it never gets there
acute_end = function(o, frac) { ip = which.max(o$PLP)
  j = which(o$PLP[ip:nrow(o)] <= frac * o$PLP[ip])[1]
  if (is.na(j)) NA_real_ else o$time[ip + j - 1L] }
## the window of a MATCH mode, from the Tregs OFF reference run
window = function(off, mode) {
  if (mode == "total") return(c(te = TMAX, reached = 1))
  te = acute_end(off, ACUTE[[mode]]); if (is.finite(te)) c(te = te, reached = 1) else c(te = TMAX, reached = 0) }
## damage of a run over [0, te]
dmg_to = function(o, te) { k = which(o$time <= te); trap(o$time[k], KE - o$Eh[k]) }
## the first time after which every state stays within tolerance of its final value
tconv = function(o) {
  X = as.matrix(o[, ..STATE_ORDER]); fin = X[nrow(X), ]
  bad = sweep(abs(sweep(X, 2, fin)), 2, pmax(CONV_REL * abs(fin), CONV_ABS), ">")
  k = which(rowSums(bad) > 0); k = if (length(k)) max(k) else 0L
  if (k >= nrow(X)) TMAX else o$time[k + 1L] }
## M1 divergence and the integrals before it
div_M1 = function(off, on0) { rel = abs(on0$M1 - off$M1) / pmax(off$M1, 1e-9); k = which(rel > DIV_FRAC)[1]
  if (is.na(k)) TMAX else off$time[k] }
pre_int = function(o, ts) { i = max(1L, findInterval(ts, o$time)); c(P = trap(o$time[1:i], o$PLP[1:i]), D = trap(o$time[1:i], KE - o$Eh[1:i])) }
## every mode's damage over the common windows W, plus the run's own acute ends
all_costs = function(o, W, sfx) {
  v = c(unlist(lapply(names(W), function(m) setNames(dmg_to(o, W[[m]][["te"]]), paste0("D_", m)))),
        unlist(lapply(names(ACUTE), function(m) setNames(acute_end(o, ACUTE[[m]]), paste0("tend_own_", m)))))
  setNames(as.list(v), paste0(names(v), sfx)) }

## ---- one unit ----------------------------------------------------------------------
one = function(i) {
  u = U[i]
  trow = tg[.(u$batch, u$row_ind), nomatch = NULL]
  if (!nrow(trow)) return(data.table(i = i, status = "no_tau_row"))
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = u$Kg; p$r_L_P = u$r_L_P; p$r_LP_P = u$r_L_P
  off = run(p, u$Eh0, u$b_off, 0); on0 = run(p, u$Eh0, u$b_off, 1)
  if (is.null(off) || is.null(on0)) return(data.table(i = i, status = "fail_reference"))
  ## the gate on the Tregs-OFF run itself: nothing to gain where it already eliminates
  if (SKIP_OFF_ELIM && (any(off$PLP <= EXT_TOL) || off$PLP[nrow(off)] < ELIM_TOL))
    return(cbind(data.table(i = i), u[, .(tau_id, batch, row_ind, Kg, r_L_P, Eh0, b_off)],
                 data.table(status = "off_eliminates", t_clear_off = off$time[which(off$PLP <= EXT_TOL)[1]], P_end_off = off$PLP[nrow(off)])))
  ## when the arms diverge at the same level: the split for every run of the unit
  rel = abs(on0$PLP - off$PLP) / pmax(off$PLP, 1e-9)
  iv = which(rel > DIV_FRAC)[1]
  t_div = if (is.na(iv)) TMAX else on0$time[iv]
  S_off = summ(off, t_div, u$Kg); S_on0 = summ(on0, t_div, u$Kg)
  t_div_M1 = div_M1(off, on0); tc_off = tconv(off); tc_on0 = tconv(on0)
  pm_off = pre_int(off, t_div_M1); pm_on0 = pre_int(on0, t_div_M1)
  Wm = setNames(lapply(c("total", names(ACUTE)), function(m) window(off, m)), c("total", names(ACUTE)))
  ref = cbind(data.table(i = i, tau_id = u$tau_id, batch = u$batch, row_ind = u$row_ind, Kg = u$Kg,
                         r_L_P = u$r_L_P, Eh0 = u$Eh0, b_off = u$b_off, t_div = t_div, t_arr_same = arrival(on0),
                         t_div_M1 = t_div_M1, tconv_off = tc_off, tconv_onB = tc_on0,
                         P_preM1_off = pm_off[["P"]], D_preM1_off = pm_off[["D"]], P_preM1_onB = pm_on0[["P"]], D_preM1_onB = pm_on0[["D"]]),
              setNames(as.data.table(S_off), paste0(names(S_off), "_off")),
              setNames(as.data.table(S_on0), paste0(names(S_on0), "_onB")),
              as.data.table(all_costs(off, Wm, "_off")), as.data.table(all_costs(on0, Wm, "_onB")),
              as.data.table(setNames(as.list(unlist(lapply(names(Wm), function(m) Wm[[m]]))),
                                     unlist(lapply(names(Wm), function(m) paste0(c("win_", "win_reached_"), m))))))
  rows = list(); profs = list()
  for (mode in MATCH) {
    ## the window of this mode for a candidate Tregs ON run o; the cost is the
    ## RELATIVE damage difference of o and the Tregs OFF run over that window
    win_of = function(o) if (mode == "conv") max(tc_off, tconv(o)) else if (mode %in% names(Wm)) Wm[[mode]][["te"]] else TMAX
    rd = function(o) { te = win_of(o); Doff = dmg_to(off, te); (dmg_to(o, te) - Doff) / Doff }
    f = function(b) { o = run(p, u$Eh0, b, 1); if (is.null(o)) NA_real_ else rd(o) }
    n_eval = 0L
    ## the same level first: which way does ON have to go?
    f0 = rd(on0)
    diag = NULL
    ## 1. scan the whole range on a log grid through beta_OFF
    kk = seq(ceiling(log(BETA_LO / u$b_off) / log(GRID_STEP)), floor(log(BETA_HI / u$b_off) / log(GRID_STEP)))
    B = u$b_off * GRID_STEP^kk
    Fv = vapply(seq_along(B), function(j) if (kk[j] == 0L) f0 else f(B[j]), numeric(1))
    n_grid = length(B); n_eval = n_grid - 1L
    ## 2. probe the near-target local minima between their neighbours
    ok = which(is.finite(Fv)); n_probe = 0L
    if (length(ok) >= 3L) {
      fo = Fv[ok]; jm = which(fo[-c(1, length(fo))] <= fo[-(length(fo) - 0:1)] & fo[-c(1, length(fo))] <= fo[-(1:2)]) + 1L
      jm = jm[fo[jm] > 0 & fo[jm] < NEAR]
      for (j in jm) {
        bp = exp(seq(log(B[ok[j - 1L]]), log(B[ok[j + 1L]]), length.out = N_PROBE + 2L))[-c(1, N_PROBE + 2L)]
        bp = bp[abs(bp / B[ok[j]] - 1) > 1e-9]
        B = c(B, bp); Fv = c(Fv, vapply(bp, f, numeric(1))); n_probe = n_probe + length(bp)
      }
      n_eval = n_eval + n_probe
    }
    o_ = order(B); B = B[o_]; Fv = Fv[o_]
    prof = data.table(i = i, match = mode, beta = B, f = Fv, probe = !(B %in% (u$b_off * GRID_STEP^kk)))
    ## crossings on the scanned profile (f <= 0 counts as "not costlier")
    ok = is.finite(Fv); Bf = B[ok]; Ff = Fv[ok]; s = Ff > 0
    up_at = which(!head(s, -1) & tail(s, -1)); dn_at = which(head(s, -1) & !tail(s, -1))
    lx = function(j) exp(log(Bf[j]) + (0 - Ff[j]) / (Ff[j + 1L] - Ff[j]) * (log(Bf[j + 1L]) - log(Bf[j])))
    up_ge = up_at[Bf[up_at + 1L] > u$b_off]
    diag = data.table(n_grid = n_grid, n_probe = n_probe, n_fail_grid = sum(!ok),
                      n_cross_up = length(up_at), n_cross_down = length(dn_at),
                      b_first_up = if (length(up_ge)) lx(up_ge[1]) else NA_real_)
    ## 3. the last point that is not costlier, bisection towards the next
    k = tail(which(!s), 1)
    if (!length(k)) { status = "no_match_lo"; b_star = NA_real_ }
    else if (k == length(Ff)) { status = "no_match_hi"; b_star = NA_real_ }
    else {
      ## bisect towards the crossing itself (|f| <= TOL_BIS), never stopping at a
      ## bracket end merely inside TOL: where f is flat, every beta of a wide band
      ## is within TOL and the first one found would bias beta* downwards. The
      ## closest point found is kept and accepted if it meets TOL.
      a = log(Bf[k]); z = log(Bf[k + 1L])
      best = if (abs(Ff[k]) <= abs(Ff[k + 1L])) c(a, Ff[k]) else c(z, Ff[k + 1L])
      for (it in seq_len(MAX_BIS)) {
        m = (a + z) / 2; fm = f(exp(m)); n_eval = n_eval + 1L
        if (!is.finite(fm)) break
        if (abs(fm) < abs(best[2])) best = c(m, fm)
        if (abs(fm) <= TOL_BIS) break
        if (fm <= 0) a = m else z = m
      }
      fm = best[2]
      if (!is.finite(fm) || abs(fm) > TOL) { status = "no_match_jump"; b_star = NA_real_ }
      else {
        b_star = exp(best[1])
        status = if (abs(b_star / u$b_off - 1) < 1e-9) "matched_same" else if (b_star > u$b_off) "matched_up" else "matched_down"
        ## 4. the jump test around beta*
        fj = vapply(b_star * JUMP_X, f, numeric(1)); n_eval = n_eval + length(JUMP_X)
        diag = cbind(diag, data.table(f_m3 = fj[1], f_m1 = fj[2], f_p1 = fj[3], f_p3 = fj[4],
                                      jump_pts = 100 * (fj[4] - fj[1])))
      }
    }
    row = cbind(ref, data.table(match = mode, status = status, b_star = b_star, n_eval = n_eval))
    if (!is.null(diag)) { row = cbind(row, diag); profs[[mode]] = prof }
    if (is.finite(b_star)) {
      ons = if (status == "matched_same") on0 else run(p, u$Eh0, b_star, 1)
      if (!is.null(ons)) {
        S_on = summ(ons, t_div, u$Kg); te = win_of(ons); pm_s = pre_int(ons, t_div_M1)
        row = cbind(row, data.table(t_arr_star = arrival(ons), win = te, tconv_onS = tconv(ons),
                                    D_target = dmg_to(off, te), D_win_onS = dmg_to(ons, te), match_err = rd(ons),
                                    P_win_off = { k = which(off$time <= te); trap(off$time[k], off$PLP[k]) },
                                    P_win_onS = { k = which(ons$time <= te); trap(ons$time[k], ons$PLP[k]) },
                                    P_preM1_onS = pm_s[["P"]], D_preM1_onS = pm_s[["D"]]),
                    setNames(as.data.table(S_on), paste0(names(S_on), "_onS")),
                    as.data.table(all_costs(ons, Wm, "_onS")))
      } else row[, status := "fail_star"]
    }
    rows[[mode]] = row
  }
  out = rbindlist(rows, fill = TRUE)
  if (length(profs)) attr(out, "profile") = rbindlist(profs)
  out
}

## ---- run, flushing as it goes -------------------------------------------------------
safe = function(i) tryCatch(one(i), error = function(e) data.table(i = i, status = paste("error:", conditionMessage(e))))
IDX = if (nzchar(UNITS)) as.integer(strsplit(UNITS, ",")[[1]]) else seq_len(nrow(U))
IDX = IDX[(IDX - 1L) %% ARRAY_N == ARRAY_ID]
OUT_PROF = sub("\\.rds$", "_profile.rds", OUT)
RES = list(); PROF = list()
if (RESUME && file.exists(OUT)) {
  RES = list(readRDS(OUT)); if (file.exists(OUT_PROF)) PROF = list(readRDS(OUT_PROF))
  IDX = setdiff(IDX, RES[[1]]$i)
  cat(sprintf("[fig04 exp] RESUME: %d units in %s, %d to go\n", uniqueN(RES[[1]]$i), OUT, length(IDX)))
}
cat(sprintf("[fig04 exp] beta in [%g, %g], grid x%g, task %d of %d: %d units\n", BETA_LO, BETA_HI, GRID_STEP, ARRAY_ID, ARRAY_N, length(IDX)))
t0 = Sys.time()
for (blk in split(IDX, ceiling(seq_along(IDX) / FLUSH))) {
  r = if (N_CORES > 1L) mclapply(blk, safe, mc.cores = N_CORES) else lapply(blk, safe)
  r = r[vapply(r, is.data.frame, logical(1))]
  PROF = c(PROF, lapply(r, attr, "profile"))
  RES = c(RES, r)
  saveRDS(rbindlist(RES, fill = TRUE), OUT)
  saveRDS(rbindlist(PROF, fill = TRUE), OUT_PROF)
  cat(sprintf("[fig04 exp] %d of %d units, %.1f min\n", match(tail(blk, 1), IDX), length(IDX),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
R = rbindlist(RES, fill = TRUE)
saveRDS(R, OUT)
cat(sprintf("[fig04 exp] wrote %s: %d units\n", OUT, nrow(R)))
print(dcast(R[, .N, by = .(match, status)], status ~ match, value.var = "N", fill = 0))
