## ============================================================================
## pipeline/50_harm_rescue_compute.R -- the single computation behind
## Supplementary Fig. 3 (Treg treatments on the
## combinations Tregs harm, pathogenic_harm branch).
## ----------------------------------------------------------------------------
## WHAT IT DOES. Takes the gated harmed combinations, integrates every treatment
## arm for each of them, measures the pathogen dip on every arm, drops any
## combination the solver could not integrate, and bisects, per combination and
## treatment, the level at which the combination stops being harmed ("neutral"
## boundary) and the level at which it is helped ("help" boundary). It writes
## ONE file. Nothing downstream integrates anything except the single
## demonstration combination in the second row of Supplementary Fig. 3.
##
## WHY ONE SCRIPT. The bisection midpoints are levels the coarse ladders never
## visit (the delay ladder is 0, 10, 20, 40, 80, ... and a threshold lands at,
## e.g., 84.375), so those integrations have to happen; doing them in the same
## process as the ladder integrations means the curves, the dip statistics and
## the thresholds cannot come from different runs. The settings stamp holds
## only numerical settings and cannot detect a stale companion file, which is
## why there is only one file.
##
## WHAT THE OUTPUT HOLDS. A list of four elements:
##   sm        one row per combination per treatment arm: the endpoint states,
##             the burden difference against the Tregs OFF arm, and the dip
##             statistics for that arm
##   tr        the trajectories on the common output grid, one row per
##             combination per arm per time, keyed by cid
##   thresh    one row per combination per treatment per boundary: the bracket
##             it was found in, the bisected threshold `thr`, the burden
##             difference and outcome class measured there, and whether the
##             response was non-monotone. Two boundaries -- `neutral`, where
##             harm ends, and `help`, where the arm beats no Tregs at all
##   settings  the stamp, recorded once
##
## COST. Forty-two integrations per combination for the arms (9 delays, 12 tau_10
## multipliers, 8 tau_t_d multipliers, 12 tau_t_a multipliers, 1 Tregs OFF),
## then about seven more for every (combination, treatment, boundary) that has
## a bracket. The production run (1,773 gated combinations) took 130.5 min for
## the arms and 85.9 min for the thresholds on the cluster node
## (pipeline/cluster/harm_rescue_compute.sh).
##
## in : ./data/skeleton_pathogenic_harm.rds, ./data/tau_grid_v6_<n>.rds,
##      ./lib/skeleton_defs.R
## out: ./data/harm_rescue_data.rds
## RUN: SAMPLE=0 N_CORES=12 Rscript pipeline/50_harm_rescue_compute.R
##      SAMPLE=100 Rscript pipeline/50_harm_rescue_compute.R      # a hundred at random, to test
##      SEL_N=200 SEL_BY=plp_off_traj Rscript pipeline/50_harm_rescue_compute.R   # a stratum
##      THR_TREATS=delay SAMPLE=0 Rscript pipeline/50_harm_rescue_compute.R       # delay thresholds only
##      sbatch pipeline/cluster/harm_rescue_compute.sh                             # the whole thing on della
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages({library(data.table); library(deSolve); library(parallel)})
BRANCH = "pathogenic_harm"
## The extinction floor is not a separate overlay. It lives in
## lib/extinction_floor.R, which lib/skeleton_defs.R sources and which
## pipeline/01_sweep.R sources too, so this computation and the sweep that produced
## its labels integrate the same dynamics under the same tolerances. While the
## floor existed only here, the sweep was still classifying basins on
## trajectories that could regrow a pathogen from 1e-100.
## skeleton_defs.R supplies EXT_TOL, rootfun_ext, EVENTS_EXT, ATOL, RTOL and HMAX.
source("./lib/skeleton_defs.R")

N_CORES = as.integer(Sys.getenv("N_CORES", "4"))
SAMPLE  = as.integer(Sys.getenv("SAMPLE", "100"))
SEED    = as.integer(Sys.getenv("SEED", "1"))
TMAX    = as.numeric(Sys.getenv("TMAX", "25000"))
## ATOL, RTOL and HMAX come from lib/extinction_floor.R, sourced above, so this
## cannot drift from the sweep. They were set at 1e-14 / 1e-10 for a reason that
## has not gone away: deSolve defaults to atol = rtol = 1e-6, so any state below
## about 1e-6 is numerical noise, P_LP was returning small NEGATIVE values at its
## minimum, and the time of that minimum was then the time of the largest
## undershoot rather than a property of the system. The floor does not remove
## that requirement, it sharpens it -- the root tests P_LP against
## EXT_TOL = 1e-11, so atol must resolve a state at that level or the event fires
## on noise. rtol is 1e-8 rather than 1e-10: it binds on the large states
## (E_h ~ 100) where it costs the most, and atol is what the small ones need.
## Override both together, and only together with the sweep.
EPS     = as.numeric(Sys.getenv("EPS", "0.01"))
CTRL_FRAC = as.numeric(Sys.getenv("CTRL_FRAC", "0.01"))  # "controlled" = this share of the peak
EH0_HI  = as.numeric(Sys.getenv("EH0_HI", "99"))
## ONE THRESHOLD, TWO USES: the level at which the model calls the pathogen
## extinct and the level at which the dip statistics call it zero. It is EXT_TOL
## and it is not separately settable -- calling the pathogen "cleared" at a
## different level from the one the integration made absorbing would report a
## threshold the model never applied.
PLP_FLOOR = EXT_TOL
## EFFECT_TOL divides help from harm on the composite. It is used HERE by the
## bisection, which needs a class to bisect on, and it is stamped so the
## plotting scripts can see which value d_star was located under. They are free
## to re-derive the help/harm labels at a different tolerance; d_star cannot be
## re-derived without re-integrating, which is why it is settled here.
TOL     = as.numeric(Sys.getenv("EFFECT_TOL", "1.0"))
GAP_MIN = as.numeric(Sys.getenv("GAP_MIN", "1e-9"))
DTOL    = as.numeric(Sys.getenv("DTOL", "1"))
DELAYS  = as.numeric(strsplit(Sys.getenv("DELAYS", "0,10,20,40,80,160,320,960,2400"), ",")[[1]])
## THE THREE RATE LADDERS NOW COVER THE SAME FOLD RANGE, and that is why the
## two descending ones were extended. Supplementary Fig. 3 plots the cumulative fraction
## converted against FOLD CHANGE away from the untreated arm, so what a panel
## can show is bounded by how far its ladder goes: tau_t_d reached 128-fold and
## saturated, while tau_10 and tau_t_a stopped at 20-fold with the curves still
## climbing. Their plateaus were therefore the end of the ladder rather than a
## property of the system, and the three panels were not comparable in the
## quantity their x axis carries. Both now run to 1/133.
##
## WHAT TO EXPECT AT THE NEW END. A tau_t_a multiplier near zero induces almost
## no Treg and a tau_10 multiplier near zero leaves M1 barely deactivated, so in
## both cases the treated arm approaches the untreated OFF arm and dcomp goes to
## zero -- which is NEUTRAL, not help. The neutral curves should climb toward
## 100% and the help curves should plateau for real. If a help curve is still
## rising at 133-fold, the ladder is still the binding constraint and this note
## needs revisiting rather than the result being read off it.
SCALES  = as.numeric(strsplit(Sys.getenv("SCALES",
            "1,0.75,0.5,0.35,0.25,0.15,0.1,0.05,0.03,0.02,0.0125,0.0075"), ",")[[1]])
DEACTS  = as.numeric(strsplit(Sys.getenv("DEACTS", "1,2,4,8,16,32,64,128"), ",")[[1]])
INDUCE  = as.numeric(strsplit(Sys.getenv("INDUCE",
            "1,0.75,0.5,0.35,0.25,0.15,0.1,0.05,0.03,0.02,0.0125,0.0075"), ",")[[1]])
## SELECTING A STRATUM INSTEAD OF A RANDOM SAMPLE. SEL_N > 0 keeps the N gated
## harmed combinations with the LOWEST steady-state pathogen load UNDER TREGS OFF
## and drops the rest; SAMPLE is then ignored. These are the combinations whose
## untreated system clears the infection, so whatever burden the Tregs-ON arm
## ends with was produced by the regulatory response and not inherited from a
## pathogen the system could not control anyway.
##
## SEL_BY says which pathogen load is used:
##
##   plp_off_traj (default) P_LP at t = TMAX on the Tregs OFF arm, measured by
##                integrating it. Preferred because it is the state the figures
##                actually plot and it does not depend on the sweep's recorded
##                attractor, which is known to be unreliable for part of the
##                ON_mono class. Costs one extra integration per gated
##                combination.
##   plp_off_fp   lo_PLP_off, the sweep's low-burden fixed point under Tregs
##                OFF. Free, but inherits that unreliability.
##   plp_on_traj  P_LP at t = TMAX on the untreated Tregs ON arm.
##   plp_on_fp    lo_PLP_on, the sweep's low-burden fixed point under Tregs ON.
##
## The two Tregs OFF options select on a property of the COMBINATION and cannot
## induce a selection artefact in the treatment comparison. The two Tregs ON
## options select on the treated state, which they can.
##
## THE STRATUM IS NOT THE CLASS. Every population fraction computed downstream
## then describes the selected stratum and not the harmed class. The two
## directions are not interchangeable either: over the gated set the lowest 200
## by lo_PLP_off carry a MEDIAN HARM OF 11.5 against 6.99 for the class, while
## the lowest 200 by lo_PLP_on carry 3.45. Selecting a low Tregs OFF pathogen
## load selects the combinations Tregs damage most. The output gets its own name
## so a selected run cannot overwrite the data of a full one.
SEL_N   = as.integer(Sys.getenv("SEL_N", "0"))
SEL_BY  = Sys.getenv("SEL_BY", "plp_off_traj")
TAG     = if (SEL_N > 0) sprintf("_%s_lowest%d", SEL_BY, SEL_N) else ""
OUT_RDS = Sys.getenv("OUT_RDS", sprintf("./data/harm_rescue_data%s.rds", TAG))

## THE SETTINGS STAMP, recorded once because there is now one file to record it
## in. It says what numerical regime the data was produced under; it says
## nothing about WHICH data, which is why the three-file arrangement it came
## from could not detect a stale companion file and this one does not have to.
STAMP = function() list(ATOL = ATOL, RTOL = RTOL, TMAX = TMAX, HMAX = HMAX,
                        EPS = EPS, EH0_HI = EH0_HI, PLP_FLOOR = PLP_FLOOR,
                        EXT_TOL = EXT_TOL, CTRL_FRAC = CTRL_FRAC, EFFECT_TOL = TOL,
                        DTOL = DTOL, SAMPLE = SAMPLE, SEL_N = SEL_N, SEL_BY = SEL_BY,
                        DELAYS = DELAYS, SCALES = SCALES, DEACTS = DEACTS,
                        INDUCE = INDUCE)

S = as.data.table(readRDS("./data/skeleton_pathogenic_harm.rds"))
## THE TOP OF THE PROBED RANGE IS EXCLUDED, and this is the reason, stated rather
## than left as a bare number. The sweep probes Eh0 up to EH0_HI = 99 and bisects
## within the bracket it finds, so a threshold returned at 95 or above sits in the
## last probe interval and is bounded by the end of the line rather than located
## by it. Those combinations are not wrong, they are unresolved at this
## resolution, and treating a boundary artefact as a measured threshold would put
## them at one end of every panel. Raise PROBE_BY resolution in the sweep if they
## matter; do not simply raise this cut.
EH0_EDGE = as.numeric(Sys.getenv("EH0_EDGE", "95"))
.n0 = nrow(S); S = S[is.finite(Eh0_off) & Eh0_off < EH0_EDGE]
cat(sprintf("[harm] dropped %d of %d combinations with Eh0_off >= %g (unresolved at the top of the probed range)\n",
            .n0 - nrow(S), .n0, EH0_EDGE))
S = S[is.na(dir_on) | dir_on != "decreasing"]
S = S[rescue_mode != "bistable" | (is.finite(rescue_gap) & rescue_gap >= GAP_MIN)]
GATED = copy(S)
cat(sprintf("[harm] %s harmed combinations pass the gate\n", format(nrow(S), big.mark = ",")))

TT = c(0, exp(seq(log(0.05), log(TMAX), length.out = 300)))
K_ECAP = STRUCTURAL$K_E
comp = function(Eh, PLP) sqrt((K_ECAP - Eh)^2 + max(PLP, 0)^2)
pars = function(cc) {
  trow = tg[.(cc$batch, cc$row_ind), nomatch = NULL]
  if (!nrow(trow)) return(NULL)
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = cc$Kg; p$beta_ros = cc$beta_ros; p$r_L_P = cc$r_L_P; p$r_LP_P = cc$r_L_P
  p
}
## One arm. th = 0 is Tregs off throughout; otherwise Tregs are induced from
## t = d onward, which for d > 0 is integrated as two legs from the same state.
arm = function(p, e0, th, d = 0) {
  if (th == 0 || d <= 0) {
    o = tryCatch(ode(seed(e0), TT, rhs_app, c(p, list(theta_treg = th)),
                     method = "lsodar", maxsteps = 1e5, hmax = HMAX,
                     atol = ATOL, rtol = RTOL,
                     rootfun = rootfun_ext, events = EVENTS_EXT), error = function(e) NULL)
    if (is.null(o) || nrow(o) < length(TT)) return(NULL)
    return(as.data.table(o))
  }
  t1 = c(TT[TT < d], d)
  o1 = tryCatch(ode(seed(e0), t1, rhs_app, c(p, list(theta_treg = 0)),
                    method = "lsodar", maxsteps = 1e5, hmax = HMAX,
                     atol = ATOL, rtol = RTOL,
                     rootfun = rootfun_ext, events = EVENTS_EXT), error = function(e) NULL)
  if (is.null(o1) || nrow(o1) < length(t1)) return(NULL)
  y1 = setNames(as.numeric(o1[nrow(o1), STATE_ORDER]), STATE_ORDER)
  t2 = c(d, TT[TT > d])
  if (length(t2) < 2) return(as.data.table(o1))
  o2 = tryCatch(ode(y1, t2, rhs_app, c(p, list(theta_treg = 1)),
                    method = "lsodar", maxsteps = 1e5, hmax = HMAX,
                     atol = ATOL, rtol = RTOL,
                     rootfun = rootfun_ext, events = EVENTS_EXT), error = function(e) NULL)
  if (is.null(o2) || nrow(o2) < length(t2)) return(NULL)
  rbind(as.data.table(o1), as.data.table(o2)[-1])
}
## The trajectory is returned on the common grid so medians across combinations
## are taken at the same times; a delayed run has an extra knot at t = d, which
## is dropped rather than shifting every later point by one.
## Tr and S are kept because they are the two steps of the route the treatments
## act on: Tr is induced at rate theta_treg*tau_t_a*(M1*C_rat)*T0, S is produced
## from Tr at sigma_S*Tr, and S removes M1 at tau_10*S*M1. The delay moves when
## that route starts; the tau_10 multiplier moves only its last step, so Tr and S
## are free to respond to the M1 that the multiplier leaves in place. Inside the
## frame S is the column, which shadows the skeleton table of the same name in
## the global environment.
oncommon = function(d) d[time %in% TT][match(TT, time), .(time, M1, PLP, Eh, Tr, S)]

## THE PATHOGEN DIP, ON A FLOORED TRAJECTORY. P_LP rises to a peak and falls;
## where the arm does not eliminate the pathogen it rises again and the minimum
## between the two is the dip. Where it does eliminate, the decline continues
## through values with no meaning -- 1e-110 to 1e-120 and beyond -- and the
## unfloored minimum lands wherever the integration happens to stop, which
## reports a dip time far later than the elimination it is meant to mark. So
## everything below PLP_FLOOR is set to zero first. The minimum of the floored
## series is then attained at the FIRST time the pathogen reaches the floor,
## which is the quantity wanted, and which.min returns that first index.
##
## t_zero is that time stated explicitly, and is NA where the arm never reaches
## the floor. t_ctrl is the first time after the peak at which P_LP falls below
## CTRL_FRAC of that peak, which is defined whether or not there is a rebound.
dipstat = function(d) {
  n = nrow(d); z = d$PLP; z[z < PLP_FLOOR] = 0
  k = which.max(z)
  j = k - 1L + which.min(z[k:n])
  i = which(z[k:n] < CTRL_FRAC * z[k])
  iz = which(z[k:n] <= 0)
  list(t_PLP_peak = d$time[k], PLP_peak = z[k],
       t_dip = d$time[j], PLP_dip = z[j],
       t_ctrl = if (length(i))  d$time[k - 1L + i[1]]  else NA_real_,
       t_zero = if (length(iz)) d$time[k - 1L + iz[1]] else NA_real_)
}

one = function(i) {
  cc = S[i]; p = pars(cc)
  if (is.null(p)) return(NULL)
  e0 = cc$Eh0_off + EPS * (EH0_HI - cc$Eh0_off)
  off = arm(p, e0, 0)
  if (is.null(off)) return(NULL)
  no = nrow(off); co = comp(off$Eh[no], off$PLP[no])
  kp = which.max(off$M1)
  base = list(M1_peak_off = off$M1[kp], M1_floor_off = min(off$M1[kp:no]),
              PLP_off = max(off$PLP[no], 0), Eh_off = off$Eh[no], comp_off = co)
  run = function(on, treat, level) {
    if (is.null(on)) return(NULL)
    n = nrow(on); k = which.max(on$M1); cn = comp(on$Eh[n], on$PLP[n])
    ds = dipstat(on)
    ## cid IS THE COMBINATION, AND IT HAS TO BE IN `tr` TOO. The trajectory table
    ## used to carry only treat, level, time and the states, because every panel
    ## draws a MEDIAN over all combinations and never needed to know which row
    ## came from which. Anything that draws a SUBSET does need to know -- the
    ## extinct-only time courses could not be written
    ## without it. One integer rather than the five key columns, because `tr` has
    ## a row per combination per treatment level per time point and five extra
    ## columns there would inflate the trajectory table for no gain; `sm` carries
    ## both, so cid maps back to the key whenever that is wanted.
    list(sm = data.table(cid = i,
           batch = cc$batch, row_ind = cc$row_ind, Kg = cc$Kg,
           beta_ros = cc$beta_ros, r_L_P = cc$r_L_P, treat = treat, level = level,
           M1_peak = on$M1[k], M1_floor = min(on$M1[k:n]), t_M1_peak = on$time[k],
           PLP_end = max(on$PLP[n], 0), Eh_end = on$Eh[n], comp_on = cn,
           dcomp = cn - co, M1_peak_off = base$M1_peak_off,
           M1_floor_off = base$M1_floor_off, PLP_off = base$PLP_off,
           Eh_off = base$Eh_off,
           ## per-arm quantities the main figure and the collapse test need
           t_PLP_peak = ds$t_PLP_peak, PLP_peak = ds$PLP_peak,
           t_dip = ds$t_dip, PLP_dip = ds$PLP_dip, t_ctrl = ds$t_ctrl,
           t_zero = ds$t_zero,
           S_end = on$S[n], Tr_end = on$Tr[n],
           ## PL is recorded because P_LP = 0 is not invariant on its own: the
           ## influx term PL*gamma_B_P*(1 - Eh/K_E)*(1 - PLP/K_LP_P) has no factor
           ## of PLP, so elimination is permanent only once PL is also zero
           PL_end = on$PL[n]),
         tr = cbind(data.table(cid = i, treat = treat, level = level), oncommon(on)))
  }
  R = c(lapply(DELAYS, function(d) run(arm(p, e0, 1, d), "delay", d)),
        lapply(SCALES, function(s) { q = p; q$tau_10 = p$tau_10 * s
                                     run(arm(q, e0, 1, 0), "strength", s) }),
        lapply(DEACTS, function(r) { q = p; q$tau_t_d = p$tau_t_d * r
                                     run(arm(q, e0, 1, 0), "deact", r) }),
        lapply(INDUCE, function(a) { q = p; q$tau_t_a = p$tau_t_a * a
                                     run(arm(q, e0, 1, 0), "induce", a) }),
        list(run(off, "off", NA_real_)))
  R = R[!vapply(R, is.null, logical(1))]
  list(sm = rbindlist(lapply(R, `[[`, "sm")),
       tr = rbindlist(lapply(R, `[[`, "tr")))
}
## THE SELECTION, done here because the trajectory option needs arm() and pars().
if (SEL_N > 0 && SEL_N < nrow(S)) {
  if (SAMPLE > 0) cat("[harm] SEL_N is set, so SAMPLE is ignored\n")
  ## th = 0 is the Tregs OFF arm, th = 1 the untreated Tregs ON arm; both are
  ## integrated from the same initial condition the main loop uses
  endPLP = function(th) unlist(mclapply(seq_len(nrow(S)), function(i) {
      cc = S[i]; p = pars(cc)
      if (is.null(p)) return(NA_real_)
      a = arm(p, cc$Eh0_off + EPS * (EH0_HI - cc$Eh0_off), th, 0)
      if (is.null(a)) NA_real_ else max(a$PLP[nrow(a)], 0)
    }, mc.cores = N_CORES))
  key = switch(SEL_BY,
    plp_off_traj = endPLP(0),
    plp_on_traj  = endPLP(1),
    plp_off_fp   = S$lo_PLP_off,
    plp_on_fp    = S$lo_PLP_on,
    stop("SEL_BY must be plp_off_traj, plp_on_traj, plp_off_fp or plp_on_fp"))
  S[, sel_key := key]
  S = S[is.finite(sel_key)][order(sel_key)]
  keep = min(SEL_N, nrow(S)); cut_at = S$sel_key[keep]
  ## ties at the boundary make the last places arbitrary, so they are reported
  n_tie = sum(S$sel_key == cut_at); n_in = sum(S$sel_key[seq_len(keep)] == cut_at)
  S = S[seq_len(keep)]
  cat(sprintf("[harm] SELECTED the %d gated combinations with the lowest %s\n", nrow(S), SEL_BY))
  cat(sprintf("     selection value: median %.3g, range %.3g to %.3g\n",
              median(S$sel_key), min(S$sel_key), max(S$sel_key)))
  cat(sprintf("     %d combinations sit exactly at the cut of %.3g and %d of them are kept, so the last places are arbitrary\n",
              n_tie, cut_at, n_in))
  cat(sprintf("     harm size d_lo: median %.3g in the stratum against %.3g over the whole gated set\n",
              median(S$d_lo), median(GATED$d_lo)))
  cat("     rescue_mode in the stratum: ")
  print(table(S$rescue_mode))
  cat("[harm] every fraction below describes THIS STRATUM, not the harmed class\n")
} else if (SAMPLE > 0 && SAMPLE < nrow(S)) {
  set.seed(SEED); S = S[sample(.N, SAMPLE)]
}
cat(sprintf("[harm] %s combinations; %d delays, %d strengths, %d deactivation, %d induction, %d integrations each\n",
            format(nrow(S), big.mark = ","), length(DELAYS), length(SCALES), length(DEACTS),
            length(INDUCE), 1 + length(DELAYS) + length(SCALES) + length(DEACTS) + length(INDUCE)))

t0 = Sys.time()
RES = mclapply(seq_len(nrow(S)), function(i) tryCatch(one(i), error = function(e) NULL),
               mc.cores = N_CORES)
RES = RES[!vapply(RES, is.null, logical(1))]
if (!length(RES)) stop("every combination failed")
SM = rbindlist(lapply(RES, `[[`, "sm"), fill = TRUE)
TR = rbindlist(lapply(RES, `[[`, "tr"), fill = TRUE)
cat(sprintf("[harm] %d combinations in %.1f min\n",
            uniqueN(SM[, .(batch, row_ind, Kg, beta_ros, r_L_P)]),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
## ---- ARMS THAT ARE NOT FINITE: FOUND, NAMED AND DROPPED -------------------
## arm() rejects an integration that returns too few rows. It cannot reject one
## that returns the right number of rows carrying Inf or NaN, which is what
## lsodar does when the solution diverges and the solver still reports success.
## Such an arm reaches the tables intact, and because M1_floor_off, PLP_off,
## Eh_off and comp_off are per-combination constants repeated across all 41
## levels, ONE bad untreated arm turns every population median in the file into
## NA and the figure block then fails at quantile() with "missing values and
## NaN's not allowed if 'na.rm' is FALSE".
##
## THE COMBINATION IS DROPPED AND NAMED, NOT REPAIRED, and na.rm is not used
## anywhere to hide it. A trajectory the solver could not integrate is not a
## measurement, and a summary taken with na.rm would describe a population the
## script cannot enumerate. The keys are printed so the combination can be
## looked up in skeleton_pathogenic_harm.rds and integrated on its own.
##
## THREE COLUMNS ARE LEGITIMATELY NOT FINITE and are excluded by name rather
## than by a blanket na.rm: `level` is NA on the untreated arm, which carries no
## dose; `t_ctrl` is NA where P_LP never falls below CTRL_FRAC of its peak; and
## `t_zero` is NA where the arm never reaches the extinction floor. Everything
## else in both tables is a state or a time and has to be a real number.
##
## IT RUNS BEFORE THE BISECTION, so a combination the solver could not carry is
## never handed to REFINE below, where a non-finite composite makes helps()
## return NA and the bisection breaks out on its first midpoint, leaving a
## d_star that is the bracket end rather than a threshold.
##
## The key columns are excluded as well: they are the combination's identity,
## carried through from the grid, not anything the integration produced
CHK_SM = setdiff(names(SM)[vapply(SM, is.numeric, logical(1))],
                 c("level", "t_ctrl", "t_zero",
                   "cid", "batch", "row_ind", "Kg", "beta_ros", "r_L_P"))
CHK_TR = intersect(c("time", "M1", "PLP", "Eh", "Tr", "S"), names(TR))
nonfin = function(dt, cols)
  Reduce(`|`, lapply(dt[, cols, with = FALSE], function(x) !is.finite(x)))
SM[, .bad := nonfin(SM, CHK_SM)]
TR[, .bad := nonfin(TR, CHK_TR)]
BADC = sort(unique(c(SM[.bad == TRUE, cid], TR[.bad == TRUE, cid])))
if (length(BADC)) {
  cat(sprintf("\n[harm] %d of %d combinations carry a non-finite state and are dropped\n",
              length(BADC), uniqueN(SM$cid)))
  cat("[harm] the arms involved:\n")
  print(SM[cid %in% BADC & .bad == TRUE, .(arms = .N), by = treat][order(-arms)])
  cat("[harm] which columns, per combination and arm:\n")
  BM = melt(SM[cid %in% BADC & .bad == TRUE], id.vars = c("cid", "treat", "level"),
            measure.vars = CHK_SM, variable.name = "column", value.name = "val")
  print(BM[!is.finite(val), .(columns = paste(unique(column), collapse = " ")),
           by = .(cid, treat, level)][seq_len(min(12L, .N))])
  cat("[harm] where each trajectory leaves the reals:\n")
  print(TR[cid %in% BADC & .bad == TRUE,
           .(first_bad_time = min(time), bad_points = .N),
           by = .(cid, treat, level)][order(cid, first_bad_time)][seq_len(min(12L, .N))])
  cat("[harm] the keys, for looking them up in skeleton_pathogenic_harm.rds:\n")
  KK5 = c("batch", "row_ind", "Kg", "beta_ros", "r_L_P")
  BADKEY = unique(SM[cid %in% BADC, ..KK5])
  print(unique(SM[cid %in% BADC, .(cid, batch, row_ind, Kg, beta_ros, r_L_P)]))
  SM = SM[!cid %in% BADC]; TR = TR[!cid %in% BADC]
  if (!nrow(SM)) stop("every combination carried a non-finite state")
  ## S IS RESTRICTED TOO. The bisection indexes back into S by row, so S and the
  ## integrated set have to be the same set of combinations or refine_one() will
  ## read the parameters of a different one; and every count this script reports
  ## from S would otherwise describe a population the output does not contain.
  .nS = nrow(S); S = S[!BADKEY, on = KK5]
  cat(sprintf("[harm] %d combinations remain; the gated set drops %d to %d to match\n",
              uniqueN(SM$cid), .nS - nrow(S), nrow(S)))
}
SM[, .bad := NULL]; TR[, .bad := NULL]
## THE THREE ANCHORS MUST AGREE. delay 0, strength 1 and deactivation 1 are the
## same untreated Tregs ON arm reached by three different code paths; if they
## differ, one of the treatments is not doing what it says.
KK = c("batch", "row_ind", "Kg", "beta_ros", "r_L_P")
anch = function(tt, lv, nm) {
  z = SM[treat == tt & level == lv][, c(KK, "dcomp"), with = FALSE]
  setnames(z, "dcomp", nm)
}
aa = Reduce(function(x, y) merge(x, y, by = KK),
            list(anch("delay", 0, "d1"), anch("strength", 1, "d2"),
                 anch("deact", 1, "d3"), anch("induce", 1, "d4")))
cat(sprintf("[harm] untreated arm reached four ways: max |difference| %.3g over %d combinations\n",
            max(abs(c(aa$d1 - aa$d2, aa$d1 - aa$d3, aa$d1 - aa$d4))), nrow(aa)))

KKR = c("batch", "row_ind", "Kg", "beta_ros", "r_L_P")
CO  = SM[treat == "off", c(KKR, "comp_on"), with = FALSE]
setnames(CO, "comp_on", "comp_off_val")

## ============================================================================
## THE BISECTED THRESHOLDS: TWO PER TREATMENT, FOUR TREATMENTS
## ----------------------------------------------------------------------------
## WHAT IS LOCATED. For each combination and each treatment, the level at which
## the outcome stops being harm, and the level at which it becomes help. They
## are two different questions and both are wanted:
##
##   neutral   the first level with dcomp <= +EFFECT_TOL -- the intervention no
##             longer costs the system anything measurable
##   help      the first level with dcomp <  -EFFECT_TOL -- the intervention
##             leaves the system better off than having no Tregs at all
##
## WHY BOTH. Locating only the onset of help drops every combination an
## intervention lifts out of harm without carrying it past the help threshold.
## Those combinations then have no threshold to plot and vanish from the figure
## even though the treatment did something to them. Locating only the end of
## harm loses the distinction between the two classes, which is what the
## figure's colour is. The two coincide wherever the response steps straight
## across the dead band, which is what the delay response does; they separate
## wherever it is graded, which is what the three rate responses do.
##
## WHY ALL FOUR TREATMENTS. Without this, panel A plots a bisected delay in
## time units while B, C and D plot an index into a coarse grid of multipliers.
## The grid quantises upward -- the reported level is the first one tested past
## the threshold, not the threshold -- so the two kinds of panel are not
## comparable and neither is a measurement of the same sort. Bisecting all four
## puts every panel on the same footing.
##
## HOW A MULTIPLIER IS BISECTED. The delay ladder is arithmetic and is bisected
## arithmetically to DTOL in time units. The three multiplier ladders are
## roughly geometric (1, 0.75, 0.5 ... 0.05 and 1, 2, 4 ... 128), so they are
## bisected at the GEOMETRIC midpoint and converge on a RELATIVE tolerance,
## LTOL. An absolute tolerance would be meaningless across a ladder spanning
## two orders of magnitude: 0.01 is a quarter of the gap at the bottom of the
## tau_10 ladder and a rounding error at the top of the tau_t_d one.
##
## DIRECTION IS TAKEN FROM THE LADDER, NOT FROM THE SIGN. Intervention
## increases as the delay and the tau_t_d multiplier rise and as the tau_10 and
## tau_t_a multipliers fall, so the levels are ordered by increasing
## intervention and the bracket is [previous level, first crossing level] in
## that order. The bisection never compares magnitudes, only positions in that
## order, so it cannot be broken by a ladder that runs downward.
THR_TREATS = strsplit(Sys.getenv("THR_TREATS", "delay,strength,deact,induce"), ",")[[1]]
THR_BOUNDS = strsplit(Sys.getenv("THR_BOUNDS", "neutral,help"), ",")[[1]]
LTOL       = as.numeric(Sys.getenv("LTOL", "0.01"))   # relative, for multipliers

## the levels of each treatment ordered from the untreated arm outward, and the
## arm each level builds. tau_10 scales the S-driven deactivation of M1, tau_t_d
## the spontaneous Treg deactivation and tau_t_a the Treg induction rate.
TRT_LV = list(
  delay    = sort(DELAYS),
  strength = sort(SCALES, decreasing = TRUE),
  deact    = sort(DEACTS),
  induce   = sort(INDUCE, decreasing = TRUE))
TRT_ARM = list(
  delay    = function(p, e0, v) arm(p, e0, 1, v),
  strength = function(p, e0, v) { q = p; q$tau_10  = p$tau_10  * v; arm(q, e0, 1, 0) },
  deact    = function(p, e0, v) { q = p; q$tau_t_d = p$tau_t_d * v; arm(q, e0, 1, 0) },
  induce   = function(p, e0, v) { q = p; q$tau_t_a = p$tau_t_a * v; arm(q, e0, 1, 0) })
## geometric for a multiplier, arithmetic for a time
TRT_GEO = c(delay = FALSE, strength = TRUE, deact = TRUE, induce = TRUE)
BOUND_F = list(neutral = function(z) z <= TOL, help = function(z) z < -TOL)

bad_thr = function() data.table(thr = NA_real_, dcomp_at_thr = NA_real_,
                                class_at_thr = NA_character_)

## One combination, one treatment, one boundary. `a` is the last level on the
## untreated side of the crossing and `b` the first past it, in intervention
## order; the bisection moves whichever end the midpoint belongs to and returns
## `b`, so the value reported is always on the satisfied side of the boundary.
thr_one = function(tt, bnd, a, b, p, e0, comp_off) {
  mk = TRT_ARM[[tt]]; geo = TRT_GEO[[tt]]; pf = BOUND_F[[bnd]]
  dc = function(v) {
    x = mk(p, e0, v)
    if (is.null(x)) return(NA_real_)
    n = nrow(x); comp(x$Eh[n], x$PLP[n]) - comp_off
  }
  done = function(a, b) if (geo) abs(log(b / a)) <= log1p(LTOL) else abs(b - a) <= DTOL
  while (!done(a, b)) {
    m = if (geo) sqrt(a * b) else 0.5 * (a + b)
    z = dc(m)
    if (is.na(z)) break
    if (pf(z)) b = m else a = m
  }
  z = dc(b)
  data.table(thr = b, dcomp_at_thr = z,
             class_at_thr = if (is.na(z)) NA_character_
                            else if (z < -TOL) "Harm to help"
                            else if (z > TOL)  "Harm to harm"
                            else               "Harm to neutral")
}

S[, .row := .I]
t0r = Sys.time()
THR = list()
for (TT_ in THR_TREATS) {
  if (!TT_ %in% names(TRT_LV)) { cat(sprintf("[harm] unknown treatment %s -- skipped\n", TT_)); next }
  lv  = TRT_LV[[TT_]]
  ## the levels in intervention order, as an index, so the bracket is positional
  A   = SM[treat == TT_, c(KKR, "level", "dcomp"), with = FALSE]
  A[, ord := match(level, lv)]
  setorderv(A, c(KKR, "ord"))
  ## only combinations harmed at the untreated level; without this filter a
  ## combination already helped untreated brackets on itself and returns the
  ## untreated level, which is not a level at which anything changes
  H0 = A[ord == 1L & dcomp > TOL, ..KKR]
  A  = merge(A, H0, by = KKR)
  if (!nrow(A)) { cat(sprintf("[harm] %s: nothing harmed at the untreated level\n", TT_)); next }
  ## neutral first, because the help pass can reuse its answer -- see below
  DONE_B = list()
  for (BND in intersect(c("neutral", "help"), THR_BOUNDS)) {
    pf = BOUND_F[[BND]]
    BR = A[, { i = which(pf(dcomp))
               if (!length(i) || i[1] == 1L) .(a = NA_real_, b = NA_real_, nonmono = FALSE)
               else .(a = level[i[1] - 1L], b = level[i[1]],
                      nonmono = length(i) != .N - i[1] + 1L) }, by = KKR]
    BR = merge(BR[is.finite(b)], CO, by = KKR)
    BR = merge(BR, S[, c(KKR, ".row"), with = FALSE], by = KKR)
    if (!nrow(BR)) {
      cat(sprintf("[harm] %-8s / %-7s : no combination has a bracket\n", TT_, BND)); next }
    OUT = cbind(BR[, c(KKR, "a", "b", "nonmono"), with = FALSE],
                treat = TT_, bound = BND,
                thr = NA_real_, dcomp_at_thr = NA_real_, class_at_thr = NA_character_)
    ## WHERE THE NEUTRAL THRESHOLD IS ALREADY IN HELP, IT IS THE HELP THRESHOLD.
    ## The help predicate is strictly stronger than the neutral one, so its
    ## crossing is at or after the neutral crossing; if the outcome measured AT
    ## the neutral crossing is already help, the two are the same point and
    ## bisecting again integrates a few thousand trajectories to rediscover a
    ## number already in hand. That is the usual case for the delay treatment,
    ## whose response steps across the whole dead band in one level.
    if (BND == "help" && !is.null(DONE_B[["neutral"]])) {
      RU = DONE_B[["neutral"]][class_at_thr == "Harm to help",
                               c(KKR, "thr", "dcomp_at_thr", "class_at_thr"), with = FALSE]
      if (nrow(RU)) OUT[RU, on = KKR, `:=`(thr = i.thr, dcomp_at_thr = i.dcomp_at_thr,
                                           class_at_thr = i.class_at_thr)]
    }
    todo = which(is.na(OUT$thr))
    cat(sprintf("[harm] %-8s / %-7s : %d of %d harmed have a bracket, %d non-monotone, %d bisected, %d reused\n",
                TT_, BND, nrow(BR), uniqueN(A[, ..KKR]), sum(BR$nonmono),
                length(todo), nrow(BR) - length(todo)))
    if (length(todo)) {
      RF = rbindlist(mclapply(todo, function(k) {
             r = BR[k]; cc = S[r$.row]; p = pars(cc)
             if (is.null(p)) return(bad_thr())
             e0 = cc$Eh0_off + EPS * (EH0_HI - cc$Eh0_off)
             tryCatch(thr_one(TT_, BND, r$a, r$b, p, e0, r$comp_off_val),
                      error = function(e) bad_thr())
           }, mc.cores = N_CORES), fill = TRUE)
      OUT[todo, `:=`(thr = RF$thr, dcomp_at_thr = RF$dcomp_at_thr,
                     class_at_thr = RF$class_at_thr)]
    }
    DONE_B[[BND]] = OUT
    THR[[length(THR) + 1L]] = OUT
  }
}
THRESH = if (length(THR)) rbindlist(THR, fill = TRUE) else
  data.table(batch = integer(0), row_ind = integer(0), Kg = numeric(0),
             beta_ros = numeric(0), r_L_P = numeric(0), a = numeric(0), b = numeric(0),
             nonmono = logical(0), treat = character(0), bound = character(0),
             thr = numeric(0), dcomp_at_thr = numeric(0), class_at_thr = character(0))
cat(sprintf("\n[harm] thresholds: %d located in %.1f min\n", nrow(THRESH),
            as.numeric(difftime(Sys.time(), t0r, units = "mins"))))
if (nrow(THRESH)) {
  cat("[harm] how many per treatment and boundary, and the class each lands in:\n")
  print(dcast(THRESH[, .N, by = .(treat, bound, class_at_thr)],
              treat + bound ~ class_at_thr, value.var = "N", fill = 0L))
  cat("      'Harm to harm' means the bisection converged to a level that is still\n")
  cat("      harmful, which can only happen on a non-monotone response.\n")
  ## HOW FAR THE BISECTION MOVED THE ANSWER, as a fraction of the bracket it
  ## started in. A plain ratio thr/b cannot be read across treatments: two of
  ## the ladders descend, so there the true threshold is numerically LARGER
  ## than the grid point that first crossed it and the ratio flips below 1
  ## while meaning the same thing. The fraction is direction-free -- 0 means
  ## the grid point was already the threshold, 1 means the threshold sits at
  ## the previous level -- and it is what says whether the coarse grid was
  ## good enough.
  THRESH[, moved := {
    g = TRT_GEO[treat]
    num = fifelse(g, abs(log(thr / b)), abs(thr - b))
    den = fifelse(g, abs(log(a / b)),   abs(a - b))
    fifelse(den > 0, num / den, NA_real_) }]
  cat("[harm] how far the bisection moved the threshold inside its bracket\n")
  cat("      (0 = the coarse grid point was already right, 1 = the level before it):\n")
  print(THRESH[, .(n = .N, moved_med = round(median(moved, na.rm = TRUE), 3),
                   moved_q90 = round(quantile(moved, .9, na.rm = TRUE), 3)),
               by = .(treat, bound)][order(treat, bound)])
  ## DO THE TWO BOUNDARIES COINCIDE? They do wherever the response steps across
  ## the dead band in one level, which is the signature of a switch rather than
  ## a dose. Reported per treatment, because it is the difference between the
  ## delay result and the three rate results.
  W = dcast(THRESH[, .(thr = thr[1]), by = .(batch, row_ind, Kg, beta_ros, r_L_P, treat, bound)],
            batch + row_ind + Kg + beta_ros + r_L_P + treat ~ bound, value.var = "thr")
  if (all(c("neutral", "help") %in% names(W))) {
    W = W[is.finite(neutral) & is.finite(help)]
    if (nrow(W)) {
      cat("\n[harm] combinations with BOTH thresholds, and whether they coincide:\n")
      print(W[, .(n = .N, same = sum(abs(log(help / neutral)) < 1e-6),
                  sep_med = round(median(abs(log10(help / neutral))), 4)),
              by = treat][order(treat)])
      cat("      same = the two boundaries are the same level, so the response\n")
      cat("      crosses the dead band in one step rather than passing through it.\n")
    }
  }
}
## ---- ONE FILE -------------------------------------------------------------
## sm, tr and dstar together with the stamp that covers all three. They were
## produced in the same process from the same integrations, so no consistency
## check between them is needed or possible -- which is the point.
saveRDS(list(sm = SM, tr = TR, thresh = THRESH, settings = STAMP()), OUT_RDS)
cat(sprintf("\n[harm] wrote %s\n", OUT_RDS))
cat(sprintf("      sm    %s rows, %d combinations, %d treatment arms\n",
            format(nrow(SM), big.mark = ","), uniqueN(SM$cid), uniqueN(SM$treat)))
cat(sprintf("      tr    %s rows on %d output times\n",
            format(nrow(TR), big.mark = ","), uniqueN(TR$time)))
cat(sprintf("      thresh %d rows: %d combinations over %d treatments and %d boundaries, %d non-monotone\n",
            nrow(THRESH), uniqueN(THRESH[, .(batch, row_ind, Kg, beta_ros, r_L_P)]),
            uniqueN(THRESH$treat), uniqueN(THRESH$bound), sum(THRESH$nonmono)))
cat("      settings: ")
cat(paste(sprintf("%s=%s", c("ATOL", "RTOL", "TMAX", "HMAX", "EXT_TOL", "EFFECT_TOL", "DTOL"),
                  c(ATOL, RTOL, TMAX, HMAX, EXT_TOL, TOL, DTOL)), collapse = "  "), "\n")
cat("\n[harm] the three plotting scripts read this file and integrate nothing:\n")
cat("      figure_scripts/Supplementary_Fig_3.R   thresholds against the dip time, conversion curves, one combination\n")
