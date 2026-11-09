## ============================================================================
## pipeline/61_imperfect_sweep.R -- the imperfect-Treg sweep (Fig. 5, Extended Data Figs. 6 and 7):
## one cluster array task's slice of (unit, replicate) items, each item being
## the full phi x eta grid for one unit and one noise replicate.
## ----------------------------------------------------------------------------
## THE TWO IMPERFECTIONS (lib/ODE_SYSTEM_v6.R, make_ode_system(noise_fun)):
##   phi_instability  tau_t_d_eff = tau_t_d + phi * D / (sigma_D / tau_D)
##                    Treg reversion grows with the danger signal D relative to
##                    its largest attainable steady level sigma_D / tau_D.
##   treg_eff (eta)   C_rat_perceived = b_k * C_rat + (1 - b_k) * u_k on
##                    [k, k + 1), b_k ~ Bernoulli(eta), u_k ~ U(0, 1), all
##                    independent (make_noise_fun); eta = 1 is the perfect Treg.
## Induction reads C_rat_perceived in both dT0 and dTr. The untreated arm
## (theta_treg = 0) is integrated once per unit, by the replicate-1 item.
##
## SOLVER: lib/extinction_floor.R ode_ext(): deSolve lsodar with the extinction
## root on PL/PLP (EXT_TOL 1e-11), atol 1e-14, rtol 1e-8, maxsteps 2e5, and
## hmax = TSTEP = 1 (the noise is piecewise constant on unit intervals; the jump
## times are not passed to the solver as tcrit). Output every TSTEP to TMAX.
##
## RECORDED PER RUN: mean, variance and maximum of every state over the last
## TAIL_FRAC of the output rows (tail statistics; no endpoint matching, since
## under eta < 1 the trajectory does not settle on a fixed point).
##
## in : ./data/fig05_units_<branch>.rds, ./data/tau_grid_v6_*.rds, ./lib/*
## out: $OUTDIR/fig05_raw_t<task>.csv   one row per (unit, phi, eta, rep), plus
##      one theta_treg = 0 row (rep = 0) per unit
## RUN: on the cluster via pipeline/cluster/imperfect_sweep.sh (BRANCH and WAVE required)
##      locally, a smoke test:
##      BRANCH=sterile WAVE=1 N_UNIT=1 R_REP=1 PHI_N=3 ETA_N=3 TMAX=200 \
##        Rscript pipeline/61_imperfect_sweep.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages({library(data.table); library(deSolve)})
## BRANCH IS REQUIRED AND IS CHECKED AGAINST THE UNITS FILE. It selects the
## seeding through lib/skeleton_defs.R (PL0 = PLP0 = 0 if sterile, else
## PL0 = 500, PLP0 = 10). A default would let a sterile units file be integrated
## with an inoculum, so there is none.
if (!nzchar(Sys.getenv("BRANCH")))
  stop("BRANCH must be set explicitly (sterile | pathogenic_help | pathogenic_harm) -- ",
       "it decides PL0/PLP0 via lib/skeleton_defs.R and must match the units file")
assign("BRANCH", Sys.getenv("BRANCH"), envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R"))
suppressMessages(source("./lib/extinction_floor.R"))
suppressMessages(source("./lib/ODE_SYSTEM_v6.R"))
if (!exists("make_noise_fun")) stop("ODE_SYSTEM_v6.R did not provide make_noise_fun")

## slurm gives the task id 0-based by default; both spellings are accepted so the
## script runs identically from a laptop
ARRAY_N  = as.integer(Sys.getenv("ARRAY_N",  Sys.getenv("SLURM_ARRAY_TASK_COUNT", "1")))
ARRAY_ID = as.integer(Sys.getenv("ARRAY_ID", Sys.getenv("SLURM_ARRAY_TASK_ID",    "0")))
TMAX      = as.numeric(Sys.getenv("TMAX",  "12000"))   # the sweep horizon, as CHECK_load_balance.R uses
TSTEP     = as.numeric(Sys.getenv("TSTEP", "1"))
TAIL_FRAC = as.numeric(Sys.getenv("TAIL_FRAC", "0.5"))
R_REP     = as.integer(Sys.getenv("R_REP", "10"))
RESUME    = as.integer(Sys.getenv("RESUME", "1"))
PHI_N     = as.integer(Sys.getenv("PHI_N", "11"))
ETA_N     = as.integer(Sys.getenv("ETA_N", "11"))
N_UNIT    = as.integer(Sys.getenv("N_UNIT", "0"))      # 0 = this task's full slice
## WAVE selects one of the sampler's waves. A units file can hold more units than
## one 400-task array can finish inside the walltime, so 60_imperfect_sample.R stamps a
## `wave` column and each submission runs one of them. Required whenever the file
## carries more than one wave, so a submission cannot silently run the wrong half.
WAVE      = as.integer(Sys.getenv("WAVE", "0"))        # 0 = every wave in the file
SEED_BASE = as.integer(Sys.getenv("SEED_BASE", "20260917"))
KE        = as.numeric(Sys.getenv("K_E", "100"))
OUTDIR    = Sys.getenv("OUTDIR", "./data/fig05_raw")

PHI = seq(0, 1, length.out = PHI_N)
ETA = seq(0, 1, length.out = ETA_N)
TT  = seq(0, TMAX, by = TSTEP)

U = as.data.table(readRDS(Sys.getenv("IN_UNITS", sprintf("./data/fig05_units_%s.rds", BRANCH))))
setorder(U, unit)
## the units file names the branch it was drawn for; a mismatch would seed the
## run differently from the selection that produced it
if ("branch" %in% names(U)) {
  ub = unique(U$branch)
  if (length(ub) != 1L || ub != BRANCH)
    stop("BRANCH = ", BRANCH, " but the units file carries branch = ",
         paste(ub, collapse = ", "), " -- seeding and selection would disagree")
}
cat(sprintf("[fig05] BRANCH = %s | PL0 = %g PLP0 = %g | injury = the units file's Eh0 (e_mid)\n",
            BRANCH, PL0, PLP0))
if ("wave" %in% names(U)) {
  waves = sort(unique(U$wave))
  if (WAVE == 0L && length(waves) > 1L)
    stop("this units file holds waves ", paste(waves, collapse = ", "),
         " -- set WAVE to the one this submission is for, or the array would run ",
         "all of them in one walltime")
  if (WAVE > 0L) {
    if (!WAVE %in% waves)
      stop("WAVE = ", WAVE, " is not in this units file (has ", paste(waves, collapse = ", "), ")")
    U = U[wave == WAVE]
    cat(sprintf("[fig05] WAVE %d of {%s} | %d units\n", WAVE, paste(waves, collapse = ","), nrow(U)))
  }
}
## ---- the unit of work is (unit, replicate), not the whole unit ------------
## Slicing by replicate spreads a slow unit over R_REP tasks, so the makespan of
## the array is not set by the slowest whole unit.
##
## THE UNTREATED ARM IS STILL RUN ONCE PER UNIT, by the replicate-1 item, because
## theta_treg = 0 removes both Treg terms and neither phi nor eta can reach the
## trajectory -- running it per replicate would be R_REP identical rows.
##
## SEEDING IS PER (unit, replicate): set.seed(SEED_BASE + 1000 * unit + rep)
## at the start of each item, so an item is reproducible in isolation. Within an
## item the noise vectors are drawn in a fixed order: the untreated arm first
## (replicate 1 only), then the grid with phi in the outer and eta in the inner
## loop.
## ORDERED BY REPLICATE FIRST, AND THAT ORDER IS THE POINT. Tasks take a
## CONTIGUOUS slice, so ordering by (unit, rep) would put a unit's replicates in
## adjacent items and therefore mostly in the same task -- which is the whole-unit
## tail this slicing exists to break up. Ordering by (rep, unit) puts them
## nrow(U) items apart, so a slow unit contributes one item to each of three
## well-separated tasks and no task can collect all of it.
WORK = CJ(rep = seq_len(R_REP), ui = seq_len(nrow(U)))
setorder(WORK, rep, ui)
## GAP FILL. ITEMS=<csv with columns unit,rep> restricts the run to exactly those
## (unit, replicate) items, for recovering cluster items that never returned. Seeding is
## per (unit, replicate), so a recovered item reproduces the noise the lost one would have
## had and slots into the raw file as if it had never been missing. The slicing below still
## applies, so ARRAY_N tasks split the requested items between them.
ITEMS = Sys.getenv("ITEMS", "")
if (nzchar(ITEMS)) {
  want = fread(ITEMS)[, .(unit, rep)]
  want[, ui := match(unit, U$unit)]
  if (anyNA(want$ui)) stop("ITEMS: ", sum(is.na(want$ui)), " requested unit(s) are not in this branch's units file")
  n_before = nrow(WORK)
  WORK = WORK[want[, .(rep, ui)], on = .(rep, ui), nomatch = NULL]
  if (!nrow(WORK)) stop("ITEMS: none of the requested items are in this branch")
  cat(sprintf("[fig05] ITEMS: %d of %d items requested from %s\n", nrow(WORK), n_before, ITEMS))
}
## (cut() rejects breaks = 1, so a single-task run takes every item)
idx  = if (ARRAY_N > 1L) split(seq_len(nrow(WORK)), cut(seq_len(nrow(WORK)), breaks = ARRAY_N, labels = FALSE)) else
       list(seq_len(nrow(WORK)))
mine = idx[[min(ARRAY_ID + 1L, ARRAY_N)]]
if (N_UNIT > 0L) mine = mine[WORK$ui[mine] %in% head(unique(WORK$ui[mine]), N_UNIT)]
cat(sprintf("[fig05] task %d of %d | %d (unit,rep) items over %d units | grid %d x %d x %d reps | tmax %g\n",
            ARRAY_ID, ARRAY_N, length(mine), uniqueN(WORK$ui[mine]),
            PHI_N, ETA_N, R_REP, TMAX))

## ---- one run, and what is kept from it -------------------------------------
STATES = c(STATE_ORDER, "CumC")
tail_stats = function(o) {
  n_tail = max(2L, floor(nrow(o) * TAIL_FRAC))
  M = o[(nrow(o) - n_tail + 1L):nrow(o), STATE_ORDER, drop = FALSE]
  list(mean = colMeans(M), var = apply(M, 2, var), max = apply(M, 2, max))
}
run_one = function(p, noise_fun, e0) {
  y = c(seed(e0), CumC = 0)
  t0 = Sys.time()
  o = tryCatch(ode_ext(y, TT, make_ode_system(noise_fun), p, hmax = TSTEP),
               error = function(e) NULL)
  secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok = !is.null(o) && nrow(o) == length(TT) &&
       all(is.finite(o[, STATE_ORDER])) && o[nrow(o), 1] >= TMAX
  if (!ok) return(list(ok = FALSE, secs = secs, s = NULL))
  list(ok = TRUE, secs = secs, s = tail_stats(as.matrix(o)))
}
## the columns one run contributes, flattened
flat = function(s) {
  if (is.null(s)) return(as.list(setNames(rep(NA_real_, 2 * length(STATE_ORDER) + 4),
                                          c(paste0(STATE_ORDER, "_mean"),
                                            paste0(STATE_ORDER, "_var"),
                                            "Ei_mean", "Ei_var", "Ei_max", "PLP_max"))))
  c(setNames(as.list(s$mean), paste0(STATE_ORDER, "_mean")),
    setNames(as.list(s$var),  paste0(STATE_ORDER, "_var")),
    list(Ei_mean = KE - s$mean[["Eh"]], Ei_var = s$var[["Eh"]],
         Ei_max  = KE - min(s$mean[["Eh"]], s$max[["Eh"]]), PLP_max = s$max[["PLP"]]))
}

## ---- the loop --------------------------------------------------------------
## EACH UNIT IS WRITTEN AS SOON AS IT FINISHES, appended to this task's file.
## A unit on the full grid takes over an hour, so a task holding its results
## until the end loses every finished unit when it hits the walltime -- and a
## walltime kill is the failure this array is most likely to have. Appending
## costs one file handle and makes a timed-out task worth reading. Any earlier
## file for this task id is removed first, so a resubmitted task does not append
## to its own previous attempt.
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
OUT = sprintf("%s/fig05_raw_t%04d.csv", OUTDIR, ARRAY_ID)
## RESUME RATHER THAN RESTART. At R_REP = 10 a task carries ~8 items of ~2.7
## core-hours each, so the slowest tasks approach the walltime and a kill would
## otherwise throw away everything the task had finished. Items are written as
## they complete, so on restart the ones already in this task's file are skipped
## and only the remainder runs. RESUME=0 forces a clean rebuild of the file.
## The check is on the (unit, rep) pair actually present, not on a count, so it
## is correct whatever slice this task holds.
done_items = data.table(unit = integer(), rep = integer())
if (RESUME == 1L && file.exists(OUT)) {
  prev = tryCatch(fread(OUT, select = c("unit", "rep")), error = function(e) NULL)
  if (!is.null(prev) && nrow(prev)) {
    done_items = unique(prev[rep > 0L, .(unit, rep)])
    cat(sprintf("[fig05] resuming: %d (unit, rep) items already in %s\n", nrow(done_items), OUT))
  }
} else if (file.exists(OUT)) invisible(file.remove(OUT))
k = if (nrow(done_items)) 1L else 0L    # so the first new write appends, not truncates
t_start = Sys.time()

for (w in mine) {
  ii = WORK$ui[w]; rr = WORK$rep[w]
  cc = U[ii]
  if (nrow(done_items) && nrow(done_items[unit == cc$unit & rep == rr])) next
  trow = tg[.(cc$batch, cc$row_ind), nomatch = NULL]
  if (!nrow(trow)) { cat(sprintf("[fig05] unit %d: no tau row, skipped\n", cc$unit)); next }
  p0 = c(STRUCTURAL, EFFERO)
  for (nm in KIN) p0[[nm]] = trow[[nm]][1]
  p0$Kg = cc$Kg; p0$beta_ros = cc$beta_ros
  p0$r_L_P = cc$r_L_P; p0$r_LP_P = cc$r_L_P
  set.seed(SEED_BASE + 1000L * cc$unit + rr)

  rows = vector("list", 1L + PHI_N * ETA_N); j = 0L
  ## the untreated arm, once per unit, carried by the first replicate
  if (rr == 1L) {
    a = run_one(c(p0, list(theta_treg = 0, phi_instability = 0, treg_eff = 1)),
                make_noise_fun(TMAX, 1), cc$Eh0)
    j = j + 1L
    rows[[j]] = c(list(unit = cc$unit, theta_treg = 0L, phi_instability = NA_real_,
                       treg_eff = NA_real_, rep = 0L, solver_ok = a$ok, solve_secs = a$secs),
                  flat(a$s))
  }
  ## the treated arm over the grid, this replicate only
  for (ph in PHI) for (et in ETA) {
    b = run_one(c(p0, list(theta_treg = 1, phi_instability = ph, treg_eff = et)),
                make_noise_fun(TMAX, et), cc$Eh0)
    j = j + 1L
    rows[[j]] = c(list(unit = cc$unit, theta_treg = 1L, phi_instability = ph,
                       treg_eff = et, rep = rr, solver_ok = b$ok, solve_secs = b$secs),
                  flat(b$s))
  }
  D = rbindlist(rows[seq_len(j)], fill = TRUE)
  ## the unit's identity travels with every row, so a raw file is self-describing
  for (nm in c("batch", "row_ind", "tau_id", "Kg", "r_L_P", "beta_ros", "Eh0",
               "rescue_mode", "margin", "cmp_off", "cmp_on", "help"))
    D[, (nm) := cc[[nm]]]
  k = k + 1L
  fwrite(D, OUT, append = (k > 1L))
  cat(sprintf("[fig05] unit %d rep %d (%d/%d) | %d runs | %d failed | %.1f min elapsed\n",
              cc$unit, rr, k, length(mine), nrow(D), D[solver_ok == FALSE, .N],
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

if (!k) stop("task produced no units")
cat(sprintf("[fig05] wrote %s | %d units | %.1f min\n", OUT, k,
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
