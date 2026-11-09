## ============================================================================
## pipeline/21_corridor_tzero.R -- the basin boundaries at t = 0 (Supplementary Fig. 1B, C).
##
## Measures, for one combination and one coordinate X, how far you must displace
## the state from a base point before the outcome flips -- once with Tregs OFF and
## once with Tregs ON. The two distances are the two separatrices' positions along
## X, and their difference is the corridor width:
##
##      walking OUT along a damaging axis, the label pair passes three regimes
##        (low, low)    both reach the low basin      <- the base state
##        (high, low)   only ON reaches it            <- THE CORRIDOR
##        (high, high)  both reach the high basin
##      first transition = d_off,X    second = d_on,X    w_X = d_on,X - d_off,X
##
## in : ./data/rescue_<BRANCH>.rds, ./data/skeleton_<BRANCH>.rds, ./data/tau_grid_v6_<n>.rds
## out: ./data/corridor_tzero_<BRANCH><OUT_TAG>.rds   (long: one row per combination x coordinate)
## RUN: pilot     BRANCH=sterile N_SAMP_MODE=5 N_CORES=6 OUT_TAG=_pilot Rscript pipeline/21_corridor_tzero.R
##      cluster   BRANCH=sterile Rscript pipeline/21_corridor_tzero.R <n_chunks> <this_chunk>
##                (driven by pipeline/cluster/corridor_tzero_array.sh; ~15 core-min per combination)
##
## TWO DEFECTS IN THE v5 STAGE ARE FIXED HERE, both silent:
##
## 1. THE FINDER WAS CALLED WITH NO SEEDS. v5's fpset() called
##    find_fixed_points(par, ...) with the grid alone. That is the same call that
##    made 14_eh0_star_fine.R report monostable on combinations whose Eh0 line
##    plainly reached two attractors -- the grid missed a reached attractor at
##    (M1 9.7, PLP 13.5, L 5) on every one of six tested. Here the finder is seeded
##    with the final states of the Eh0 probe, exactly as pipeline/01_sweep.R does.
##
## 2. THE LABEL WAS A COMPOSITE MIDPOINT, NOT AN ATTRACTOR. v5's comment said
##    "identity of the fixed point the trajectory actually reaches" and the code
##    tested cmp <= (min + max)/2. Those agree when exactly two attractors are
##    reached and disagree otherwise. v6 labels by membership of the lowest-burden
##    CLUSTER, the same rule the sweep uses, so a combination cannot be labelled
##    one way here and the other way there.
## ============================================================================
suppressMessages({library(data.table); library(deSolve); library(parallel)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
BRANCH = Sys.getenv("BRANCH", "sterile")
source("./lib/skeleton_defs.R")

CK        = c("tau_id","Kg","beta_ros","r_L_P")
NSCAN     = as.integer(Sys.getenv("NSCAN", "80"))     # points per outward scan
BIS_ITERS = as.integer(Sys.getenv("BIS_ITERS", "30")) # halvings per wall
EXPAND_MAX= as.integer(Sys.getenv("EXPAND_MAX", "10"))# range doublings before giving up
PAD       = 2.0
TMAX      = as.numeric(Sys.getenv("TMAX", "12000"))
HMAX      = as.numeric(Sys.getenv("HMAX", "5"))
## RAY_CLUST_TOL IS 1, NOT 2, AND IT NOW MEANS WHAT THE SWEEP MEANS BY IT.
## This script used to break clusters on a gap in the sorted composite at a
## tolerance of 2, which is the rule pipeline/01_sweep.R abandoned: the composite is
## one number standing for two coordinates, so a cleared attractor and one
## carrying 27 units of pathogen can differ by 1.1 in it and be merged. The
## definition is in lib/attractor_cluster.R, sourced by both, and the
## tolerance matches the sweep's export.
RAY_CLUST_TOL = as.numeric(Sys.getenv("RAY_CLUST_TOL", "1.0"))
source("./lib/attractor_cluster.R")
PROBE_BY  = as.numeric(Sys.getenv("PROBE_BY", "10"))
EH0_LO    = as.numeric(Sys.getenv("EH0_LO", "0")); EH0_HI = as.numeric(Sys.getenv("EH0_HI", "99"))
PROBE     = unique(c(seq(EH0_LO, EH0_HI, PROBE_BY), EH0_HI))
N_CORES   = as.integer(Sys.getenv("N_CORES", "1"))
N_SAMP    = as.integer(Sys.getenv("N_SAMP", "0"))     # 0 = every gated combination
SEED      = as.integer(Sys.getenv("SEED", "1"))
## WHERE THE BASE STATE SITS. Rungs along the resolving side of the OFF
## threshold, as fractions of the way from it to the far end of the Eh0 range.
## The FIRST rung the OFF system actually resolves at is taken, so the default
## places base 2% past the threshold -- and measured on the sterile run, all
## 33,021 bistable combinations used that first rung, making it a uniform
## convention rather than a per-combination choice.
##
## It is still a convention, and the tolerance ratio d_on/d_off depends on it:
## base near the wall makes d_off small and the ratio large; base deep in the
## basin sends the ratio toward 1. Vary this to measure that dependence, e.g.
##   BASE_FRAC=0.20,0.35,0.50,0.75,1.00 OUT_TAG=_base20 N_SAMP=300 ...
BASE_FRAC = as.numeric(strsplit(Sys.getenv("BASE_FRAC",
              "0.02,0.05,0.10,0.20,0.35,0.50,0.75,1.00"), ",")[[1]])
if (!length(BASE_FRAC) || any(!is.finite(BASE_FRAC)) ||
    any(BASE_FRAC <= 0) || any(BASE_FRAC > 1) || is.unsorted(BASE_FRAC))
  stop("BASE_FRAC must be increasing values in (0, 1]")

## The coordinates swept. PLP only where a pathogen exists; in sterile it is
## identically zero and its basin span is degenerate.
SWEEP = c("CLP","Eh","Ea","En","M1","R","D")
if (BRANCH != "sterile") SWEEP = c(SWEEP, "PLP")
## SWEEP_ONLY restricts the coordinates measured. A change confined to one
## coordinate's scan does not justify re-running the others, whose results would
## come back identical -- run the affected coordinate and join it onto the
## existing file on (tau_id, Kg, beta_ros, r_L_P, coord).
SWEEP_ONLY = Sys.getenv("SWEEP_ONLY", "")
if (nzchar(SWEEP_ONLY)) {
  .want = trimws(strsplit(SWEEP_ONLY, ",")[[1]])
  .bad  = setdiff(.want, SWEEP)
  if (length(.bad)) stop("SWEEP_ONLY names coordinates not swept in BRANCH=", BRANCH,
                         ": ", paste(.bad, collapse = ", "))
  SWEEP = .want
}
## Physical ceilings, so an outward scan cannot run away. R and D are set per
## combination from the parameters that bound them.
## M1 is capped at M_TOT, not by a chosen number: setx() converts resting to
## activated, so x > M_TOT would require negative M0.
HARD0 = list(CLP=300, Eh=100, Ea=600, En=800, M1=M_TOT, PLP=300)
FP_MATCH_VARS = c("Eh","Ea","En","PLP","CLP","M1","R","D","Tr","S","PL")

snap_fp = function(endpt, fps) {
  if (is.null(fps) || !nrow(fps) || is.na(endpt["Eh"])) return(NA_integer_)
  vars = intersect(FP_MATCH_VARS, names(fps)); ev = endpt[vars]
  dd = vapply(seq_len(nrow(fps)), function(i) {
    fv = unlist(fps[i, vars]); sum(((ev - fv)/(abs(fv) + 1))^2, na.rm = TRUE) }, numeric(1))
  if (!length(dd) || all(!is.finite(dd))) return(NA_integer_)
  which.min(dd) }

tail_state = function(y, th, par) {
  ## ode_ext: lsodar + the extinction root + the sweep's atol/rtol
  ## (lib/extinction_floor.R). The corridor decides which side of a boundary a
  ## displaced state falls on, so it must run the same dynamics as the sweep that
  ## located the boundary.
  o = tryCatch(ode_ext(y, c(0, TMAX), rhs_app, c(par, list(theta_treg = th)),
                       maxsteps = 1e5), error = function(e) NULL)
  if (is.null(o) || nrow(o) < 2) return(NULL)
  if (o[nrow(o), 1] < TMAX * (1 - 1e-9)) return(NULL)
  setNames(as.numeric(o[nrow(o), STATE_ORDER]), STATE_ORDER) }

## PROGRESS. A run of this length must not go silent. mclapply FORKS, so several
## children writing to stdout interleave mid-line and the log is unreadable; each
## child appends one short atomic line to a shared file instead, and the caller
## tails it.
PROG = sprintf("./data/corridor_progress_%s%s.log", BRANCH, Sys.getenv("OUT_TAG", ""))
progress = function(i, cc, secs, note) {
  cat(sprintf("[%s] %5d/%-5d %-8s tau_id %-7s Kg %-3g beta %-7g rlp %-5g %7.1fs  %s\n",
              format(Sys.time(), "%H:%M:%S"), i, nrow(U), cc$rescue_mode, cc$tau_id,
              cc$Kg, cc$beta_ros, cc$r_L_P, secs, note),
      file = PROG, append = TRUE) }

one = function(i) {
  .t0 = Sys.time()
  cc = U[i]
  trow = tg[.(cc$batch, cc$row_ind), nomatch = NULL]
  if (!nrow(trow)) return(NULL)
  trow = trow[1]
  par = c(STRUCTURAL, EFFERO); for (nm in KIN) par[[nm]] = trow[[nm]]
  par$Kg = cc$Kg; par$beta_ros = cc$beta_ros; par$r_L_P = cc$r_L_P; par$r_LP_P = cc$r_L_P
  el  = function() as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  bad = function(nt) { progress(i, cc, el(), nt)
                       data.table(cc[, ..CK], coord = NA_character_, note = nt) }

  ## ---- fixed points, SEEDED FROM THE Eh0 PROBE (defect 1) -------------------
  Z  = lapply(PROBE, function(e0) tail_state(seed(e0), 0L, par))
  Z1 = lapply(PROBE, function(e0) tail_state(seed(e0), 1L, par))
  mkseeds = function(L) { L = L[!vapply(L, is.null, logical(1))]
    if (!length(L)) return(NULL)
    m = do.call(rbind, lapply(L, function(z) c(M1 = unname(z["M1"]), PLP = unname(z["PLP"]),
                                               L = unname(z["Ea"] + z["En"]))))
    m = m[is.finite(rowSums(m)), , drop = FALSE]; if (nrow(m)) unique(m) else NULL }
  fpset = function(th, sd) {
    z = tryCatch(as.data.frame(finder(c(par, list(theta_treg = th)), cc$Kg, cc$beta_ros,
                                      cc$r_L_P, th, extra_seeds = sd)), error = function(e) NULL)
    if (is.null(z) || !nrow(z)) return(NULL)
    z$cmp = cmpf(z$Eh, z$PLP); z = z[order(z$cmp), , drop = FALSE]; rownames(z) = NULL
    st = which(z$type == "stable"); if (!length(st)) return(NULL)
    ## the LOWEST-BURDEN CLUSTER, not the single lowest point: two stable points
    ## fp_same() calls one attractor are one attractor, and labelling by bare
    ## index would let a pair inside one cluster manufacture a wall.
    ##
    ## CLUSTERED IN (Eh, P_LP), NOT ALONG THE COMPOSITE. fp_cluster() returns
    ## ids ordered by increasing composite, so cluster 1 is the lowest-burden
    ## attractor and lo_ids is every stable point in it. The composite still
    ## orders the clusters; it no longer decides which points are in one.
    cs = z$cmp[st]
    .cl = fp_cluster(z$Eh[st], z$PLP[st], cs)
    list(fp = z, lo_ids = st[.cl == 1L], n_stable = length(st),
         n_attractor = max(.cl),
         gap = cs[length(cs)] - cs[1], lo = z[st[1], ], hi = z[st[length(st)], ]) }
  F0 = fpset(0L, mkseeds(Z)); F1 = fpset(1L, mkseeds(Z1))
  if (is.null(F0) || is.null(F1)) return(bad("no_fixed_points"))
  FS = list(`0` = F0, `1` = F1)

  ## ---- the label: attractor identity, never a cutoff (defect 2) -------------
  ## 1 = low basin, 0 = some other stable attractor, NA = the final state matched
  ## nothing stable (or the integration failed). NA is its own answer: a state
  ## resting on a saddle is in no basin and must not be read as the high one.
  lab = function(y, th) {
    z = tail_state(y, th, par); if (is.null(z)) return(NA_integer_)
    f = FS[[as.character(th)]]; j = snap_fp(z, f$fp)
    if (is.na(j) || !identical(f$fp$type[j], "stable")) return(NA_integer_)
    as.integer(j %in% f$lo_ids) }

  ## per-coordinate basin span: the OFF system's own high-minus-low separation in
  ## that coordinate. Every distance below is divided by it, so coordinates whose
  ## raw units differ by orders of magnitude become comparable.
  span = vapply(SWEEP, function(v) { d = abs(F0$hi[[v]] - F0$lo[[v]])
    if (length(d) != 1 || !is.finite(d) || d <= 0) NA_real_ else d }, numeric(1))
  HARD = HARD0
  HARD$R = 3 * (cc$beta_ros / trow$mu_R) * M_TOT
  HARD$D = 3 * (trow$sigma_D / trow$tau_D)

  ## ---- BASE STATE ----------------------------------------------------------
  ## A point on the Eh0 line that the OFF system takes to the LOW basin. Found by
  ## stepping from the OFF threshold toward the low basin -- which way that is
  ## depends on dir_off, so it is read rather than assumed. Verified by
  ## integration; no epsilon, no "barely reached" rule.
  up = !identical(cc$dir_off, "decreasing")
  e_lo = if (up) cc$Eh0_off else EH0_LO
  e_hi = if (up) EH0_HI     else cc$Eh0_off
  frac = BASE_FRAC
  e_try = if (up) e_lo + (e_hi - e_lo) * frac else e_hi - (e_hi - e_lo) * frac
  e0 = NA_real_; e_frac = NA_real_
  for (.k in seq_along(e_try)) if (identical(lab(seed(e_try[.k]), 0L), 1L)) {
    e0 = e_try[.k]; e_frac = BASE_FRAC[.k]; break }
  if (!is.finite(e0)) return(bad("no_base_in_low_basin_under_OFF"))
  ybase = seed(e0)

  ## ---- DISPLACING ONE COORDINATE ------------------------------------------
  ## For an unconstrained state, set it and leave the rest alone. M1 is NOT
  ## unconstrained: dM0 = -dM1, so M0 + M1 is conserved at whatever the initial
  ## condition sets it to. Writing y["M1"] = x while M0 stays at M_TOT creates a
  ## pool of M_TOT + x macrophages, up to double at the ceiling, and that
  ## inflated total then persists for the whole trajectory. The scan would be
  ## varying activation state and pool size together, which is not a
  ## displacement along the M1 axis. Converting resting to activated at fixed
  ## pool size is.
  setx = function(v, x) {
    y = ybase
    if (v == "M1") { y["M1"] = x; y["M0"] = max(M_TOT - x, 0) } else y[v] = x
    y }

  ## ---- ONE OUTWARD SCAN PER COORDINATE, both conditions ---------------------
  walls = function(vv) {
    hard = if (is.null(HARD[[vv]])) Inf else HARD[[vv]]
    hi = min(max(PAD * max(c(F0$hi[[vv]], F0$lo[[vv]], ybase[[vv]]), na.rm = TRUE), 1), hard)
    for (k in seq_len(EXPAND_MAX)) {
      g  = seq(as.numeric(ybase[[vv]]), hi, length.out = NSCAN)
      a0 = vapply(g, function(x) lab(setx(vv, x), 0L), integer(1))
      a1 = vapply(g, function(x) lab(setx(vv, x), 1L), integer(1))
      ok = !is.na(a0) & !is.na(a1); g2 = g[ok]; b0 = a0[ok]; b1 = a1[ok]
      if (length(g2) < 3) return(NULL)
      t0 = which(diff(b0) != 0); t1 = which(diff(b1) != 0)
      ## The two walls are refined INDEPENDENTLY. A clean single transition under
      ## one condition is a usable measurement even if the other condition's scan
      ## found none -- discarding both would throw away good ON walls because the
      ## OFF scan ran out of range.
      ref = function(j, th, aa) { lo = g2[j]; hi2 = g2[j+1]; bl = aa[j]
        for (q in seq_len(BIS_ITERS)) { m = (lo+hi2)/2; b = lab(setx(vv, m), th)
          if (is.na(b)) return(NA_real_); if (identical(b, bl)) lo = m else hi2 = m }
        (lo+hi2)/2 }
      got0 = length(t0) == 1; got1 = length(t1) == 1
      if ((got0 && got1) || hi >= hard)
        return(list(x_off = if (got0) ref(t0, 0L, b0) else NA_real_,
                    x_on  = if (got1) ref(t1, 1L, b1) else NA_real_,
                    n_both_low = sum(b0 == 1 & b1 == 1), n_corridor = sum(b0 == 0 & b1 == 1),
                    n_both_high = sum(b0 == 0 & b1 == 0),
                    n_t0 = length(t0), n_t1 = length(t1), hi_used = hi, hard = hard))
      hi = min(hi * 2, hard)
    }
    NULL }

  out = rbindlist(lapply(SWEEP, function(vv) {
    w = walls(vv); if (is.null(w)) return(NULL)
    data.table(cc[, ..CK], e0_base = e0, base_frac = e_frac, coord = vv,
               span = span[[vv]], base = as.numeric(ybase[[vv]]),
               lo_off = as.numeric(F0$lo[[vv]]), hi_off = as.numeric(F0$hi[[vv]]),
               x_off = w$x_off, x_on = w$x_on, n_t0 = w$n_t0, n_t1 = w$n_t1,
               hi_used = w$hi_used, hard = w$hard, n_both_low = w$n_both_low,
               n_corridor = w$n_corridor, n_both_high = w$n_both_high) }), fill = TRUE)
  if (!nrow(out)) return(bad("no_walls_on_any_coordinate"))
  out[, `:=`(d_off = abs(x_off - base)/span, d_on = abs(x_on - base)/span,
             d_ceil = abs(hi_used - base)/span)]
  out[, w := d_on - d_off]
  out[, `:=`(rescue_gap = cc$rescue_gap, rescue_gap_kind = cc$rescue_gap_kind,
             rescue_mode = cc$rescue_mode, gap_censored = cc$gap_censored,
             n_stable_off = F0$n_stable, n_stable_on = F1$n_stable, note = "ok")]
  ## Report WHAT was found, not just that it finished: a log full of "on-walls 0/7"
  ## is the mono_B signature and simultaneously explains the wall clock, because a
  ## combination with no ON wall never satisfies the early return and runs the
  ## range doubling to exhaustion on every coordinate.
  progress(i, cc, el(), sprintf("ok  on-walls %d/%d", sum(is.finite(out$d_on)), nrow(out)))
  out }

## ============================================================================
## THE POOL
## ============================================================================
R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", BRANCH)))
S = as.data.table(readRDS(sprintf("./data/skeleton_%s.rds", BRANCH)))
R = merge(R, unique(S[, .(tau_id, Kg, beta_ros, r_L_P, batch, row_ind)]), by = CK, all.x = TRUE)
## bi-bi AND mono_B, deliberately. They ask different questions of the same
## measurement: for bi-bi, w_X is the corridor width and is correlated against the
## rescue gap; for mono_B there is no ON threshold on the Eh0 line at all, so the
## question is whether an ON wall exists in ANY swept coordinate -- a missing d_on
## there is a finding, not a failure.
U = R[rescue_mode %in% c("bistable", "mono_B")]
## EXCLUDE names an .rds of rows already measured -- the merged round-1 output of
## an array that timed out (cluster/corridor_merge_partials.sh). Its combinations are
## dropped here, BEFORE the permutation, so every task of the new array sees the
## same remaining pool and takes a disjoint slice of it. Before sampling too, so
## N_SAMP counts what is left.
EXCLUDE = Sys.getenv("EXCLUDE", "")
if (nzchar(EXCLUDE)) {
  if (!file.exists(EXCLUDE)) stop("EXCLUDE file not found: ", EXCLUDE)
  XCL = unique(as.data.table(readRDS(EXCLUDE))[, ..CK])
  n0 = nrow(U); U = U[!XCL, on = CK]
  cat(sprintf("[corridor_tzero] EXCLUDE %s: %d combinations already measured, %d of %d remain\n",
              EXCLUDE, n0 - nrow(U), nrow(U), n0))
}
## SAMPLING. N_SAMP draws that many at random from the whole pool, which is
## proportional and therefore right for a production subsample but wrong for a
## pilot: the pool is ~87% bistable, so N_SAMP=30 yields ~4 mono_B and the mono_B
## question cannot be answered from it. N_SAMP_MODE draws that many from EACH
## rescue_mode instead. Use N_SAMP_MODE for pilots, N_SAMP or nothing for real runs.
N_SAMP_MODE = as.integer(Sys.getenv("N_SAMP_MODE", "0"))
if (N_SAMP_MODE > 0) {
  set.seed(SEED)
  U = U[, .SD[sample(.N, min(N_SAMP_MODE, .N))], by = rescue_mode]
  cat(sprintf("[corridor_tzero] SAMPLED %d per rescue_mode -> %d of %d gated combinations (seed %d)\n",
              N_SAMP_MODE, nrow(U), nrow(R), SEED))
} else if (N_SAMP > 0 && N_SAMP < nrow(U)) {
  set.seed(SEED); U = U[sample(.N, N_SAMP)]
  cat(sprintf("[corridor_tzero] SAMPLED %d of %d gated combinations at random (seed %d) -- proportional,\n     so rare modes are thin; use N_SAMP_MODE for a pilot\n",
              nrow(U), nrow(R), SEED))
}
## PERMUTE BEFORE CHUNKING. The pool arrives in sweep order -- (batch, row_ind,
## Kg, beta_ros, r_L_P) -- and rescue_mode tracks the parameters, so mono_B is
## CLUSTERED in that order rather than spread through it. A contiguous slice then
## inherits the clustering: measured on the first real array, chunk 1 drew 82
## bistable and 0 mono_B against a 12.4% base rate, which is a 2-in-100,000 event
## by chance.
##
## That is not merely untidy. mono_B is the EXPENSIVE class -- it has no ON wall,
## so walls() never satisfies its early return and runs the range doubling to
## exhaustion on every coordinate, about 27 core-min against 5 for a bistable
## combination. A task that lands on the mono_B block is ~37 h against a 24 h
## wall, and dies with nothing to show.
##
## The permutation is seeded and applied identically in every task, so each still
## takes a disjoint slice and the union is one serial run. CHUNK_SEED must be the
## same across the whole array; it is not a tuning knob.
CHUNK_SEED = as.integer(Sys.getenv("CHUNK_SEED", "20260820"))
set.seed(CHUNK_SEED); U = U[sample(.N)]

.args = commandArgs(trailingOnly = TRUE)
if (length(.args) >= 2) {
  n1 = as.integer(.args[1]); n2 = as.integer(.args[2])
  U = U[split(seq_len(nrow(U)), cut(seq_len(nrow(U)), n1, labels = FALSE))[[n2]]]
  cat(sprintf("[corridor_tzero] chunk %d/%d -> %d combinations (pool permuted, seed %d)\n",
              n2, n1, nrow(U), CHUNK_SEED))
  cat(sprintf("[corridor_tzero] this chunk: %d bistable, %d mono_B -- expect roughly the pool's mix\n",
              sum(U$rescue_mode == "bistable"), sum(U$rescue_mode == "mono_B")))
}
if (!nrow(U)) stop("empty pool for BRANCH=", BRANCH)
cat(sprintf("[corridor_tzero] base state rungs: %s\n", paste(BASE_FRAC, collapse = ", ")))
cat(sprintf("[corridor_tzero] %s | %d combinations (%s bistable, %s mono_B) | %d coordinates | %d cores\n",
            BRANCH, nrow(U), sum(U$rescue_mode == "bistable"), sum(U$rescue_mode == "mono_B"),
            length(SWEEP), N_CORES))

## PARTIAL RESULTS. A task killed at the wall used to lose everything, because
## the only write happened after the last combination. The serial path (which is
## what the cluster runs, one core per task) now writes what it has every
## FLUSH_EVERY combinations -- to a .partial file, NEVER to the final name.
## That distinction is load-bearing: cluster/corridor_merge_chunks.sh checks the CHUNK COUNT,
## so a partial file under the final name would pass that check and silently
## under-cover the pool. A .partial left behind is a task that did not finish.
## Default 1: at ~15 min per combination on della the write is free, and a task
## killed at the wall then loses at most the combination in flight.
FLUSH_EVERY = as.integer(Sys.getenv("FLUSH_EVERY", "1"))
if (file.exists(PROG)) file.remove(PROG)
cat(sprintf("[corridor_tzero] progress -> %s\n     watch with:  tail -f %s\n", PROG, PROG))
t0 = Sys.time()
## Braces are load-bearing: at top level a bare `else` on a new line is a parse
## error, because the `if` is already a complete expression by then.
OUT  = sprintf("./data/corridor_tzero_%s%s.rds", BRANCH, Sys.getenv("OUT_TAG", ""))
PART = paste0(OUT, ".partial")
L = if (N_CORES > 1) {
  rbindlist(mclapply(seq_len(nrow(U)), one, mc.cores = N_CORES), fill = TRUE)
} else {
  acc = vector("list", nrow(U))
  for (i in seq_len(nrow(U))) {
    acc[[i]] = one(i)
    if (i %% FLUSH_EVERY == 0)
      saveRDS(rbindlist(acc[seq_len(i)], fill = TRUE), PART)
  }
  rbindlist(acc, fill = TRUE)
}
el = as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(L, OUT)
if (file.exists(PART)) file.remove(PART)   # completed: the partial is now noise
## CORE-time, not wall/n. With n combinations on n cores every core does one and
## wall/n reports a tenth of the truth; core-minutes per combination is the number
## that scales to a different core count or to the cluster.
cat(sprintf("\n[corridor_tzero] %.1f min wall on %d cores for %d combinations\n", el, N_CORES, nrow(U)))
cat(sprintf("[corridor_tzero] %.1f core-min per combination -- multiply by the pool size, then\n     divide by the cores you will actually have\n",
            el * N_CORES / max(nrow(U), 1)))
cat(sprintf("[corridor_tzero] wrote %s (%s rows)\n", OUT, format(nrow(L), big.mark=",")))

if (nrow(L) && "note" %in% names(L)) {
  cat("\n[corridor_tzero] per-combination outcome:\n"); print(L[, .N, by = note][order(-N)])
  D = L[note == "ok"]
  if (nrow(D)) {
    cat("\n[corridor_tzero] wall coverage by coordinate (of the combinations that got a base state):\n")
    print(D[, .(n = .N, off_wall = sum(is.finite(d_off)), on_wall = sum(is.finite(d_on)),
                both = sum(is.finite(w)), w_pos = sum(w > 0, na.rm = TRUE),
                med_w = round(median(w, na.rm = TRUE), 5)), by = coord][order(-both)])
    cat("\n[corridor_tzero] mono_B: does an ON wall exist off the Eh0 line?\n")
    print(D[rescue_mode == "mono_B", .(n = .N, on_wall = sum(is.finite(d_on)),
                pct = round(100*mean(is.finite(d_on)), 1)), by = coord][order(-pct)])
  }
}
