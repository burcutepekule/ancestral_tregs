## ============================================================================
## pipeline/06_effero_knockout.R -- the efferocytosis-knockout experiment
## behind Figure 2. Local, parallel (forked with mclapply); not a cluster job.
##
## WHAT IT DOES
## The main sweep has only one toggle, theta_treg; efferocytosis is never switched
## off there. This experiment adds the missing arm. Tregs are OFF throughout
## (theta_treg = 0).
##
## One UNIT is a (kinetic vector, beta_ros, r_L_P) triple: a kinetic vector is a
## row (batch, row_ind) of the five tau grids ./data/tau_grid_v6_<n>.rds. Each
## unit is integrated under 8 ARMS -- efferocytosis intact at each Kg in KG_GRID
## (7 arms), and one knockout arm -- from each of 21 initial conditions (the EH0_FAN
## line), so one unit = 8 x 21 = 168 integrations.
##
## THE KNOCKOUT is eff_ko = 0 in rhs_app (lib/skeleton_defs.R): it multiplies the
## clearance rate g(L), and only g. It is NOT a0 = a_floor = 0: K_En is defined as
## (K_E*tau_P_B + mu_Ea*K_Ea)/(a0*K_M), so a0 = 0 sends K_En to infinity, f_En to
## 0, and also switches off DAMP production from necrotic debris. K_Ea and K_En
## keep their intact values under knockout. The knockout arm does not depend on
## Kg (g is identically 0); Kg is set to KG_GRID[1] as a placeholder and the
## stored Kg is NA.
##
## SAMPLING. For each r_L_P level, N_PER kinetic vectors are drawn without
## replacement (N_PER < 10,000) from all rows of the five grids, and each is paired
## with a beta_ros level drawn uniformly with replacement. Equal weight per r_L_P
## stratum. Sterile = r_L_P == 0 (no pathogen in the initial condition either).
## The unit order is then shuffled (only for load balance across cores) and cut
## into chunks of CHUNK_SIZE that are written as they finish (RESUME skips chunks
## already on disk). Both the draw and the chunk assignment depend on N_PER and
## SEED: never rerun with a different N_PER on an existing OUT_STEM.
##
## GRIDS must match the main sweep: Kg 15..45 by 5, beta_ros 0.0015..0.0095 by
## 0.001, r_L_P 0..0.48 by 0.08.
##
## SOLVER. deSolve::ode, method lsodar with the extinction root and event of
## lib/extinction_floor.R (PL, PLP set to 0 when max(PL, PLP) falls to 1e-11),
## atol 1e-14, rtol 1e-8, hmax 5, maxsteps 1e5, from t = 0 to TMAX = 12000. The
## stored final state is the state at t = TMAX (not a tail mean). A run that
## errors or stops before TMAX is dropped (no row written).
##
## in : ./data/tau_grid_v6_<n>.rds, ./lib/skeleton_defs.R
## out: ./data/effero_knockout.rds          (the assembled table, read by figure_scripts/Fig_2.R)
##      ./data/effero_knockout_units.rds    (the sampled unit table, before integration)
##      ./data/effero_knockout_partNNN.rds  (chunk checkpoints)
## RUN (the production run, 10,500 units, 1,764,000 integrations; measured
##      19.5 h on 10 cores in the original run):
##      N_PER=1500 SEED=1 N_CORES=6 Rscript pipeline/06_effero_knockout.R
## CHECKS:
##      UNITS_ONLY=true  ... -> writes only the unit table and stops (seconds)
##      UIDS=12,345 OUT_STEM=./data/effero_knockout_check ... -> integrates only
##      those units of the same draw (to spot-check stored final states)
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
## skeleton_defs reads BRANCH to set PL0/PLP0 for its own seed(). This script
## builds its initial states explicitly, per unit, so that setting is unused here.
BRANCH = "sterile"
suppressMessages({library(data.table); library(deSolve); library(parallel)})
source("./lib/skeleton_defs.R")

KG_GRID       = c(15, 20, 25, 30, 35, 40, 45)
BETA_ROS_VALS = seq(0.0015, 0.0095, by = 0.001)
RLP_VALS      = seq(0.00, 0.48, 0.08)
PLP0_HERE     = as.numeric(Sys.getenv("PLP0", "10"))
EH0_FAN       = as.numeric(strsplit(Sys.getenv("EH0_FAN", ""), ",")[[1]])
if (!length(EH0_FAN)) EH0_FAN = seq(0, 99, length.out = 21)

TMAX     = as.numeric(Sys.getenv("TMAX", "12000"))
HMAX     = as.numeric(Sys.getenv("HMAX", "5"))
N_PER    = as.integer(Sys.getenv("N_PER", "1500"))
N_CORES  = as.integer(Sys.getenv("N_CORES", "6"))
SEED     = as.integer(Sys.getenv("SEED", "1"))
OUT_STEM = Sys.getenv("OUT_STEM", "./data/effero_knockout")
CHUNK_SIZE = as.integer(Sys.getenv("CHUNK_SIZE", "250"))
RESUME     = !identical(tolower(Sys.getenv("RESUME", "true")), "false")
UNITS_ONLY = identical(tolower(Sys.getenv("UNITS_ONLY", "false")), "true")
UIDS       = as.integer(strsplit(Sys.getenv("UIDS", ""), ",")[[1]])
## UNITS_FILE: integrate a unit table written elsewhere (pipeline/06b_effero_fg_extra_units.R)
##   instead of drawing one; N_PER and SEED are then unused for the draw.
## DEADLINE_H: wall-clock limit in hours; no chunk is started that the mean chunk
##   time so far (x 1.5) says would end after it. Finished chunks stay on disk.
UNITS_FILE = Sys.getenv("UNITS_FILE", "")
DEADLINE_H = as.numeric(Sys.getenv("DEADLINE_H", "Inf"))

## ---- the unit draw -----------------------------------------------------------
## TAU is ordered by (batch, row_ind) because tg is keyed on it; the draw is
## therefore reproducible from SEED and N_PER alone.
if (nzchar(UNITS_FILE)) {
units = as.data.table(readRDS(UNITS_FILE))
cat(sprintf("[effero-KO] read %s (%s units)\n", UNITS_FILE, format(nrow(units), big.mark = ",")))
} else {
set.seed(SEED)
TAU = tg[is.finite(row_ind), .(batch, row_ind)]
units = rbindlist(lapply(RLP_VALS, function(rlp) {
  i = sample(nrow(TAU), N_PER, replace = N_PER > nrow(TAU))
  data.table(batch = TAU$batch[i], row_ind = TAU$row_ind[i],
             beta_ros = sample(BETA_ROS_VALS, N_PER, replace = TRUE),
             r_L_P = rlp) }))
units[, uid := .I]
## shuffle so the strata are interleaved across cores (load balance only)
units = units[sample(.N)]
saveRDS(units, sprintf("%s_units.rds", OUT_STEM))
cat(sprintf("[effero-KO] wrote %s_units.rds (%s units)\n", OUT_STEM, format(nrow(units), big.mark = ",")))
}
if (UNITS_ONLY) quit(save = "no")
if (length(UIDS)) {
  if (OUT_STEM == "./data/effero_knockout") stop("UIDS needs its own OUT_STEM")
  units = units[uid %in% UIDS]
}
units[, chunk := ((seq_len(.N) - 1L) %/% CHUNK_SIZE) + 1L]
NARM = length(KG_GRID) + 1L
cat(sprintf("[effero-KO] %s units (%s) | arms %d intact Kg + 1 KO | Eh0 fan %d | Tregs OFF\n",
            format(nrow(units), big.mark=","),
            if (nzchar(UNITS_FILE)) UNITS_FILE else sprintf("%d strata x %s", length(RLP_VALS), format(N_PER, big.mark=",")),
            length(KG_GRID), length(EH0_FAN)))
cat(sprintf("            Kg %s | beta_ros %.4f..%.4f | PLP0 %g | TMAX %g\n",
            paste(KG_GRID, collapse=","), min(BETA_ROS_VALS), max(BETA_ROS_VALS), PLP0_HERE, TMAX))
cat(sprintf("            %s integrations total on %d core(s)\n",
            format(nrow(units) * NARM * length(EH0_FAN), big.mark=","), N_CORES))

## ---- initial condition: the sweep's line ----------------------------------------
## Eh = e, En = K_E - e, M0 = K_M, T0 = K_T, everything else 0; pathogenic units
## also start with PL = 500 and PLP = PLP0_HERE (10).
build_state = function(Eh0, pathogenic) {
  y = setNames(numeric(length(STATE_ORDER)), STATE_ORDER)
  y["Eh"] = Eh0; y["En"] = K_E - Eh0; y["M0"] = M_TOT; y["T0"] = T_TOT
  if (pathogenic) { y["PL"] = 500; y["PLP"] = PLP0_HERE }
  y }

one_unit = function(i) {
  u = units[i]
  trow = tg[.(u$batch, u$row_ind), nomatch = NULL]
  if (!nrow(trow)) return(NULL)
  trow = trow[1]
  patho = u$r_L_P > 0
  mkp = function(KG, ko) {
    p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]]
    p$Kg = KG; p$beta_ros = u$beta_ros; p$r_L_P = u$r_L_P; p$r_LP_P = u$r_L_P
    p$theta_treg = 0L; p$eff_ko = ko; p }
  arms = c(lapply(KG_GRID, function(k) list(arm = sprintf("Kg%02d", k), Kg = k, ko = 1L)),
           list(list(arm = "KO", Kg = NA_integer_, ko = 0L)))
  rbindlist(lapply(arms, function(a) {
    P = mkp(if (is.na(a$Kg)) KG_GRID[1] else a$Kg, a$ko)
    rbindlist(lapply(EH0_FAN, function(e0) {
      o = tryCatch(ode(build_state(e0, patho), c(0, TMAX), rhs_app, P,
                       method = "lsodar", rootfun = rootfun_ext, events = EVENTS_EXT,
                       maxsteps = 1e5, hmax = HMAX, atol = ATOL, rtol = RTOL),
                   error = function(e) NULL)
      ok = !is.null(o) && nrow(o) >= 2 && o[nrow(o), 1] >= TMAX * (1 - 1e-9)
      if (!ok) return(NULL)
      z = setNames(as.numeric(o[nrow(o), STATE_ORDER]), STATE_ORDER)
      data.table(batch = u$batch, row_ind = u$row_ind, beta_ros = u$beta_ros,
                 r_L_P = u$r_L_P, pathogenic = patho, uid = u$uid,
                 Eh0 = e0, arm = a$arm, Kg = a$Kg, eff_ko = a$ko,
                 comp = cmpf(z[["Eh"]], z[["PLP"]]),
                 Eh = z[["Eh"]], PLP = z[["PLP"]], En = z[["En"]],
                 M1 = z[["M1"]], R = z[["R"]], Ea = z[["Ea"]])
    }), fill = TRUE) }), fill = TRUE) }

t0 = Sys.time()
NCH = max(units$chunk)
N_DONE = 0L
for (k in seq_len(NCH)) {
  fk = sprintf("%s_part%03d.rds", OUT_STEM, k)
  if (RESUME && file.exists(fk)) { cat(sprintf("[effero-KO] chunk %d/%d already on disk, skipping\n", k, NCH)); next }
  el_h = as.numeric(difftime(Sys.time(), t0, units = "hours"))
  if (N_DONE > 0L && el_h + 1.5 * el_h / N_DONE > DEADLINE_H) {
    cat(sprintf("[effero-KO] DEADLINE_H = %g: stopping before chunk %d/%d (%.2f h elapsed)\n", DEADLINE_H, k, NCH, el_h))
    break }
  idx = which(units$chunk == k)
  tk = Sys.time()
  rk = rbindlist(if (N_CORES > 1)
    mclapply(idx, function(i) tryCatch(one_unit(i), error = function(e) NULL),
             mc.cores = N_CORES, mc.preschedule = FALSE) else
    lapply(idx, one_unit), fill = TRUE)
  saveRDS(rk, fk)
  N_DONE = N_DONE + 1L
  el = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[effero-KO] chunk %2d/%d | %s units | %s rows | %.1f min this chunk | %.1f min elapsed | eta %.1f min\n",
              k, NCH, format(length(idx), big.mark=","), format(nrow(rk), big.mark=","),
              as.numeric(difftime(Sys.time(), tk, units = "mins")), el, el/k * (NCH - k)))
}
RES = rbindlist(lapply(sprintf("%s_part%03d.rds", OUT_STEM, seq_len(NCH)),
                       function(f) if (file.exists(f)) readRDS(f) else NULL), fill = TRUE)
cat(sprintf("\n[effero-KO] %s rows in %.1f min\n", format(nrow(RES), big.mark=","),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (!nrow(RES)) stop("nothing integrated")

cat("\n[effero-KO] rows per arm:\n"); print(RES[, .N, by = arm][order(arm)])
COMPLETE = RES[, .(n = .N), by = uid][n == NARM * length(EH0_FAN)]
cat(sprintf("\n[effero-KO] units returning every arm and every Eh0: %s of %s\n",
            format(nrow(COMPLETE), big.mark=","), format(nrow(units), big.mark=",")))
BAD = RES[!is.finite(comp)]
cat(sprintf("[effero-KO] non-finite final states: %s of %s rows\n",
            format(nrow(BAD), big.mark=","), format(nrow(RES), big.mark=",")))
if (nrow(BAD)) print(BAD[, .N, by = .(arm, pathogenic)][order(arm, pathogenic)])

saveRDS(RES, sprintf("%s.rds", OUT_STEM))
cat(sprintf("\n[effero-KO] wrote %s.rds\n", OUT_STEM))
