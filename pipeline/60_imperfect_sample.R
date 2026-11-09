## ============================================================================
## pipeline/60_imperfect_sample.R -- draw the units of the imperfect-Treg sweep
## (Fig. 5, Extended Data Figs. 6 and 7).
## ----------------------------------------------------------------------------
## A UNIT is one parameter combination (kinetic row tau_id = (batch, row_ind),
## Kg, r_L_P, beta_ros) together with one initial condition Eh0 on the
## initial-condition line.
##
## ELIGIBILITY (the load-balance population, CHECK_load_balance.R):
##   both modes : largest per-probe deviation from the assigned fixed point
##                <= MAP_DEV_TOL (1e-2) in both the Tregs-OFF and Tregs-ON line
##                (dev_max_on, dev_max_off from selected_big.rds; EXCLUDE_DEV=1)
##   bistable   : rescue_gap >= GAP_MIN (5 sterile, 2 pathogenic_help), both
##                directions increasing, gap not censored
##   mono_B     : Tregs-OFF direction increasing, Eh0_off finite
## The initial condition is e_mid, the midpoint of the rescued interval:
##   Eh0 = (Eh0_on + Eh0_off) / 2 (bistable),  Eh0 = Eh0_off / 2 (mono_B).
##
## THE DRAW. Strata are (Kg, beta_ros, r_L_P) within each branch. Within each
## stratum the eligible rows (both modes pooled) are permuted (set.seed(SEED))
## and the first min(n, K_PER) are taken; then bistable rows not yet drawn are
## added until the stratum holds BISTABLE_FLOOR bistable units, where available.
## Units are numbered 1..N in (branch, stratum, rescue_mode) order.
##
## ONE BRANCH PER FILE: lib/skeleton_defs.R seeds PL0 = PLP0 = 0 on the sterile
## branch and PL0 = 500, PLP0 = 10 otherwise; every eligible sterile row has
## r_L_P = 0 and every eligible pathogenic_help row r_L_P > 0, so each branch's
## units are integrated under their own BRANCH.
##
## WAVES. One cluster array (400 tasks, 23:59 walltime) cannot hold every
## (unit, replicate) item, so units are split into waves of at most WAVE_ITEMS
## items; a branch is never mixed within a wave.
##
## TOPUP=1 reads the ledger, excludes rows already drawn, tops each stratum up
## to K_PER and continues the unit numbering, so unit ids remain stable.
##
## in : ./data/rescue_<branch>.rds, ./data/skeleton_<branch>.rds,
##      ./data/selected_big.rds          (upstream fixed-point sweep products)
## out: ./data/fig05_units_<branch>.rds  one per branch, carrying `wave`
##      ./data/fig05_units.csv           all units
##      ./data/fig05_sample_ledger.csv   every unit ever drawn
## RUN: Rscript pipeline/60_imperfect_sample.R
## Verified 2026-10-01: with the defaults this reproduces the stored units files
## exactly (977 units; all.equal TRUE for both branches).
## ============================================================================
SRC         = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages(library(data.table))
BRANCHES    = strsplit(Sys.getenv("BRANCHES", "sterile,pathogenic_help"), ",")[[1]]
K_PER       = as.integer(Sys.getenv("K_PER", "5"))     # units per stratum, capped
BISTABLE_FLOOR = as.integer(Sys.getenv("BISTABLE_FLOOR", "2"))  # per stratum, where available
R_REP       = as.integer(Sys.getenv("R_REP", "10"))    # replicates per cell, for sizing only
WAVE_ITEMS  = as.integer(Sys.getenv("WAVE_ITEMS", "3200"))
LEDGER      = Sys.getenv("LEDGER", "./data/fig05_sample_ledger.csv")
TOPUP       = as.integer(Sys.getenv("TOPUP", "0"))
SEED        = as.integer(Sys.getenv("SEED", "1"))
EXCLUDE_DEV = !identical(Sys.getenv("EXCLUDE_DEV", "1"), "0")
MAP_DEV_TOL = as.numeric(Sys.getenv("MAP_DEV_TOL", "1e-2"))
CK = c("tau_id", "Kg", "beta_ros", "r_L_P")

BG = as.data.table(readRDS("./data/selected_big.rds"))[
  , .(dev_max_on = max(dev_max_on, na.rm = TRUE),
      dev_max_off = max(dev_max_off, na.rm = TRUE)), by = CK]

eligible = function(br) {
  GAP_MIN = as.numeric(Sys.getenv("GAP_MIN", if (br == "sterile") "5" else "2"))
  R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", br)))
  S = as.data.table(readRDS(sprintf("./data/skeleton_%s.rds", br)))
  R = merge(R, unique(S[, .(tau_id, Kg, beta_ros, r_L_P, batch, row_ind)]), by = CK)
  if (EXCLUDE_DEV) {
    n0 = nrow(R); R = merge(R, BG, by = CK, all.x = TRUE)
    R = R[is.finite(dev_max_on) & dev_max_on <= MAP_DEV_TOL &
          is.finite(dev_max_off) & dev_max_off <= MAP_DEV_TOL]
    cat(sprintf("[sample] %s EXCLUDE_DEV: %d of %d rows kept\n", br, nrow(R), n0))
  }
  BI = R[rescue_mode == "bistable" & rescue_gap >= GAP_MIN &
         dir_off == "increasing" & dir_on == "increasing" & !gap_censored]
  MB = R[rescue_mode == "mono_B" & dir_off == "increasing" & is.finite(Eh0_off)]
  U = rbindlist(list(BI, MB), fill = TRUE)
  if (!nrow(U)) return(U)
  ## the injury e_mid: midpoint of the rescued interval (Eh0_off / 2 for mono_B)
  U[, Eh0 := fifelse(rescue_mode == "bistable", (Eh0_on + Eh0_off) / 2, Eh0_off / 2)]
  U[, branch := br]
  U[is.finite(Eh0)]
}
E = rbindlist(lapply(BRANCHES, eligible), fill = TRUE)
if (!nrow(E)) stop("no eligible rows on ", paste(BRANCHES, collapse = ", "))
E[, stratum := paste(Kg, beta_ros, r_L_P, sep = "_")]
cat("\n[sample] eligible pool:\n")
print(E[, .(rows = .N, strata = uniqueN(stratum), med_Eh0 = round(median(Eh0), 1)),
        by = .(branch, rescue_mode)][order(branch, rescue_mode)])

## ---- the draw: K_PER per stratum, capped, honouring the ledger --------------
set.seed(SEED)
PREV = if (TOPUP == 1L && file.exists(LEDGER)) fread(LEDGER) else NULL
if (!is.null(PREV)) {
  E[, .drawn := paste(branch, tau_id, Kg, beta_ros, r_L_P) %in%
                paste(PREV$branch, PREV$tau_id, PREV$Kg, PREV$beta_ros, PREV$r_L_P)]
  cat(sprintf("\n[sample] TOPUP: ledger holds %d units; %d eligible rows already drawn\n",
              nrow(PREV), E[.drawn == TRUE, .N]))
} else {
  E[, .drawn := FALSE]
  if (TOPUP == 1L) cat("\n[sample] TOPUP requested but no ledger found -- drawing fresh\n")
}
## already-drawn units count toward the stratum's depth, so a top-up adds only
## the shortfall rather than K_PER more on top of what is there
HAVE = if (!is.null(PREV)) PREV[, .(have = .N), by = .(branch, stratum)] else
       E[0, .(branch, stratum)][, have := integer()]
NEW = E[.drawn == FALSE][, .SD[sample(.N)], by = .(branch, stratum)]
NEW = merge(NEW, HAVE, by = c("branch", "stratum"), all.x = TRUE)
NEW[is.na(have), have := 0L]
NEW = NEW[, .SD[seq_len(max(0L, min(.N, K_PER - have[1])))], by = .(branch, stratum)]

## the bistable floor, on top of the pooled draw and only where the stratum has
## them; already-drawn bistable units (pooled draw or ledger) count toward it
if (BISTABLE_FLOOR > 0L) {
  HAVE_BI = rbind(
    NEW[rescue_mode == "bistable", .(n = .N), by = .(branch, stratum)],
    if (!is.null(PREV)) PREV[rescue_mode == "bistable", .(n = .N), by = .(branch, stratum)]
    else NEW[0, .(branch, stratum)][, n := integer()])[, .(have_bi = sum(n)), by = .(branch, stratum)]
  taken = paste(NEW$branch, NEW$tau_id, NEW$Kg, NEW$beta_ros, NEW$r_L_P)
  POOL = E[.drawn == FALSE & rescue_mode == "bistable"][
    !paste(branch, tau_id, Kg, beta_ros, r_L_P) %in% taken][, .SD[sample(.N)], by = .(branch, stratum)]
  POOL = merge(POOL, HAVE_BI, by = c("branch", "stratum"), all.x = TRUE)
  POOL[is.na(have_bi), have_bi := 0L]
  EXTRA = POOL[, .SD[seq_len(max(0L, min(.N, BISTABLE_FLOOR - have_bi[1])))], by = .(branch, stratum)]
  if (nrow(EXTRA)) {
    cat(sprintf("\n[sample] BISTABLE_FLOOR = %d adds %d units over %d strata (%s)\n",
                BISTABLE_FLOOR, nrow(EXTRA), EXTRA[, uniqueN(paste(branch, stratum))],
                paste(sprintf("%s %d", EXTRA[, .N, by = branch]$branch,
                              EXTRA[, .N, by = branch]$N), collapse = ", ")))
    EXTRA[, have_bi := NULL]; EXTRA[, have := 0L]
    NEW = rbind(NEW, EXTRA, fill = TRUE)
  }
}
if (!nrow(NEW)) stop("nothing left to draw at K_PER = ", K_PER)
U = NEW
## ---- the columns 61_imperfect_sweep.R stamps onto every raw row ----------------
## cmp_off / cmp_on / help are NOT recomputed here: the sweep's own composites
## are the ones the fixed-point sweep recorded, so they are read from selected_big and
## left NA where it has no row rather than being invented from the rescue table.
LO = as.data.table(readRDS("./data/selected_big.rds"))[
  , .(lo_off = lo_off[1], hi_off = hi_off[1], lo_on = lo_on[1], hi_on = hi_on[1],
      margin = if ("margin" %in% names(.SD)) margin[1] else NA_real_), by = CK]
U = merge(U, LO, by = CK, all.x = TRUE)
U[, cmp_off := fifelse(is.finite(lo_off), lo_off, hi_off)]
U[, cmp_on  := fifelse(is.finite(lo_on),  lo_on,  hi_on)]
U[, help    := cmp_off - cmp_on]
setorder(U, branch, stratum, rescue_mode)
U[, unit := .I + (if (!is.null(PREV)) max(PREV$unit) else 0L)]

## ---- waves: BRANCH cannot be mixed in one job, so sterile goes first --------
U[, items := R_REP]
U[, wave := NA_integer_]
w = 0L
for (br in c("sterile", "pathogenic_help")) {
  idx = U[branch == br, which = TRUE]
  if (!length(idx)) next
  nw = max(1L, ceiling(length(idx) * R_REP / WAVE_ITEMS))
  ## alternate down the stratum-ordered list so each wave spans the design
  U[idx, wave := w + (seq_along(idx) - 1L) %% nw + 1L]
  w = w + nw
}

KEEP = c("unit", "branch", "wave", "batch", "row_ind", "tau_id", "Kg", "r_L_P", "beta_ros",
         "Eh0", "Eh0_off", "Eh0_on", "rescue_mode", "rescue_gap", "margin",
         "cmp_off", "cmp_on", "help", "stratum")
U = U[, ..KEEP]
cat("\n[sample] drawn:\n")
print(U[, .(units = .N, strata = uniqueN(stratum), med_Eh0 = round(median(Eh0), 1),
            med_help = round(median(help, na.rm = TRUE), 1)),
        by = .(branch, rescue_mode)][order(branch, rescue_mode)])
## the realised mode mix (bistable share per branch)
cat(sprintf("\n[sample] draw: K_PER = %d per stratum (capped), BISTABLE_FLOOR = %d where available\n",
            K_PER, BISTABLE_FLOOR))
cat("[sample] realised rescue_mode mix:\n")
print(U[, .(units = .N, bistable = sum(rescue_mode == "bistable"),
            mono_B = sum(rescue_mode == "mono_B"),
            pct_bistable = round(100 * mean(rescue_mode == "bistable"))), by = branch][order(branch)])
cat("[sample] per (branch x mode):\n")
print(U[, .(units = .N, strata = uniqueN(stratum),
            per_stratum = round(.N / uniqueN(stratum), 1)),
        by = .(branch, rescue_mode)][order(branch, rescue_mode)])
if (min(U[, .N, by = .(branch, rescue_mode)]$N) < 20)
  cat("[sample] WARNING: a (branch x mode) cell has fewer than 20 units; this cell will be thin\n")
cat("\n[sample] waves (", R_REP, " replicates per cell, 400 tasks each):\n", sep = "")
print(U[, .(branch = branch[1], units = .N, items = .N * R_REP,
            items_per_task = round(.N * R_REP / 400, 1),
            est_h_per_task = round(.N * R_REP / 400 * 2.70, 1)), by = wave][order(wave)])

for (br in unique(U$branch)) {
  f = sprintf("./data/fig05_units_%s.rds", br)
  saveRDS(U[branch == br], f)
  cat(sprintf("[sample] wrote %s (%d units, %d strata, waves %s)\n", f,
              U[branch == br, .N], U[branch == br, uniqueN(stratum)],
              paste(sort(unique(U[branch == br, wave])), collapse = ",")))
}
## the ledger: append, so it is the record of everything ever simulated
LED = U[, .(unit, branch, wave, tau_id, Kg, beta_ros, r_L_P, rescue_mode, stratum,
            Eh0, K_PER, drawn_on = as.character(Sys.Date()))]
fwrite(LED, LEDGER, append = file.exists(LEDGER) && TOPUP == 1L)
cat(sprintf("[sample] ledger %s now holds %d units\n", LEDGER,
            if (file.exists(LEDGER)) nrow(fread(LEDGER)) else nrow(LED)))
fwrite(U, "./data/fig05_units.csv")
cat(sprintf("[sample] %d units total -> ./data/fig05_units.csv\n", nrow(U)))
