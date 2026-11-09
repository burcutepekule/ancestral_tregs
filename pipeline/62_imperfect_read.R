## ============================================================================
## pipeline/62_imperfect_read.R -- merge the imperfect-Treg array's csv, drop
## what did not solve, keep the units that carry a complete grid.
## ----------------------------------------------------------------------------
## RUN FILTER (per run): dropped if solver_ok is FALSE, Eh_mean or PLP_mean is
## not finite, a tail mean of PLP or CLP exceeds BUG_MAX (1e4), or a tail mean
## of Eh, Ei = K_E - Eh, M0, M1 or Tr exceeds CELL_MAX (101).
## UNIT FILTER: replicates are averaged per (unit, phi, eta); a unit is kept
## only if all PHI_N x ETA_N cells have at least one surviving replicate, and
## (through the merge with the untreated reference) only if its untreated run
## (theta_treg = 0, rep = 0) survived.
## The across-replicate spread of the composite res is printed before the
## replicates are collapsed, and its sd and se are carried into the cells file.
##
## in : $RAWDIR/**/fig05_raw_t*.csv          (cluster scratch)
## out: fig05_raw_all.rds   every surviving run, replicates intact
##      fig05_cells.rds     one row per (unit, phi, eta), replicates averaged
##      (under $SCRATCH_BASE on the cluster, else ./data; copy both to ./data)
## RUN: Rscript pipeline/62_imperfect_read.R     (on the cluster)
## ============================================================================
SRC    = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages(library(data.table))
## SYNTHETIC DATA LIVES ON SCRATCH, NOT IN HOME. The array writes the raw csv to
## $SCRATCH_BASE/fig05_raw/<branch>/ and the merged products go back beside them;
## only the scripts and the units files stay in the home directory. On a laptop
## neither path exists, so both fall back to the working directory.
SCRATCH_BASE = Sys.getenv("SCRATCH_BASE", "~/scratch")
OUT_BASE = Sys.getenv("OUT_BASE", if (dir.exists(SCRATCH_BASE)) SCRATCH_BASE else "./data")
RAWDIR = Sys.getenv("RAWDIR", file.path(OUT_BASE, "fig05_raw"))
## RAW_SUBDIRS merges only the named wave directories, e.g. "sterile_w1,sterile_w2".
## Needed while other waves are still running: their csv are being appended to as
## tasks finish, so a full merge would read half-written files and, worse, would
## look complete. Incomplete units are dropped by the full-grid test below anyway,
## but reading only the finished waves is the honest way to look early.
## OUT_TAG suffixes the outputs so a partial merge cannot overwrite the real one.
RAW_SUBDIRS = Sys.getenv("RAW_SUBDIRS", "")
OUT_TAG     = Sys.getenv("OUT_TAG", "")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[read] raw from %s | merged products to %s\n", RAWDIR, OUT_BASE))
KE       = as.numeric(Sys.getenv("K_E", "100"))
CELL_MAX = as.numeric(Sys.getenv("CELL_MAX", "101"))
BUG_MAX  = as.numeric(Sys.getenv("BUG_MAX",  "1e4"))

## RECURSIVE, because the two branches write to fig05_raw/<branch>/ and their
## task ids would collide in one flat directory. Unit ids are globally unique
## across branches (60_imperfect_sample.R numbers the combined table), so the files
## merge without a branch key.
fs = if (nzchar(RAW_SUBDIRS)) {
  dirs = file.path(RAWDIR, trimws(strsplit(RAW_SUBDIRS, ",")[[1]]))
  miss = dirs[!dir.exists(dirs)]
  if (length(miss)) stop("RAW_SUBDIRS names directories that do not exist: ",
                         paste(basename(miss), collapse = ", "))
  cat(sprintf("[read] merging only: %s\n", paste(basename(dirs), collapse = ", ")))
  unlist(lapply(dirs, list.files, pattern = "^fig05_raw_t[0-9]+\\.csv$",
                full.names = TRUE, recursive = TRUE))
} else list.files(RAWDIR, "^fig05_raw_t[0-9]+\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(fs)) stop("no fig05_raw_t*.csv under ", RAWDIR)
## A FLAT FILE UNDER RAWDIR IS FROM A RETIRED SAMPLE. The rebuilt arrays write
## to fig05_raw/<branch>/, so anything sitting directly in RAWDIR predates them.
## Merging it would combine two samples whose unit ids mean different things.
.flat = list.files(RAWDIR, "^fig05_raw_t[0-9]+\\.csv$", full.names = TRUE, recursive = FALSE)
if (length(.flat))
  stop(length(.flat), " raw csv sit directly in ", RAWDIR, " rather than in a per-branch\n",
       "  subdirectory. These are from a retired sample and their unit ids collide with the\n",
       "  current one. Remove them (rm -rf ", RAWDIR, ") and rerun the arrays.")
B = rbindlist(lapply(fs, fread), fill = TRUE)
## NO RUN MAY BE PRESENT TWICE, CHECKED BEFORE ANYTHING JOINS ON IT. A partial
## run left on disk and merged with a rerun gives duplicate
## (unit, theta_treg, phi, eta, rep) rows; the cell mean would then average a run
## against a stale copy of itself, and the joins below fail with a cartesian
## blow-up that says nothing about the cause. The id guards above do not catch
## this -- leftovers sit in a valid per-branch directory with valid unit ids.
## Nothing here is a fixed count: it compares rows to distinct keys, so it holds
## whatever the grid, the replicate count or the sample size are.
.key = c("unit", "theta_treg", "phi_instability", "treg_eff", "rep")
.dup = B[, .N, by = .key][N > 1L]
if (nrow(.dup))
  stop(nrow(.dup), " of ", nrow(B), " runs share a (unit, theta_treg, phi, eta, rep) key with\n",
       "  another -- output from two different runs is mixed under ", RAWDIR, ".\n",
       "  Clear it and rerun the arrays; a mean over a run and a stale copy of itself is not a mean.")
cat(sprintf("[read] %s runs from %d files, all keys distinct\n",
            format(nrow(B), big.mark = ","), length(fs)))
cat(sprintf("[read] %d files | %d rows | %d units\n", length(fs), nrow(B), uniqueN(B$unit)))

## ---- bad runs --------------------------------------------------------------
CAP_COLS = intersect(c("Eh_mean", "Ei_mean", "M0_mean", "M1_mean", "Tr_mean"), names(B))
BUG_COLS = intersect(c("PLP_mean", "CLP_mean"), names(B))
B[, bad := solver_ok == FALSE |
      !is.finite(Eh_mean) | !is.finite(PLP_mean) |
      Reduce(`|`, lapply(.SD[, BUG_COLS, with = FALSE], function(z) !is.finite(z) | z > BUG_MAX)) |
      Reduce(`|`, lapply(.SD[, CAP_COLS, with = FALSE], function(z) !is.finite(z) | z > CELL_MAX))]
cat(sprintf("[read] bad runs: %d of %d (%.2f%%)\n", B[bad == TRUE, .N], nrow(B),
            100 * B[bad == TRUE, .N] / nrow(B)))
if (B[bad == TRUE, .N]) print(B[bad == TRUE, .N, by = .(theta_treg, treg_eff)][order(-N)][1:min(5, .N)])
B = B[bad == FALSE][, bad := NULL]

## ---- the untreated reference, one per unit ---------------------------------
OFFD = B[theta_treg == 0]
ONE  = B[theta_treg == 1]
cat(sprintf("[read] untreated rows: %d for %d units\n", nrow(OFFD), uniqueN(OFFD$unit)))

## ---- how far apart the replicates of a cell are ----------------------------
## REPORTED BEFORE IT IS COLLAPSED. The replicate count is justified only if
## this number is small; if it is not, the figures need more replicates
## and this is where that shows.
ONE[, res := sqrt((Ei_mean / KE)^2 + (PLP_mean / 100)^2)]
SP = ONE[, .(n = .N, m = mean(res), s = sd(res)), by = .(unit, phi_instability, treg_eff)][n > 1]
## THE ACROSS-REPLICATE SD IS CARRIED, NOT JUST REPORTED. It used to be computed
## here for the console and then dropped, which left every downstream error bar
## on the cell mean calibrated from a docstring instead of measured. With it in
## CELL, the standard error of a cell mean is s / sqrt(n) directly.
SPK = ONE[, .(res_sd = sd(res), res_n = .N), by = .(unit, phi_instability, treg_eff)]
cat(sprintf("[read] replicate CV within a cell: median %.4f, 95th %.4f, max %.4f over %d cells\n",
            median(SP$s / SP$m, na.rm = TRUE), quantile(SP$s / SP$m, .95, na.rm = TRUE),
            max(SP$s / SP$m, na.rm = TRUE), nrow(SP)))

## ---- average the replicates, then require a complete grid ------------------
NUM = names(ONE)[sapply(ONE, is.numeric)]
## phi and eta are the grouping, and the unit's identity is re-attached from ID
## below -- averaging Kg over replicates and then merging Kg back would leave
## data.table holding Kg.x and Kg.y and nothing downstream would find either
NUM = setdiff(NUM, c("unit", "rep", "theta_treg", "phi_instability", "treg_eff",
                     "batch", "row_ind", "tau_id", "Kg", "r_L_P", "beta_ros",
                     "Eh0", "margin", "cmp_off", "cmp_on", "help"))
CELL = ONE[, lapply(.SD, mean), by = .(unit, phi_instability, treg_eff), .SDcols = NUM]
CELL = merge(CELL, SPK, by = c("unit", "phi_instability", "treg_eff"), all.x = TRUE)
CELL[, res_se := res_sd / sqrt(pmax(res_n, 1L))]
NCELL = uniqueN(ONE[, .(phi_instability, treg_eff)])
full = CELL[, .N, by = unit][N == NCELL, unit]
cat(sprintf("[read] complete %d-cell grid: %d of %d units\n", NCELL, length(full), uniqueN(CELL$unit)))
CELL = CELL[unit %in% full]
OFFD = OFFD[unit %in% full]
## the untreated reference travels with every cell
OFFC = OFFD[, .(unit, res_off = sqrt(((KE - Eh_mean) / KE)^2 + (PLP_mean / 100)^2),
                Ei_mean_off = KE - Eh_mean, PLP_mean_off = PLP_mean,
                Ei_var_off = Eh_var, Eh_mean_off = Eh_mean)]
CELL = merge(CELL, OFFC, by = "unit")
## and the unit's identity
ID = unique(B[, .(unit, batch, row_ind, tau_id, Kg, r_L_P, beta_ros, Eh0,
                  rescue_mode, margin, cmp_off, cmp_on, help)])
CELL = merge(CELL, ID, by = "unit")

## EVERY UNIT READ MUST BE A UNIT WE DREW. This is the generic form of the check
## above: it catches a stale csv anywhere under RAWDIR, not just a flat one, and
## it counts nothing -- it compares two id sets.
.uf = list.files(Sys.getenv("UNITS_DIR", "./data"), "^fig05_units_.*\\.rds$", full.names = TRUE)
if (length(.uf)) {
  .known = unique(unlist(lapply(.uf, function(f) as.data.table(readRDS(f))$unit)))
  .extra = setdiff(unique(B$unit), .known)
  if (length(.extra))
    stop(length(.extra), " unit ids in the raw csv are not in any fig05_units_*.rds ",
         "(e.g. ", paste(head(sort(.extra), 5), collapse = ", "), "). The raw directory ",
         "holds output from a sample that is no longer the current one -- clear it and rerun.")
  cat(sprintf("[read] all %d raw unit ids are in the current sample\n", uniqueN(B$unit)))
}
saveRDS(B,    file.path(OUT_BASE, sprintf("fig05_raw_all%s.rds", OUT_TAG)))
saveRDS(CELL, file.path(OUT_BASE, sprintf("fig05_cells%s.rds", OUT_TAG)))
cat(sprintf("[read] wrote %s/fig05_raw_all%s.rds (%d rows) and %s/fig05_cells%s.rds (%d rows, %d units)\n",
            OUT_BASE, OUT_TAG,
            nrow(B), OUT_BASE, OUT_TAG, nrow(CELL), uniqueN(CELL$unit)))
