## ============================================================================
## pipeline/checks/sweep_smoke_test.R -- does the sweep, run here, give the stored
## result?
##
## The full sweep is a cluster job (about 1.7 core-hours per kinetic vector, 6,000
## vectors). This script runs the SAME code path on a handful of parameter
## combinations that are in the gated set, and compares what comes out with the
## stored rows of data/selected_big.rds:
##   pipeline/01_sweep.R     in cell-list mode (CELLS), all other settings default
##   pipeline/02_assemble.R  on the chunk it wrote
## It takes a few minutes and writes only to a temporary folder and to
## ./data/smoke_fp_batch_<b>*.rds, which it removes at the end.
##
## in : ./data/selected_big.rds, ./data/tau_grid_v6_<b>.rds
## out: a table on the console: per combination, the largest difference in any
##      numeric column, and whether the class labels agree
## RUN: Rscript pipeline/checks/sweep_smoke_test.R           (from the repository root)
##      N_STERILE=2 N_PATHO=2 BATCH=1 Rscript pipeline/checks/sweep_smoke_test.R
## ============================================================================
suppressMessages(library(data.table))
SRC = Sys.getenv("SRC", "."); setwd(SRC)
BATCH     = as.integer(Sys.getenv("BATCH", "1"))
N_STERILE = as.integer(Sys.getenv("N_STERILE", "1"))
N_PATHO   = as.integer(Sys.getenv("N_PATHO", "1"))
CK = c("tau_id", "Kg", "beta_ros", "r_L_P")

S = as.data.table(readRDS("./data/selected_big.rds"))[batch == BATCH]
## one ROS level for all cells (cell-list mode runs every listed cell at every level given)
b0 = S[, .N, by = beta_ros][order(-N)]$beta_ros[1]
set.seed(1)
pick = rbind(S[beta_ros == b0 & r_L_P == 0][sample.int(.N, N_STERILE)],
             S[beta_ros == b0 & r_L_P > 0 & rescue_gap_kind == "measured"][sample.int(.N, N_PATHO)])
cat(sprintf("[smoke] %d combinations of grid %d at beta_ros = %s:\n", nrow(pick), BATCH, format(b0, digits = 17)))
print(pick[, .(tau_id, Kg, beta_ros, r_L_P, rescue_mode, treg_effect, Eh0_off, Eh0_on)])

tmp = tempfile("smoke_"); dir.create(tmp)
cells = file.path(tmp, "cells.rds")
saveRDS(pick[, .(batch, row_ind, Kg, r_L_P)], cells)
env = c(sprintf("CELLS=%s", cells), sprintf("BETA_ROS_VALS=%s", sprintf("%.17g", b0)),
        sprintf("OUT_PREFIX=%s/fp_batch_", tmp), sprintf("SRC=%s", getwd()))
t0 = Sys.time()
st = system2("Rscript", c("pipeline/01_sweep.R", "1", "1", BATCH), env = env, stdout = file.path(tmp, "sweep.log"), stderr = file.path(tmp, "sweep.log"))
if (st != 0) { cat(readLines(file.path(tmp, "sweep.log")), sep = "\n"); stop("pipeline/01_sweep.R failed") }
cat(sprintf("[smoke] sweep done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
st = system2("Rscript", "pipeline/02_assemble.R",
             env = c(sprintf("DATA_DIR=%s/fp_batch_%d", tmp, BATCH), sprintf("BATCH=%d", BATCH), "OUT_STEM=smoke_fp",
                     sprintf("SRC=%s", getwd())), stdout = file.path(tmp, "assemble.log"), stderr = file.path(tmp, "assemble.log"))
if (st != 0) { cat(readLines(file.path(tmp, "assemble.log")), sep = "\n"); stop("pipeline/02_assemble.R failed") }
out = sprintf("./data/smoke_fp_batch_%d.rds", BATCH)
N = as.data.table(readRDS(out))
invisible(file.remove(Sys.glob(sprintf("./data/smoke_fp_batch_%d*.rds", BATCH))))

## the ROS level went through text on its way into the sweep, so match it by nearness
N[, beta_ros := b0]
M = merge(pick, N, by = CK, suffixes = c(".stored", ".new"))
if (nrow(M) != nrow(pick)) stop("the sweep returned ", nrow(M), " of the ", nrow(pick), " combinations")
both = intersect(setdiff(names(pick), CK), names(N))
num = both[vapply(both, function(n) is.numeric(pick[[n]]) && !is.list(pick[[n]]), NA)]
num = setdiff(num, c("batch", "row_ind"))
chr = both[vapply(both, function(n) is.character(pick[[n]]), NA)]
res = M[, {
  d = vapply(num, function(n) { a = get(paste0(n, ".stored")); b = get(paste0(n, ".new"))
    if (is.na(a) && is.na(b)) 0 else if (is.na(a) != is.na(b)) Inf else abs(a - b) / max(1, abs(a)) }, 0)
  same = vapply(chr, function(n) identical(get(paste0(n, ".stored")), get(paste0(n, ".new"))), NA)
  .(numeric_columns = length(num), largest_difference = max(d), in_column = num[which.max(d)],
    label_columns = length(chr), labels_agree = all(same),
    labels_differing = paste(chr[!same], collapse = ","))
}, by = CK]
cat("\n[smoke] re-run against the stored rows (difference = |stored - new| / max(1, |stored|)):\n")
print(res)
cat(if (all(res$labels_agree) && all(res$largest_difference < 1e-6))
  "\n[smoke] PASS: the sweep reproduces the stored rows.\n" else
  "\n[smoke] CHECK: something differs by more than 1e-6 or a label changed; see the table.\n")
