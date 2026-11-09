## ============================================================================
## pipeline/06b_effero_fg_extra_units.R -- extra pathogenic units for Figure 2 F/G
##
## WHY. Panels F and G use a stratum-balanced median: at each beta_ros (F) every
## r_L_P stratum has equal weight, and at each r_L_P (G) every beta_ros stratum.
## The production draw (pipeline/06_effero_knockout.R, N_PER = 1500, SEED = 1)
## leaves some (beta_ros, r_L_P) cells with only 3-8 units that have a
## low-burden basin, because the fraction with a basin falls with r_L_P. This
## script draws extra pathogenic units aimed at those cells.
##
## DESIGN. For each of the 54 pathogenic cells (9 beta_ros x 6 r_L_P) the
## fraction of drawn units with a low-burden basin at >= 1 Kg (spread of the
## 21 final composites > GAP_CUT = 20, as in figure_scripts/Fig_2.R) is estimated from the
## production run as p = (b + 1) / (n + 2). N_EXTRA draws are then allocated one
## at a time to the cell with the smallest expected count b + p * (draws so far),
## i.e. the allocation raises the smallest cells first. Within a cell, kinetic
## vectors are drawn uniformly without replacement from the rows of the five tau
## grids that the cell does not already contain, so a cell holds exchangeable
## random draws from the grid, old and new alike. The extra units are used only in
## panels F and G (figure_scripts/Fig_2.R, FG_EXTRA_RDS); every other panel and number of
## Figure 2 rests on the production draw alone, whose design they describe.
##
## in : ./data/effero_knockout.rds, ./data/effero_knockout_units.rds, tau grids
## out: ./data/effero_fg_extra_units.rds   (batch, row_ind, beta_ros, r_L_P, uid,
##                                          and the cell's p and planned draws)
## RUN: N_EXTRA=4000 SEED=2 Rscript pipeline/06b_effero_fg_extra_units.R   (batch 1)
##      N_EXTRA=2600 SEED=3 UID0=200000 PRIOR_UNITS=./data/effero_fg_extra_units.rds \
##        OUT=./data/effero_fg_extra2_units.rds Rscript pipeline/06b_effero_fg_extra_units.R   (batch 2)
##      then N_CORES=6 DEADLINE_H=10.5 CHUNK_SIZE=100 \
##        UNITS_FILE=./data/effero_fg_extra_units.rds OUT_STEM=./data/effero_fg_extra \
##        Rscript pipeline/06_effero_knockout.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
BRANCH = "sterile"
suppressMessages(library(data.table))
source("./lib/skeleton_defs.R")                       # tg: the five tau grids

N_EXTRA = as.integer(Sys.getenv("N_EXTRA", "3000"))
SEED    = as.integer(Sys.getenv("SEED", "2"))
GAP_CUT = 20
OUT     = Sys.getenv("OUT", "./data/effero_fg_extra_units.rds")
## PRIOR_UNITS: comma-separated unit tables of earlier extra batches. Their units
##   are excluded from the new draw, and their planned draws count towards each
##   cell's expected count (b + p * draws), so a new batch fills what is still thin.
PRIOR_UNITS = Sys.getenv("PRIOR_UNITS", "")
UID0        = as.integer(Sys.getenv("UID0", "100000"))   # new uids are UID0 + 1, 2, ...

KD = as.data.table(readRDS("./data/effero_knockout.rds"))
KD = KD[r_L_P > 0 & eff_ko == 1L & is.finite(comp)]
has_lower = function(v) length(v) >= 2L && (max(v) - min(v)) > GAP_CUT
B = KD[, .(b = has_lower(comp)), by = .(batch, row_ind, beta_ros, r_L_P, Kg)
       ][, .(b = any(b)), by = .(batch, row_ind, beta_ros, r_L_P)]
U0 = as.data.table(readRDS("./data/effero_knockout_units.rds"))[r_L_P > 0]
C = merge(U0[, .(n = .N), by = .(beta_ros, r_L_P)],
          B[, .(b = sum(b)), by = .(beta_ros, r_L_P)], by = c("beta_ros", "r_L_P"), all.x = TRUE)
C[is.na(b), b := 0L]
C[, p := (b + 1) / (n + 2)]
C[, draws := 0L]
PRIOR = if (nzchar(PRIOR_UNITS)) rbindlist(lapply(trimws(strsplit(PRIOR_UNITS, ",")[[1]]),
          function(f) as.data.table(readRDS(f))[, .(batch, row_ind, beta_ros, r_L_P, uid)])) else
        data.table(batch = integer(), row_ind = integer(), beta_ros = numeric(), r_L_P = numeric(), uid = integer())
if (nrow(PRIOR)) {
  pc = PRIOR[, .(prior = .N), by = .(beta_ros, r_L_P)]
  C = merge(C, pc, by = c("beta_ros", "r_L_P"), all.x = TRUE); C[is.na(prior), prior := 0L]
  C[, draws := prior]
  stopifnot(UID0 >= max(PRIOR$uid))
  cat(sprintf("[06b] %d prior extra units counted and excluded\n", nrow(PRIOR)))
} else C[, prior := 0L]
for (i in seq_len(N_EXTRA)) {                         # raise the smallest expected count
  j = C[, which.min(b + p * draws)]
  C[j, draws := draws + 1L]
}
C[, expected := b + p * draws]
C[, new := draws - prior]
cat(sprintf("[06b] %d extra pathogenic units over %d cells\n", N_EXTRA, nrow(C)))
cat(sprintf("      units with a basin per cell: now %d-%d; expected after %.1f-%.1f\n",
            min(C$b), max(C$b), min(C$expected), max(C$expected)))
print(dcast(C[, .(beta_ros, r_L_P, v = sprintf("%2d+%3d+%3d", b, prior, new))], beta_ros ~ r_L_P, value.var = "v"))

set.seed(SEED)
TAU = tg[is.finite(row_ind), .(batch, row_ind)]
setorder(C, r_L_P, beta_ros)
units = rbindlist(lapply(seq_len(nrow(C)), function(j) {
  cj = C[j]; if (cj$new == 0L) return(NULL)
  have = rbind(U0[beta_ros == cj$beta_ros & r_L_P == cj$r_L_P, .(batch, row_ind)],
               PRIOR[beta_ros == cj$beta_ros & r_L_P == cj$r_L_P, .(batch, row_ind)])
  pool = TAU[!have, on = .(batch, row_ind)]
  i = sample(nrow(pool), cj$new)
  data.table(batch = pool$batch[i], row_ind = pool$row_ind[i],
             beta_ros = cj$beta_ros, r_L_P = cj$r_L_P, cell_p = cj$p) }))
units[, uid := UID0 + .I]                             # distinct from production (<= 10,500) and prior batches
units = units[sample(.N)]                             # interleave cells across cores and chunks
stopifnot(!nrow(merge(units, rbind(U0[, .(batch, row_ind, beta_ros, r_L_P)], PRIOR[, .(batch, row_ind, beta_ros, r_L_P)]),
                       by = c("batch", "row_ind", "beta_ros", "r_L_P"))))
saveRDS(units, OUT)
fwrite(C, sub("\\.rds$", "_plan.csv", OUT))
cat(sprintf("[06b] wrote %s (%s units) and the per-cell plan\n", OUT, format(nrow(units), big.mark = ",")))
