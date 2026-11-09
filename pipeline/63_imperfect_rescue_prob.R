## ============================================================================
## pipeline/63_imperfect_rescue_prob.R -- the probability of rescue in every
## cell of the phi x eta grid, from the replicates (input of Fig. 5 and
## Extended Data Figs. 6 and 7).
## ----------------------------------------------------------------------------
## A replicate is RESCUED when its tail-mean dead-cell load L = Ea + En is below
## the efferocytosis half-saturation point Kg (the Hill midpoint of g(L)):
##     p_open(unit, phi, eta) = #{rep : L_tail < Kg} / n_rep
## where n_rep is the number of surviving replicates of that cell.
##
## PER UNIT: p_mean = mean of p_open over the 121 cells. Units with p_mean = 0
## (never rescued) or p_mean = 1 (rescued in every replicate of every cell) are
## excluded by the figure scripts; the rest are the CROSSING units.
## Also written (not used by any figure): the composite benefit
## res_off - res_mean per cell, marginal ranges and normalised areas of p and of
## benefit along each axis, and a switch weight.
##
## in : ./data/fig05_raw_all.rds  (per-replicate rows, from 62_imperfect_read.R)
##      ./data/fig05_cells.rds    (only the untreated reference res_off)
## out: ./data/fig05_switch_prob_cells.rds   one row per (unit, phi, eta)
##      ./data/fig05_switch_prob_units.rds   one row per unit
## RUN: Rscript pipeline/63_imperfect_rescue_prob.R     (about a minute)
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages(library(data.table))
RAW = Sys.getenv("RAW", "./data/fig05_raw_all.rds")
if (!file.exists(RAW))
  stop(RAW, " not found -- copy fig05_raw_all.rds from the cluster scratch to ./data first")
KE = 100
B = as.data.table(readRDS(RAW))
setnames(B, c("phi_instability", "treg_eff"), c("phi", "eta"), skip_absent = TRUE)
B = B[theta_treg == 1L & rep > 0L]              # the treated grid; the untreated arm is rep 0
cat(sprintf("[prob] %s treated runs | %d units | %d replicates per cell\n",
            format(nrow(B), big.mark = ","), uniqueN(B$unit), B[, uniqueN(rep)]))
B[, `:=`(L = Ea_mean + En_mean, stab = 1 - phi)]
B[, open := L < Kg]
B[, res := sqrt((Ei_mean / KE)^2 + (PLP_mean / KE)^2)]

## ---- per cell: the probability of rescue, and the benefit as before ---------
OFF = as.data.table(readRDS("./data/fig05_cells.rds"))[, .(unit, res_off = res_off[1]), by = unit][, .(unit, res_off)]
CELL = B[, .(n_rep = .N, n_open = sum(open), p_open = mean(open),
             res_mean = mean(res), L_mean = mean(L)), by = .(unit, phi, eta, stab, Kg)]
CELL = merge(CELL, unique(OFF), by = "unit")
CELL[, benefit := res_off - res_mean]
cat(sprintf("[prob] cells: %d | split cells (0 < p < 1): %d (%.1f%%)\n",
            nrow(CELL), CELL[p_open > 0 & p_open < 1, .N], 100 * mean(CELL$p_open > 0 & CELL$p_open < 1)))
saveRDS(CELL, "./data/fig05_switch_prob_cells.rds")

## ---- per unit: marginals of p and of benefit, on each axis -----------------
area_id = function(x, y) { o = order(x); x = x[o]; y = y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) - 0.5 }
norm_area = function(x, f) { r = max(f) - min(f); if (!is.finite(r) || r <= 0) return(NA_real_)
  area_id(x, (f - min(f)) / r) }
U = CELL[, {
  ps = .SD[, .(v = mean(p_open)),  by = stab][order(stab)]
  pe = .SD[, .(v = mean(p_open)),  by = eta ][order(eta)]
  bs = .SD[, .(v = mean(benefit)), by = stab][order(stab)]
  be = .SD[, .(v = mean(benefit)), by = eta ][order(eta)]
  dlt = if (any(p_open > 0.5) && any(p_open <= 0.5))
          mean(benefit[p_open > 0.5]) - mean(benefit[p_open <= 0.5]) else NA_real_
  .(p_mean = mean(p_open), n_split = sum(p_open > 0 & p_open < 1),
    p_range_stab = max(ps$v) - min(ps$v), p_range_eta = max(pe$v) - min(pe$v),
    parea_stab = norm_area(ps$stab, ps$v), parea_eta = norm_area(pe$eta, pe$v),
    area_stab  = norm_area(bs$stab, bs$v), area_eta  = norm_area(be$eta, be$v),
    ## raw (unnormalised) area of the p marginal against the identity: p is
    ## already on [0,1], so this is meaningful where the normalised one is not
    parea_stab_raw = area_id(ps$stab, ps$v), parea_eta_raw = area_id(pe$eta, pe$v),
    delta = dlt,
    w_stab = if (is.finite(dlt)) dlt * (max(ps$v)-min(ps$v)) / (max(bs$v)-min(bs$v)) else NA_real_,
    w_eta  = if (is.finite(dlt)) dlt * (max(pe$v)-min(pe$v)) / (max(be$v)-min(be$v)) else NA_real_)
}, by = unit]
U[, `:=`(w_stab = pmin(pmax(w_stab, 0), 1), w_eta = pmin(pmax(w_eta, 0), 1))]
U = merge(U, unique(B[, .(unit, Kg, rescue_mode)]), by = "unit")
saveRDS(U, "./data/fig05_switch_prob_units.rds")
cat(sprintf("[prob] units: %d | never rescued %d | always rescued %d | crossing %d\n",
            nrow(U), U[p_mean == 0, .N], U[p_mean == 1, .N], U[p_mean > 0 & p_mean < 1, .N]))
X = U[p_mean > 0 & p_mean < 1]
cat(sprintf("[prob] on crossing units, p-marginal area vs benefit area: stability rho %+.3f | eta rho %+.3f\n",
            cor(X$parea_stab, X$area_stab, method = "spearman", use = "complete.obs"),
            cor(X$parea_eta,  X$area_eta,  method = "spearman", use = "complete.obs")))
cat("[prob] wrote ./data/fig05_switch_prob_cells.rds and ./data/fig05_switch_prob_units.rds\n")
