## ============================================================================
## pipeline/22_corridor_attractors.R -- extract the OFF system's low- and high-burden attractors,
## per coordinate, from the sweep output.
##
## WHY THIS EXISTS. 21_corridor_tzero.R stored span = |hi_off - lo_off| but not lo_off itself, so
## every corridor distance in corridor_tzero_<BRANCH>.rds is measured from the
## BASE STATE. Measuring from lo_off instead needs that one number back.
##
## NOTHING IS RECOMPUTED. pipeline/01_sweep.R writes the full 13-state of every
## landmark as plain columns -- lo_off_<state>, hi_off_<state>, one per state in
## STATE_ORDER -- and they survive into selected_big.rds. pipeline/03_select.R already
## reads them (gk("lo_off_Eh"), gk("lo_off_M1")). This is a reshape, not a
## computation.
##
## NOT skeleton_<BRANCH>.rds: that file is a hand-picked column subset written by
## pipeline/03_select.R and carries only the Eh, M1 and PLP landmarks.
##
## THE SPAN CHECK IS THE POINT. 21_corridor_tzero.R computed span from ITS OWN finder call. If
## that call and the sweep's found the same attractors, |hi_off - lo_off| here
## reproduces the stored span exactly. Where it does not, the two stages
## disagreed about the fixed points and those rows must not be used.
##
## in : ./data/corridor_tzero_<BRANCH><CORR_TAG>.rds, ./data/selected_big.rds
## out: ./data/corridor_lo_off_<BRANCH><CORR_TAG>.rds   (long: one row per combination x coordinate)
##      ./data/corridor_lo_off_<BRANCH><CORR_TAG>_disagree.rds (combinations whose span differs > 1%)
## RUN: BRANCH=sterile CORR_TAG=_round1 Rscript pipeline/22_corridor_attractors.R
##      BRANCH=pathogenic_help Rscript pipeline/22_corridor_attractors.R
##
## selected_big.rds is NOT branch-suffixed -- pipeline/03_select.R overwrites it on every
## run. If it was last written for a different branch it will not cover this
## branch's combinations, and the coverage check below stops rather than joining
## a partial answer.
## ============================================================================
suppressMessages(library(data.table))
SRC = Sys.getenv("SRC", "."); setwd(SRC)
BRANCH = Sys.getenv("BRANCH", "sterile")
## CORR_TAG names a corridor file that is not the final one -- "_round1" for the
## uniform random subsample recovered from a timed-out array. It is appended to
## the branch in every corridor input and output name, so preliminary products
## never overwrite final ones.
CTAG = Sys.getenv("CORR_TAG", "")
CK = c("tau_id","Kg","beta_ros","r_L_P")

SWEEP = c("CLP","Eh","Ea","En","M1","R","D")
if (BRANCH != "sterile") SWEEP = c(SWEEP, "PLP")

CORR = sprintf("./data/corridor_tzero_%s%s.rds", BRANCH, CTAG)
if (!file.exists(CORR)) stop("no ", CORR)
L = as.data.table(readRDS(CORR))[note == "ok"]
if (all(c("lo_off","hi_off") %in% names(L)))
  cat("[corridor_attractors] NOTE: the corridor file already carries lo_off -- this run is redundant\n")
U = unique(L[, ..CK])
cat(sprintf("[corridor_attractors] %s combinations in %s\n", format(nrow(U), big.mark = ","), CORR))

SB = Sys.getenv("SEL_FILE", "./data/selected_big.rds")
if (!file.exists(SB)) stop("no ", SB, " -- run pipeline/03_select.R for BRANCH=", BRANCH, " first")
S = as.data.table(readRDS(SB))
NEED = c(paste0("lo_off_", SWEEP), paste0("hi_off_", SWEEP))
MISS = setdiff(NEED, names(S))
if (length(MISS)) stop(SB, " is missing ", paste(MISS, collapse = ", "))

S = merge(U, S[, c(CK, NEED), with = FALSE], by = CK, all.x = TRUE)
if (nrow(S) != nrow(U)) stop("the join changed the row count: ", nrow(U), " -> ", nrow(S),
                             " -- duplicate keys in ", SB)
n_un = sum(is.na(S[[NEED[1]]]))
if (n_un)
  stop(format(n_un, big.mark = ","), " of ", format(nrow(U), big.mark = ","),
       " corridor combinations are absent from ", SB, ".\n",
       "  That file is not branch-suffixed -- it was probably last written by\n",
       "  pipeline/03_select.R for a different branch. Re-run:  BRANCH=", BRANCH,
       " Rscript pipeline/03_select.R")

## ---- wide to long -----------------------------------------------------------
LO = rbindlist(lapply(SWEEP, function(v)
  data.table(S[, ..CK], coord = v,
             lo_off = as.numeric(S[[paste0("lo_off_", v)]]),
             hi_off = as.numeric(S[[paste0("hi_off_", v)]]))), fill = TRUE)

OUT = sprintf("./data/corridor_lo_off_%s%s.rds", BRANCH, CTAG)
saveRDS(LO, OUT)
cat(sprintf("[corridor_attractors] wrote %s (%s rows, %s combinations)\n", OUT,
            format(nrow(LO), big.mark = ","), format(uniqueN(LO[, ..CK]), big.mark = ",")))

## ---- CONSISTENCY: does this reproduce the span 21_corridor_tzero.R stored? -------------------
CH = merge(unique(L[, c(CK, "coord", "span"), with = FALSE]), LO, by = c(CK, "coord"))
CH[, span_here := abs(hi_off - lo_off)]
CH[, rel := abs(span_here - span) / pmax(abs(span), 1e-12)]
cat("\n[corridor_attractors] span reproduced from the stored attractors (relative difference):\n")
print(CH[, .(n = .N, med = signif(median(rel, na.rm = TRUE), 3),
             p99 = signif(quantile(rel, 0.99, na.rm = TRUE), 3),
             max = signif(max(rel, na.rm = TRUE), 3),
             n_over_1pct = sum(rel > 0.01, na.rm = TRUE)), by = coord][order(-max)])
cat("[corridor_attractors] medians at machine precision mean the sweep and 21_corridor_tzero.R found the same\n")
cat("      attractors. Rows in n_over_1pct did not agree and must be excluded.\n")
BAD = unique(CH[rel > 0.01, ..CK])
if (nrow(BAD)) {
  BF = sprintf("./data/corridor_lo_off_%s%s_disagree.rds", BRANCH, CTAG)
  saveRDS(BAD, BF)
  cat(sprintf("[corridor_attractors] %s combinations disagree on at least one coordinate -> %s\n",
              format(nrow(BAD), big.mark = ","), BF))
}
