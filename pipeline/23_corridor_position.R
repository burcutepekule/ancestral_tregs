## ============================================================================
## pipeline/23_corridor_position.R -- the corridor as a position, not a width;
## also classifies every (combination, coordinate) row by q_kind, which Supplementary Fig. 1B/C filter on.
##
## THE MEASURE. Distances are taken from the base state -- the injury initial
## condition the sweep actually perturbs -- out to each wall along coordinate X:
##
##      d_res_off = |x_off - base|        the OFF wall
##      d_res_on  = |x_on  - base|        the ON wall
##      q         = d_res_on / (d_res_on + d_res_off)
##
## q = 0.5 is no Treg effect, above is help, below is harm. Equivalently
## q - 1/2 = w / (2 (d_res_on + d_res_off)) with w = d_res_on - d_res_off, so the
## numerator is the corridor width as before and only the denominator changes.
##
## WHY NOT THE BASIN SPAN. 21_corridor_tzero.R divided every distance by span = |hi_off - lo_off|,
## the separation between the OFF system's two attractors. That is a property of
## the endpoints; once an initial condition is on one side of the separatrix the
## outcome is fixed and where the attractor sits no longer bears on the question.
## The span also differs between coordinates for reasons that have nothing to do
## with the boundary -- En accumulates at the high-burden attractor because
## nothing clears it, giving span_En / span_Ea = 10.7 while the two raw wall
## distances differ by 5.6% -- and dividing two axes of the same physical unit by
## different numbers rotates the measured normal. q needs no external ruler: it
## is normalised by the corridor's own extent, which the sweep measured.
##
## WHY base AND NOT THE LOW ATTRACTOR. Measuring from lo_off was tried and is
## wrong here. The scan line is base + t*e_X: one coordinate varied, every other
## held at the base state. lo_off is not on that line -- it differs from base in
## every coordinate at once -- so |x_off_X - lo_off_X| is a difference of two
## numbers, not a distance from the attractor to the wall along anything that was
## measured. It also destroys the range: with lo_off in the denominator the whole
## of q lands within 6e-5 of 0.5.
##
## base is not an arbitrary point either. Eh = e0, En = K_E - e0, everything else
## zero is the injury initial condition, and rescue_gap is measured on that same
## line, so the two share an origin.
##
## THE SPAN CANCELS. d_off and d_on as 21_corridor_tzero.R stored them are both divided by span,
## and the ratio d_on / (d_on + d_off) is unchanged by that. q therefore needs no
## ruler at all -- the normalisation question does not arise for it.
##
## CENSORING. Where no ON wall exists within the searched range, d_res_on is
## unbounded and q -> 1. That is reported as q = 1 with q_kind = "on_absent",
## and q_lb gives the lower bound obtained by putting the ON wall AT the search
## ceiling. The pair brackets the truth. Under the width w these combinations
## were infinite and had to be dropped.
##
## in : ./data/corridor_tzero_<BRANCH><CORR_TAG>.rds
##      ./data/corridor_lo_off_<BRANCH><CORR_TAG>.rds   (optional -- only for the geometry
##                                        diagnostic at the end)
## out: ./data/corridor_q_<BRANCH><CORR_TAG>.rds
## RUN: BRANCH=sterile CORR_TAG=_round1 Rscript pipeline/23_corridor_position.R
##      BRANCH=pathogenic_help Rscript pipeline/23_corridor_position.R
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

IN = sprintf("./data/corridor_tzero_%s%s.rds", BRANCH, CTAG)
if (!file.exists(IN)) stop("no ", IN)
M = as.data.table(readRDS(IN))[note == "ok"]
cat(sprintf("[corridor_position] %s: %s rows, %s combinations\n", IN, format(nrow(M), big.mark = ","),
            format(uniqueN(M[, ..CK]), big.mark = ",")))

## ---- lo_off: from the corridor file if it carries it, else from pipeline/22_corridor_attractors.R ---------
HAVE_LO = TRUE
if (!"lo_off" %in% names(M)) {
  LOF = sprintf("./data/corridor_lo_off_%s%s.rds", BRANCH, CTAG)
  if (!file.exists(LOF)) {
    HAVE_LO = FALSE
    cat("[corridor_position] no ", LOF, " -- q does not need it; the geometry diagnostic is skipped\n")
  }
  LO = if (HAVE_LO) as.data.table(readRDS(LOF)) else NULL
  if (HAVE_LO) {
    n0 = nrow(M)
    M = merge(M, LO[, c(CK, "coord", "lo_off", "hi_off"), with = FALSE],
              by = c(CK, "coord"), all.x = TRUE)
    if (nrow(M) != n0) stop("the lo_off join changed the row count: ", n0, " -> ", nrow(M),
                            " -- duplicate keys in ", LOF)
    cat(sprintf("[corridor_position] lo_off joined from %s | %s of %s rows matched\n", LOF,
                format(sum(!is.na(M$lo_off)), big.mark = ","), format(nrow(M), big.mark = ",")))
  }
} else {
  cat("[corridor_position] lo_off read from the corridor file directly\n")
}

## ---- the three distances ----------------------------------------------------
M[, d_res_off  := abs(x_off   - base)]
M[, d_res_on   := abs(x_on    - base)]
M[, d_res_ceil := abs(hi_used - base)]

## ---- q, and what kind of q it is -------------------------------------------
## The OFF wall is required in every case: without it there is no corridor to
## position. The ON wall is either measured, or provably absent within the search
## range, or ambiguous -- and those three are not the same answer.
M[, q_kind := fifelse(
    !is.finite(d_res_off),                     "no_off_wall",
    fifelse(is.finite(d_res_on),               "measured",
    fifelse(is.na(x_on) & n_t1 == 0L,          "on_absent",
                                               "on_ambiguous")))]
M[, q    := NA_real_]
M[, q_lb := NA_real_]
M[q_kind == "measured",  q    := d_res_on / (d_res_on + d_res_off)]
M[q_kind == "measured",  q_lb := q]
## ON wall not found anywhere in the searched range: d_res_on is unbounded, so q
## is 1. Placing it at the ceiling instead gives the lower bound.
M[q_kind == "on_absent", q    := 1]
M[q_kind == "on_absent", q_lb := d_res_ceil / (d_res_ceil + d_res_off)]

## ---- GEOMETRY CHECK ---------------------------------------------------------
## Not a validity test for q -- q is measured entirely along the scan line and
## does not depend on where the attractor sits. This reports whether the OFF wall
## lies beyond the attractor's coordinate or short of it, which says something
## about the shape of the basin and is worth seeing. Reported, never filtered on.
if ("lo_off" %in% names(M)) {
  M[, wall_past_attractor := (x_off - base) * (x_off - lo_off) > 0]
  cat("\n[corridor_position] is the OFF wall further from base than the attractor is?\n")
  print(M[is.finite(d_res_off) & is.finite(lo_off),
          .(n = .N, past = sum(wall_past_attractor, na.rm = TRUE),
            pct = round(100 * mean(wall_past_attractor, na.rm = TRUE), 1)),
          by = coord][order(-pct)])
}

## ---- report -----------------------------------------------------------------
cat("\n[corridor_position] q by kind:\n")
print(M[, .N, by = q_kind][order(-N)])

cat("\n[corridor_position] q by coordinate (measured rows only):\n")
print(M[q_kind == "measured", .(n = .N,
        med_q = round(median(q, na.rm = TRUE), 4),
        q25   = round(quantile(q, 0.25, na.rm = TRUE), 4),
        q75   = round(quantile(q, 0.75, na.rm = TRUE), 4),
        frac_help = round(mean(q > 0.5, na.rm = TRUE), 3)), by = coord][order(-med_q)])

cat("\n[corridor_position] q against the rescue gap, by coordinate and class:\n")
RH = function(d) if (sum(is.finite(d$q) & is.finite(d$rescue_gap)) > 2)
       suppressWarnings(cor(d$q, d$rescue_gap, method = "spearman",
                            use = "complete.obs")) else NA_real_
print(M[q_kind == "measured" & !gap_censored,
        .(n = .N, rho = round(RH(.SD), 3)), by = .(coord, rescue_mode),
        .SDcols = c("q","rescue_gap")][order(coord, rescue_mode)])

OUT = sprintf("./data/corridor_q_%s%s.rds", BRANCH, CTAG)
saveRDS(M, OUT)
cat(sprintf("\n[corridor_position] wrote %s (%s rows)\n", OUT, format(nrow(M), big.mark = ",")))
