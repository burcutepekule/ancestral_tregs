## ============================================================================
## lib/plp_clean.R -- the ONE rule for negative P_LP. Sourced by every reader.
##
## WHY THIS EXISTS
## P_LP = 0 is an invariant manifold of the system: every term in dPLP carries a
## factor of PLP or of PL, so a cleared pathogen approaches zero and never leaves
## it. lsoda does not know that -- it straddles the manifold on approach, and
## 8.8% of converged endpoints come back very slightly negative (median -1.6e-21,
## 99.7% inside 1e-9). That is floating point, and it is floored to 0.
##
## A small tail is NOT floating point. dPLP contains the logistic self-limitation
##      PLP * r_LP_P * (1 - PLP / K_LP_P)
## which at PLP = -x evaluates to  -r*x*(1 + x/K_LP_P) ~ -r*x^2/K_LP_P: negative
## and accelerating. That is finite-time blow-up to -Inf, and nothing in the
## equations pulls it back -- the influx term is the only restoring force and it
## is bounded. Measured worst case: P_LP = -1.9e11 at t = TMAX. Those rows are
## failed integrations wearing the costume of a fixed point.
##
## They must be removed AT INGEST, not in a figure script, because comp =
## sqrt((K_E - Eh)^2 + P_LP^2) inherits the garbage and comp is the ruler that
## every attractor count, rescue gap and basin classification downstream uses. A
## single poisoned endpoint inside a system makes that system look like it has an
## extra distinct attractor.
##
## THE PERMANENT FIX is structural: clamp PLP at 0 inside the RHS so the
## quadratic can never engage. That requires regenerating the sweeps, so it goes
## in with the next run; this file makes the existing data safe meanwhile.
##
## Measured impact of the filter (effero_knockout_off.rds, 3.28M rows):
##   pathogenic bistable fraction 0.2852417 -> 0.2852167, sterile unchanged.
## i.e. it is a numerical hygiene fix, not a result-changing one. Keep it that
## way -- if a future run shows this filter MOVING a result, the run is wrong.
## ============================================================================

PLP_TOL = 1e-9        # below -PLP_TOL = diverged; between -PLP_TOL and 0 = roundoff

## Two modes, because the two data shapes need different remedies:
##
##   mode = "row"  -- LONG tables where one row IS one endpoint
##                    (effero_knockout_off.rds: row_ind x arm x Eh0). A diverged
##                    endpoint is a dead row; drop it.
##   mode = "cell" -- WIDE tables where one row is a whole combo and the P_LP
##                    column is a summary of it (the sweep output). Dropping the
##                    row would throw away the combo, and NA in these columns is
##                    already meaningful ("no attractor of this kind"), so a
##                    diverged value is set to NA instead and the row survives.
##
## In both modes: values in [-tol, 0) are roundoff and are floored to 0, and comp
## is recomputed afterwards if it is present.
##
##   D     : data.table
##   plp   : name of the P_LP column
##   comp  : name of the composite-burden column to recompute (NULL to skip)
##   eh    : name of the epithelium column comp is built from
##   K_E   : epithelial carrying capacity used by comp
##   label : printed tag so the log line says which table/column was cleaned
plp_clean = function(D, plp = "PLP", comp = "comp", eh = "Eh", K_E = 100,
                     tol = PLP_TOL, mode = c("row", "cell"), label = "") {
  mode = match.arg(mode)
  tag  = if (nzchar(label)) paste0(" ", label) else ""
  if (!inherits(D, "data.table")) D = as.data.table(D)
  if (!plp %in% names(D)) { cat(sprintf("[plp_clean%s] no '%s' column, nothing to do\n", tag, plp)); return(D) }
  n0    = nrow(D)
  v     = D[[plp]]
  nblow = length(which(v < -tol))                 # which() drops NA
  if (mode == "row") {
    nnaP = sum(!is.finite(v))
    D = D[which(is.finite(v) & v >= -tol)]
    data.table::set(D, j = plp, value = pmax(D[[plp]], 0))
    msg = sprintf("dropped %d diverged + %d non-finite rows | %s kept",
                  nblow, nnaP, format(nrow(D), big.mark = ","))
  } else {
    w = which(v < -tol)
    if (length(w)) data.table::set(D, i = w, j = plp, value = NA_real_)
    data.table::set(D, j = plp, value = pmax(D[[plp]], 0))   # NA stays NA
    msg = sprintf("blanked %d diverged cells | %s rows kept", nblow, format(nrow(D), big.mark = ","))
  }
  if (!is.null(comp) && comp %in% names(D) && eh %in% names(D))
    data.table::set(D, j = comp, value = sqrt((K_E - D[[eh]])^2 + D[[plp]]^2))
  cat(sprintf("[plp_clean%s] %s rows in | %s\n", tag, format(n0, big.mark = ","), msg))
  D[]
}
