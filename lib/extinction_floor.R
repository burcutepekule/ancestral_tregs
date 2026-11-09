## ============================================================================
## lib/extinction_floor.R -- the pathogen extinction floor and the solver
## settings that make it meaningful. ONE definition, sourced by everything that
## integrates: pipeline/01_sweep.R, lib/skeleton_defs.R, and through the latter every
## figure script. Requires STATE_ORDER to be defined already; it is asserted
## rather than assumed, because a mis-ordered state vector would put the root on
## the wrong compartments silently.
## ----------------------------------------------------------------------------
## WHY IT EXISTS. The state variables are continuous, so nothing in the system
## stops a pathogen population from regrowing out of an arbitrarily small value.
## In the pathogenic branch P_LP is driven below 1e-100 in the combinations where
## the Tregs OFF arm clears the infection, and it can still return to a load of
## ten or more once M1 is suppressed. Without a floor the outcome is decided by
## the sign of the integrator's rounding error: a state of +1e-112 regrows over
## 111 decades and one of -1e-113 is held at zero by the pmax(y, 0) in rhs_app.
## In one combination two Treg-induction delays 35 time units apart separated
## "harm" from "help" on that basis alone, with final E_h of 81.9 against 99.4.
##
## THE VALUE, DERIVED RATHER THAN CHOSEN. One conversion is used everywhere in
## this model: 1 state unit = 1e11 bacteria. K_LP_P = K_LP_C = 100 units is then
## the 1e13 organisms a mucosal compartment can carry, K_L_P = 1000 units is 1e14
## in the lumen, and one bacterium is 1e-11 units in every compartment. A
## compartment holding less than one organism holds none, so that is the floor.
## State the conversion once, in the methods, and both compartments follow from
## it; do not give PL a floor of its own derived from its own capacity.
##
## EXT_TOL AND atol ARE NOT INDEPENDENT. lsoda accepts a local error of
## rtol*|y| + atol, so at deSolve's default atol = 1e-6 the solver is permitted an
## absolute error a hundred thousand times larger than the threshold this root
## tests: the event would fire on noise, or fail to fire because noise held the
## value above the floor. atol must sit well below EXT_TOL, and the assertion at
## the bottom of this file enforces it. Raising EXT_TOL instead is not an option
## that survives contact with the units -- 1e-4 would be declaring extinction at
## ten million organisms.
## ============================================================================
if (!exists("STATE_ORDER")) stop("extinction_floor.R: STATE_ORDER must be defined first")

EXT_TOL = as.numeric(Sys.getenv("EXT_TOL", "1e-11"))

## positional, not by name: rhs_app reads its state through
## setNames(as.list(y), STATE_ORDER), so it is the order that carries meaning.
EXT_IDX = match(c("PL", "PLP"), STATE_ORDER)
stopifnot(!anyNA(EXT_IDX))

## HOW IT IS APPLIED. Not by clamping inside the derivative: a hard cut there is a
## discontinuity lsoda chases, and it fails with an excessive amount of work at
## the crossing. It is a solver ROOT instead. The root is on the LARGER of the two
## compartments, so it is crossed once, when BOTH are below the floor; the event
## then sets both to exactly zero. Every term of dPL carries a factor of PL and
## every term of dPLP carries a factor of PL or PLP, so from that point both
## derivatives are zero, the root function rests at a constant -EXT_TOL, and the
## pathogen stays extinct without any modification to rhs_app. Rooting on either
## compartment alone would chatter: setting P_LP to zero while PL is still above
## the floor leaves the influx term PL*gamma_B_P*(1-Eh/K_E)*(1-PLP/K_LP_P)
## positive, P_LP rises back through the floor, and the event fires again.
##
## THE COMMENSAL COMPARTMENT IS DELIBERATELY NOT FLOORED. A state needs a floor
## exactly when every term of its derivative carries a factor of itself, and only
## PL and PLP satisfy that. CLP does not: its influx carries CL, a fixed reservoir
## at 1000, so it has no autonomous growth and no extinction to represent, and
## flooring it would freeze a compartment the model intends to resupply. Nor does
## anything else -- Eh regenerates through kappa_E*(K_E - Eh), Ea's source is
## Eh*tau_r*R, En's is Eh*tau_P_B*(PLP/K_LP_P) + mu_Ea*Ea, R's is beta_ros*M1,
## D's is sigma_D*(...), S's is sigma_S*Tr, and M0+M1 and T0+Tr are conserved
## exactly since dM0 = -dM1 and dT0 = -dTr.
##
## THE FINDER IS UNAFFECTED. PL = PLP = 0 zeroes every term of both derivatives,
## so it is a fixed point of the UNFLOORED algebra too, and FP_finder_v6.R already
## returns it: the sterile solutions reappear in the pathogenic finder as the
## PLP = 0 branch, and endpoint seeding tries plp_try = 1e-12 for exactly the case
## where lsoda has undershot P_LP below zero at a pathogen-free attractor. Floored
## endpoints therefore have fixed points to match against, and snap_fp does not
## need changing.
##
## OUTPUT ROWS. lsodar appends a row at the root time, so a run can return more
## rows than the requested output grid. Every consumer in this repo reads the LAST
## row (.N, or nrow(out)) rather than row 2 -- check any new one that does not.
rootfun_ext = function(t, y, p) max(y[EXT_IDX]) - EXT_TOL
event_ext   = function(t, y, p) { y[EXT_IDX] = 0; y }
EVENTS_EXT  = list(func = event_ext, root = TRUE)

## ---- SOLVER SETTINGS, ONE DEFINITION FOR EVERY CALLER ----------------------
## atol is what the small states need and what the root above depends on. rtol is
## looser because it binds on the large states (E_h ~ 100), where it costs the
## most and buys the least. hmax IS NOT OPTIONAL -- see the long note at the ode()
## call in pipeline/01_sweep.R: deSolve defaults it to max(diff(times)), which on a
## two-point output grid is the whole horizon, and lsoda then strides the stiff
## tail and returns silently wrong states in every compartment.
if (!exists("ATOL"))     ATOL     = as.numeric(Sys.getenv("ATOL", "1e-14"))
if (!exists("RTOL"))     RTOL     = as.numeric(Sys.getenv("RTOL", "1e-8"))
if (!exists("HMAX"))     HMAX     = as.numeric(Sys.getenv("HMAX", "5"))
if (!exists("MAXSTEPS")) MAXSTEPS = as.integer(Sys.getenv("MAXSTEPS", "200000"))

if (!(ATOL < EXT_TOL / 100))
  stop(sprintf("extinction_floor.R: atol = %g is not at least two decades below EXT_TOL = %g; the root would fire on integrator noise",
               ATOL, EXT_TOL))

## ---- THE CALL EVERY ANALYSIS SCRIPT SHOULD USE ------------------------------
## pipeline/01_sweep.R writes its ode() call out in full, because its long comment block
## is the record of why each setting is what it is. Everything downstream should
## use this instead, so that a figure cannot quietly integrate different dynamics
## from the sweep that produced its labels. That is not hypothetical: after the
## floor was folded in, six call sites were still on plain lsoda at deSolve's
## 1e-6 defaults with no root, including the efferocytosis knockout -- whose whole
## point is to remove clearance and see whether the pathogen returns.
##
## IT RESTRICTS THE OUTPUT TO THE REQUESTED TIMES. lsodar appends a row at the
## root time, so a run on seq(0, TMAX, by = 2) can come back one row longer than
## the grid. Endpoint readers do not care, but anything pairing the result with
## its own time vector would silently misalign. Trimming here means every caller
## behaves exactly as it did before, with only the dynamics corrected. A run that
## stopped early still returns short, so existing nrow() and final-time guards
## keep working.
ode_ext = function(y, times, func, parms, hmax = HMAX, maxsteps = MAXSTEPS, ...) {
  o = deSolve::ode(y, times, func, parms, method = "lsodar",
                   rootfun = rootfun_ext, events = EVENTS_EXT,
                   atol = ATOL, rtol = RTOL, hmax = hmax, maxsteps = maxsteps, ...)
  if (is.null(o) || !nrow(o)) return(o)
  keep = match(times, o[, 1]); keep = keep[!is.na(keep)]
  if (length(keep)) o[keep, , drop = FALSE] else o
}

cat(sprintf("[extinction_floor] root on %s at EXT_TOL = %g | atol %g rtol %g hmax %g maxsteps %d\n",
            paste(STATE_ORDER[EXT_IDX], collapse = ", "), EXT_TOL,
            ATOL, RTOL, HMAX, MAXSTEPS))
