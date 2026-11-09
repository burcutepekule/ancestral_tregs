SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## 01_sweep.R -- one row per (row_ind, Kg, beta_ros, r_L_P), carrying the
## FULL fixed-point set as list-columns instead of two hand-picked landmarks.
##
## WHAT v6 CHANGES, AND WHY
## v5 probed Eh0 at spacing 5 and reported a threshold quantised to that spacing;
## a separate later stage re-derived the fine threshold by bisection, using its own
## call to the finder. That stage is gone. Everything it did happens here, on the
## same integrations, and the two stages can no longer disagree.
##
## THE FOUR STEPS, in order. Each one exists because the alternative was measured
## to fail.
##
## 1. PROBE.  Integrate Eh0 = 0, 10, ... 90, 99 and keep the final states.
##    Not 100: in the sterile system Eh0 = 100 gives Eh = K_E, Ea = En = 0, which
##    IS a fixed point -- measured max|x(T) - x(0)| = 0 over TMAX = 12000. It is a
##    saddle (one positive eigenvalue, +0.059, along CLP), so it is never the
##    low-burden attractor, but a final state resting on it is not in any basin and
##    must not be labelled as though it were. 99 escapes it long before TMAX.
##    Spacing 10 rather than 5: compared directly on 20 combinations / 36 arms
##    across both branches, the two spacings returned the same crossing count every
##    time, and 11 probe points plus ~7 bisection halvings costs less than the 19
##    probe points v5 spent.
##
## 2. SEED THE FINDER FROM THOSE FINAL STATES ONLY.  A final state sits ON the
##    attractor it reached, so Newton started there converges to it exactly. The
##    grid pass runs separately and answers a different question (what EXISTS,
##    saddles included); it must never decide whether a threshold is looked for.
##    Called with the grid alone, the pathogenic finder missed the reached
##    attractor on every one of six combinations tested -- returning fixed points
##    only at PLP = 0 while the line settled at PLP = 13.5 -- and each was recorded
##    as monostable. Endpoint seeding also makes "monostable" mean exactly "no
##    second attractor was reached from this line", which is the statement the
##    analysis needs, and makes it equivalent to "no crossing".
##
## 3. THREE LABELS, NOT TWO.  snap_fp returns the nearest fixed point of ANY type,
##    so the match must then be checked for stability:
##        low basin        matched a stable fixed point, in the lowest-burden cluster
##        high basin       matched a stable fixed point, in any other cluster
##        unassigned ('.')  no match, or the match is not stable, or no integration
##    Without the third, a final state resting on an unstable fixed point reads as
##    high basin, and dev cannot catch it -- the state is exactly at the point it
##    matched, so the deviation is zero. Clustering with RAY_CLUST_TOL comes first,
##    so two attractors closer than the tolerance cannot manufacture a crossing.
##
## 4. BISECT ONLY WHERE A CROSSING WAS SEEN.
##        exactly one crossing -> bisect that bracket to EH0_TOL. Eh0_star is real.
##        more than one        -> n_cross recorded, Eh0_star = NA. There is no
##                                single threshold; returning one is what testing
##                                only the two bracket ends used to do silently.
##        none                 -> Eh0_star = NA, note says which of the two cases.
##    The direction is read off the labels, never assumed: an arm can cross going
##    up in Eh0 or going down, and both were observed.
##
## NAMING. resolved/chronic is retired in favour of low basin / high basin, which
## is what the composite ruler actually orders. hi_lo_gap_off keeps its name.
##
## WHAT IS UNCHANGED FROM v5, AND MUST STAY SO
##   * fp_<cond>_<state> are LIST-COLUMNS -- one numeric vector per combo, length
##     = however many fixed points the finder returned. Nothing is dropped, and
##     clustering stays a downstream choice you can change without re-running.
##   * map_<cond> indexes into those vectors, so every probed Eh0 resolves to real
##     coordinates.
##   * PLP0 is fixed at PLP0_MAIN.
##   * chunks are written as .rds, not .csv -- list-columns cannot survive CSV.
##
## THE DEVIATION CHECK
## Endpoints are matched to fixed points by nearest-neighbour (snap_fp). That is
## only meaningful if the endpoint actually IS near the fixed point. dev is the
## largest RELATIVE deviation over all probed Eh0 and all 13 states,
##      max |integrated - matched| / (|matched| + 1)
## relative because the states span three orders of magnitude (PL ~ 1000, R ~ 1)
## and an absolute tolerance would be dominated by PL alone. Stored as a raw
## number so the tolerance is applied downstream, never baked in. NOTE its blind
## spot: a state resting exactly on an unstable fixed point has dev = 0, which is
## why step 3 checks stability and does not rely on dev.
##
## RUN: Rscript pipeline/01_sweep.R <n_chunks> <this_chunk> <batch>
##      (normally through pipeline/cluster/sweep_array.sh or sweep_local.sh)
##
## TWO MODES, ONE CODE PATH.
##   full grid (default)  every tau row of ./data/tau_grid_v6_<batch>.rds x the
##                        full 7 Kg x 9 beta_ros x 7 r_L_P block (441 per row).
##   cell list            CELLS=<rds with batch, row_ind, Kg, r_L_P> restricts the
##                        sweep to the tau rows of this batch that appear in it and,
##                        per row, to that row's (Kg, r_L_P) cells, each at every
##                        level in BETA_ROS_VALS (a re-sweep of chosen cells at other
##                        levels; not used by any figure of the paper).
## Everything else -- probe, bisection, finder, tolerances, output schema -- is the
## same in both modes, so a cell-list row is the same measurement as a grid row.
## ============================================================================

## ============================================================================
## THE TWO FINDER PASSES  (inherited design, still current)
##   * ENDPOINT-seeded (grid_seed = FALSE): the ray-REACHABLE attractors. Each
##     endpoint converges back to its own attractor. The IC -> FP map is built
##     from this, and it cannot see anything off-ray.
##   * GRID-seeded (grid_seed = TRUE): what EXISTS in the field, SADDLES
##     INCLUDED. Endpoint seeding can never return a saddle -- trajectories do
##     not end on them. Computed ONCE per (combo, condition) by the caller:
##     fixed points do not depend on initial conditions.
##   * FALLBACK: if endpoint seeding returns nothing, the ray snaps against the
##     grid set instead of discarding the combo. seed_mode_off/on records which
##     path each row took: endpoint / grid_fallback / none.
##   Reported as separate columns (n_ray_* vs n_field_*) so ray-reachability and
##   existence are never conflated.
##
## Output is ONE ROW PER COMBO. The IC -> FP map is a string with one character
## per probed Eh0 (attractor id, '.' = no stable match).
##
## COST. The claim that used to sit here -- that ode() with c(0, TMAX) is cheap
## because lsoda strides the converged tail -- was TRUE about the cost and WRONG
## about the correctness, and the two are the same fact. deSolve sets hmax to
## max(diff(times)) when you do not, so a two-point output grid removes the step
## cap entirely. That is what made it fast, and it is also what produced
## PLP = -4.7e155, Eh = 101.3 (above K_E), and silent basin flips. hmax = 5 is
## now set explicitly; the striding is gone and so is the saving. Budget ~5x the
## per-integration cost of the old runs and RE-MEASURE with TIMEIT=1 before
## submitting -- do not carry the old numbers forward. As a starting point,
## measured on v5: 19 probed Eh0 values at TMAX = 12000 with hmax = 5 is ~6.4 s
## per combo per Treg condition, so ~13 s per combo. v6 integrates 11 probe values
## plus ~7 bisection halvings per arm, i.e. 18 against 19, so the per-combo cost
## should be close to unchanged -- but the combo count per tau row is now
## 5 Kg x 9 beta_ros x 7 r_L_P = 315, against 280. RE-MEASURE with TIMEIT=1 before
## submitting; do not carry any of these numbers forward.
##
## STALE CLAIMS REMOVED FROM THIS HEADER, so nobody acts on them:
##   * "PLP0 as a 4-level factor", "168 integrations per pathogenic combo" --
##     one PLP0 level since v5.
##   * "REQUIRES the negative-PLP fix in lib/FP_finder_v6.R ... lsoda
##     undershoots PLP below zero at a pathogen-cleared attractor" -- that
##     diagnosis was wrong. PLP does not undershoot when the step cap is set:
##     measured min PLP over a whole trajectory is +5.5e-21. The undershoot was
##     the missing hmax, not the model and not the finder.
##   * "Usage: Rscript FP_sweep_rescue_cluster_big.R ... TMAX=4000" -- that
##     script is retired (_to_delete/) and TMAX = 4000 leaves 13.8% of the
##     model's eigen-modes short of the dev gate. See the TMAX note below.
## ============================================================================

suppressMessages({library(data.table); library(deSolve)})
source("./lib/ODE_SYSTEM_v6.R")
source("./lib/FP_finder_v6.R")

## ----------------------------------------------------------------------------
## The finder is lib/FP_finder_v6.R's own dispatcher, find_fixed_points(): it
## routes on r_L_P and forwards extra_seeds to both branches. The shadow
## dispatcher that used to live here existed only because the finder file was
## off-limits at the time; the sterile correction is now in the finder itself.
## ----------------------------------------------------------------------------

## ---- sweep axes -------------------------------------------------------------
## CORRECTED beta_ros: the v4 denormalisation scaled M1 up by K_M = 100 but
## never rescaled beta_ros, which multiplies M1 in dR = beta_ros*M1 - mu_R*R.
## R therefore ran ~100x its intended magnitude (measured: median 20, max 959
## against a nominal K_R = 1), and every term linear in R inherited it --
## tau_r*R (epithelial damage) and alpha_1_r*R (pathogen clearance) alike.
## Measured consequence: Eh* = K_E*kappa_E/(kappa_E + tau_r*R + ...) collapsed
## to a median of 1.97, i.e. 2% of capacity, with 95% of attractors below 26.
## These values are the old grid / 100, which restores R to O(1) and puts both
## balances at their transition instead of ~100x on one side.
## 0.0005 (was 0.0002): below this ROS is simply absent -- DAMAGE = tau_r*R/kappa_E
## ~ 0.06 and CLEAR = alpha_1_r*R/r_L_P ~ 0.03 -- so those levels spent grid on a
## limit that is already understood. The TOP stays at 0.0050 deliberately: R is
## linear in beta_ros and lo_off_R already reaches ~9.7 there against K_R = 1.
## Raising it to 0.05 would put R near 100, i.e. back in the regime the scaling
## bug created (R median 20, max 959), where 30_ros_balance.R measures DAMAGE > 1
## and CLEAR > 1 in ~99% of sets -- one saturated corner, not a tradeoff.
## Spacing stays LINEAR by request.
## beta_ros grid. FLOOR RAISED 0.0005 -> 0.0015.
## The sterile algebra collapses to one dimensionless group,
##      Lambda = tau_r * beta_ros * K_M / (mu_R * kappa_E),
## and the lowest epithelium a parameter set can EVER reach, at full inflammation
## M1 = K_M, is  Eh_min = K_E / (1 + Lambda).  Rows with a high Eh_min cannot
## produce damage no matter what the immune system does -- they are not testing
## anything, they are pinned.
##
## At beta_ros = 0.0005, 21.2% of tau rows cannot push Eh below 80 and 70.2%
## cannot push it below 50. That corner is where the "fully inflamed, intact
## barrier" fixed points came from (e.g. row_ind 1188: M1 = 91.6 with Eh = 89.1,
## a genuine stable equilibrium -- verified by integrating from it, drift 5e-8 --
## but only because Eh_min for that row is 89).
##
##   grid                     dead (Eh_min>80)   saturated (Eh_min<5)   useful
##   0.0005-0.0050 x10 (old)       3.03%               4.78%            92.2%
##   0.0015-0.0050 x8  (this)      0.36%               5.97%            93.7%
##   0.0020-0.0050                 0.14%               6.76%            93.1%
##   0.0015-0.0060                 0.29%               7.87%            91.8%
## Raising the floor further just trades dead rows for saturated ones; 0.0015
## minimises the total. EIGHT levels, not ten: 0.0035/7 = 5e-4 exactly, so the
## grid is 0.0015 0.0020 0.0025 0.0030 0.0035 0.0040 0.0045 0.0050 -- round
## values at the same 5e-4 spacing the old grid had, just with the dead bottom
## two removed. Lambda p05 0.61 | med 3.49 | p95 20.6 (unchanged by the level
## count). Combos per tau row: 5 Kg x 8 beta x 7 r_L_P = 280, was 350, so the
## sweep is 20% cheaper as well.
##
## NOTE the axis caveat this does NOT fix: beta_ros is not what the model feels.
## Lambda is, and tau_r / mu_R / kappa_E vary by orders of magnitude across tau
## rows, so the same beta_ros is different physics in different rows. Raising the
## floor removes the pinned corner; it does not make beta_ros an orthogonal axis.
## NINE levels at 0.001 spacing, 0.0015 .. 0.0095. The floor is unchanged and is
## still what minimises dead rows (see the Lambda note above); the ceiling is
## raised from 0.0050 to open the high-ROS end, which the old grid did not reach.
## BETA_ROS_VALS may be overridden (comma-separated) for a cell-list re-sweep; unset,
## it is the main grid below, built with seq() exactly as the stored sweep was.
.bv = Sys.getenv("BETA_ROS_VALS", "")
BETA_ROS_VALS = if (nzchar(.bv)) as.numeric(strsplit(.bv, ",")[[1]]) else seq(0.0015, 0.0095, by = 0.001)
## Kg: SEVEN levels at uniform spacing 5, 15 .. 45. The previous grid was
## c(20,30,40,50,60), and the change is evidence-led. Kg enters the efferocytosis
## gain as a_floor + (a0-a_floor)*Kg^n_g/(Kg^n_g + L^n_g) with n_g = 4 and
## L = Ea + En, so it is a threshold on the dead-cell load with half-max at
## L = Kg, and a level far outside the range of L the model actually visits sits
## saturated at a0 or floored at a_floor and does no work. That is what the old
## grid showed: of 67,486 gated rows, Kg = 20 held 22,929, 30 held 25,557,
## 40 held 15,862, 50 held 3,057 and 60 held just 81 (0.12%), with no Kg = 60
## rows at all in the sterile branch. The population is concentrated between 20
## and 40, so the resolution goes there. Uniform spacing is deliberate: a
## hand-picked set invites the question of why each level was chosen.
##
## Kg = 15 is the one level with no prior evidence either way -- a pilot at
## Kg = 10 returned a median of 3 fixed points against 7 at Kg = 30-60, i.e.
## below the onset of the multiple-attractor structure this work is about, so
## confirm 15 is above that onset before spending a seventh of the run on it.
## Kg = 50 and 60 are dropped; the claim that the phenomenon thins out at high Kg
## rests on the previous sweep, which covered them, and Kg = 45 shows the decline
## has begun. SEVEN levels is 441 combinations per tau row against the previous
## 315, i.e. +40% on the wall clock -- the array sizing in cluster/sweep_array.sh assumes it.
KG_VALS       = c(15, 20, 25, 30, 35, 40, 45)   # 7
RLP_VALS      = seq(0.00, 0.48, 0.08)      # 7  (0 => sterile)
## THE Eh0 LINE. Initial conditions Eh = Eh0, En = K_E - Eh0, Ea = 0, everything
## else at rest (plus PL and PLP when pathogenic). Two ends are excluded and the
## reason is measured, not stylistic:
##   Eh0 = 100 gives Eh = K_E, Ea = En = 0, which IS an equilibrium of the sterile
##     system -- max|x(T) - x(0)| = 0 over TMAX = 12000. Its composite is 0, below
##     every attractor. It is a saddle (leading eigenvalue +0.059, eigenvector
##     dominated by CLP: drop Eh a little, commensals translocate, M1 activates,
##     Eh falls further), so it is never the low-burden attractor -- but a final
##     state resting on it is in no basin at all and must be labelled unassigned
##     rather than high basin. 99 escapes it long before TMAX and is used instead.
##   Eh0 = 0 is kept: it moves, and it is the most-damaged end of the line.
## Spacing 10, not 5. Directly compared on 20 combinations / 36 arms across both
## branches, spacing 5 and spacing 10 gave the same crossing count in every case,
## and the threshold is now bisected rather than read off the spacing, so the
## probe only has to FIND the bracket, not resolve it.
## CHECK THIS IN THE FIRST OUTPUT: rescue_gap should no longer pile up at exactly
## 0 the way it did when the threshold was quantised to the probe spacing. If it
## still does, the gaps really are below EH0_TOL and that is a finding, not an
## artefact -- but confirm it against Eh0_off_coarse before believing it.
## THE PROBE. See step 1 in the header for why the spacing is 10 and why the top
## is 99 and not 100. EH0_SCAN is the set of Eh0 values integrated before anything
## is decided; every label, every crossing count and every bisection bracket comes
## from it, and PROBE_BY is written into the output so that any "no crossing" claim
## carries the resolution it was established at.
EH0_LO   = as.numeric(Sys.getenv("EH0_LO", "0"))
EH0_HI   = as.numeric(Sys.getenv("EH0_HI", "99"))
PROBE_BY = as.numeric(Sys.getenv("PROBE_BY", "10"))
EH0_SCAN = unique(c(seq(EH0_LO, EH0_HI, PROBE_BY), EH0_HI))   # 0 10 ... 90 99
## Bisection tolerance for Eh0_star. Cost is logarithmic in it: a PROBE_BY-wide
## bracket needs ceiling(log2(PROBE_BY / EH0_TOL)) halvings, so 7 at these values.
## 0.1 is well inside what the integrator resolves (dev on a correct run is ~1e-5),
## so it is measured rather than invented.
EH0_TOL  = as.numeric(Sys.getenv("EH0_TOL", "0.1"))
## PLP0 is a second axis of insult severity, not a fixed constant. 13_plp0_
## sensitivity.R showed rescue_gap grows with PLP0 (up to 21x) and that ray
## bistability dies out above it (4/10 combos lost it at 50, 9/10 at 90), so
## pinning one value would both understate rescue and hide that dependence.
## 10 is the MAIN level -- it matches every existing result, and it is a small
## inoculum that the damage-gated translocation term amplifies on its own
## rather than a pre-loaded established infection. The others are reported as
## sensitivity. One output row per (combo x PLP0).
PLP0_MAIN     = as.numeric(Sys.getenv("PLP0_MAIN", "10"))   # single level; see header
DEV_STATES    = NULL          # NULL = all 13; set a subset to narrow the dev check
M_TOT = FP_FIXED$K_M; T_TOT = FP_FIXED$K_T; K_E = FP_FIXED$K_E
OUT_PREFIX = Sys.getenv("OUT_PREFIX", "./data/sweep_raw/fp_batch_")

## TMAX = 12000, and here is the argument, because 20000 was a guess and 4000
## was demonstrably too short.
##
## The relevant timescale is the slowest eigenvalue of the fixed point the ray is
## approaching: a mode decays as exp(-|Re(lambda)| * t), so reaching a relative
## residual of eps takes t = ln(1/eps)/|Re(lambda)|. The gate downstream is
## dev_med <= 1e-2, so t = 4.6/|Re(lambda)| is the time to ARRIVE AT the gate --
## not to be safely inside it.
##
## Over all 3,403 finite eigenvalues stored in batch 1 (these come from the
## finder's Jacobian, not from integration, so the hmax bug never touched them):
##
##   |Re(lambda)| quantile   0.1%      1%       5%      25%      50%
##                         5.29e-4  1.07e-3  1.09e-3  1.18e-3  1.55e-3
##   t for residual 1e-2      8703     4307     4210     3907     2967
##   t for residual 1e-4     17407     8613     8421     7814     5934
##
##   fraction of modes still above 1e-2 residual at TMAX:
##     4000  13.811%   |   5000  0.852%   |   8000  0.441%
##    12000   0.000%   |  20000  0.000%
##
## The 13.8% at 4000 is the independent check: it matches the 13.0%/13.1% of
## combos that actually came back with dev_max > 1e-2 in CHECK 3. 12000 is the
## smallest horizon at which no mode in the sample is still above the gate; the
## slowest one (|Re| = 5.29e-4) sits at a 1.7e-3 residual there, ~6x margin.
## 5000 would leave ~0.9% of modes failing and 8000 ~0.44%, i.e. both re-create
## the exclusion problem we just removed. Direct measurement agrees: worst
## deviation from the TMAX=20000 answer over 12 combos x 19 probed Eh0 values was 5.5e-3 at
## 4000, 3.6e-5 at 8000 and 2.4e-7 at 12000.
##
## This matters more than it used to. With hmax capped the cost is roughly LINEAR
## in TMAX (measured 32/55/76/121 s for 4000/8000/12000/20000), where the old
## unbounded-step runs made the horizon nearly free -- and that "free" horizon
## was the bug.
TMAX    = as.numeric(Sys.getenv("TMAX", "12000"))
MAXSTEPS      = as.integer(Sys.getenv("MAXSTEPS", "200000"))
## Maximum internal step. See the long note at the ode() call in ray_one(): with
## the deSolve default (max(diff(times)) = TMAX here) lsoda strides the stiff
## tail and returns wrong answers, sometimes silently. 5 is the largest value
## that was correct on every test seed; 10 already fails.
HMAX          = as.numeric(Sys.getenv("HMAX", "5"))
## ERROR TOLERANCES, SET EXPLICITLY. These used to be deSolve's defaults,
## atol = rtol = 1e-6, and that was wrong for two reasons. First, the basin
## decision rests on P_LP endpoints and the composite squares P_LP, so any state
## below 1e-6 was being classified out of numerical noise -- P_LP was returning
## small NEGATIVE values at its own minimum. Second, and now decisive, the
## extinction root in lib/skeleton_defs.R tests max(PL, PLP) against
## EXT_TOL = 1e-11: lsoda accepts a local error of rtol*|y| + atol, so at
## atol = 1e-6 the tolerance is a hundred thousand times coarser than the
## threshold being tested and the event would fire on noise or never fire.
## atol must sit well below EXT_TOL. rtol is left looser because it binds on the
## large states (E_h ~ 100), where it costs the most and buys the least.
## MEASURE THE COST before submitting: run BATCH=1 TIMEIT=1 TIME_N=5 and compare
## against the old defaults, because the array sizing depends on it.
ATOL          = as.numeric(Sys.getenv("ATOL", "1e-14"))
RTOL          = as.numeric(Sys.getenv("RTOL", "1e-8"))
GRID_FP = as.integer(Sys.getenv("GRID_FP", "1"))   # also run the grid pass (saddles)
ROW_N   = as.integer(Sys.getenv("ROW_N", "0"))
ROW_SEED= as.integer(Sys.getenv("ROW_SEED", "20260808"))
FLUSH_EVERY = max(1L, as.integer(Sys.getenv("FLUSH_EVERY", "1")))

## ---- WHAT MAKES TWO FIXED POINTS ONE ATTRACTOR ------------------------------
## This was a separation on the composite burden sqrt((K_E-Eh)^2 + PLP^2), and
## that was a defect, because the composite is a PROJECTION: two states that
## trade epithelial damage against pathogen load can share a value. The measured
## case is Eh 38.7 / P_LP 0 against Eh 43.7 / P_LP 26.9 -- composites 61.30 and
## 62.40, 1.10 apart against a tolerance of 2, so they merged into one attractor,
## while being 27 apart in P_LP and 5 apart in E_h. RECLASS_mono_pathogenic.R was
## written to measure the damage and found that only 9 of 29 helped and 23 of 31
## harmed ON_mono combinations survived as genuinely monostable once endpoints
## were clustered on the state instead.
##
## THE RULE NOW: two fixed points are the same attractor when they agree on BOTH
## coordinates, |dEh| <= tol AND |dPLP| <= tol -- a Chebyshev ball, not a
## Euclidean one, and one rule with no special cases. A single tolerance for both
## coordinates is not a hidden weighting here, because K_E = K_LP_P = 100 puts
## them on the same scale.
##
## THE TOLERANCE IS 1, NOT THE 2 THE COMPOSITE USED. Chebyshev is the stricter
## test at equal tolerance -- it requires agreement on each coordinate separately,
## where the composite only required agreement on their combination -- so keeping
## 2 would have widened the ball in the coordinate that happened to dominate. 1 is
## also 1e11 organisms in P_LP under the model's conversion, which is the scale at
## which two pathogen loads stop being the same outcome.
##
## WHAT THIS MEANS AT P_LP NEAR ZERO, stated because the methods should say it: a
## pathogen-free fixed point and one carrying P_LP <= 1 merge. There is no
## separate extinction clause -- the tolerance alone decides, in P_LP as in E_h.
## 02_assemble.R counts how often a merge joins a cleared fixed point to a
## carrying one and prints it, so whether this ever happens is a measurement off
## the first batch rather than an assumption either way.
##
## IDENTITY ONLY. This decides how many attractors there are -- bistable against
## monostable. It does NOT decide which outcome is better: that is an ORDERING
## question, it needs a scalar by construction, and it stays on the unweighted
## composite. Help and harm are computed downstream from the recorded endpoint
## components, not baked in here.
RAY_CLUST_TOL = as.numeric(Sys.getenv("RAY_CLUST_TOL", "1.0"))

## Same attractor?  Chebyshev on (Eh, PLP). No special cases.
## The clustering rule now lives in lib/attractor_cluster.R, because
## pipeline/21_corridor_tzero.R needs the same one and had a different rule of its own --
## a gap in the sorted composite, which merges a cleared attractor with one
## carrying 27 units of pathogen. Sourced rather than duplicated so the two
## stages cannot drift apart again. RAY_CLUST_TOL is set above and the functions
## read it as a default argument, so this file keeps control of the tolerance.
source("./lib/attractor_cluster.R")

STATE_ORDER   = FP_STATE_ORDER

## ---- THE EXTINCTION FLOOR ---------------------------------------------------
## Sourced rather than defined here, and sourced by lib/skeleton_defs.R too, so
## the figures cannot integrate under different dynamics or different tolerances
## than the sweep that produced their labels. Provides EXT_TOL, rootfun_ext and
## EVENTS_EXT for the ode() call in ray_one(), and asserts that ATOL is far
## enough below EXT_TOL for the root to mean anything. Needs STATE_ORDER, set
## immediately above.
source("./lib/extinction_floor.R")

FP_MATCH_VARS = c("Eh","Ea","En","PLP","CLP","M1","R","D","Tr","S","PL")
STRUCTURAL = list(K_E=100,K_M=100,K_T=100,K_R=1,K_L_P=1000,K_LP_P=100,K_LP_C=100,CL=1000)
EFFERO     = list(a0=0.0045,a_floor=0.0002,n_g=4)
KIN = c("kappa_E","tau_r","tau_P_B","mu_Ea","mu_R","tau_01","tau_10","tau_1_d","tau_eff",
        "tau_t_a","tau_t_d","alpha_0","alpha_1_e","alpha_1_r","alpha_AMP","gamma_B_P",
        "gamma_B_C","sigma_D","sigma_S","tau_D","tau_S","w_E_over_M")

## NON-NEGATIVITY CLAMP. Same rule as lib/ODE_SYSTEM_v6.R, applied here because
## the sweep integrates THIS closure, not that file's. PLP = 0 is an invariant
## manifold; lsoda straddles it by roundoff; at PLP < 0 the logistic term is
## -r*x^2/K_LP_P, i.e. finite-time escape to -Inf, which is where the 1e154
## entries in dev_max came from. Evaluating the field at the clamped state is a
## no-op for any trajectory that stays in the positive orthant.
rhs_app = function(t,y,p) { y = pmax(y, 0); with(as.list(c(setNames(as.list(y),STATE_ORDER),p)),{
  f_P=PLP/K_LP_P;f_C=CLP/K_LP_C;f_PC=(PLP+CLP)/(K_LP_P+K_LP_C);C_rat=ifelse((f_P+f_C)>0,f_C/(f_P+f_C),0)
  K_Ea=(K_E*tau_r*K_R)/(mu_Ea+a0*K_M);K_En=(K_E*tau_P_B+mu_Ea*K_Ea)/(a0*K_M)
  L=max(Ea+En,0);g=a_floor+(a0-a_floor)*Kg^n_g/(Kg^n_g+L^n_g);effA=g*(M0+M1)*Ea;effN=g*(M0+M1)*En;fb=tau_eff*g*Ea*M1
  dEh=kappa_E*(K_E-Eh)-Eh*(tau_r*R+tau_P_B*(PLP/K_LP_P));dEa=Eh*tau_r*R-mu_Ea*Ea-effA;dEn=Eh*tau_P_B*(PLP/K_LP_P)+mu_Ea*Ea-effN
  dM0=-tau_01*D*M0+tau_10*S*M1+tau_1_d*M1+fb;dM1=tau_01*D*M0-tau_10*S*M1-tau_1_d*M1-fb;f_En=En/(K_En+En)
  dD=sigma_D*(w_E_over_M*f_En+(1-w_E_over_M)*f_PC)-tau_D*D;dS=sigma_S*Tr-tau_S*S
  dT0=-theta_treg*tau_t_a*(M1*C_rat)*T0+tau_t_d*Tr;dTr=theta_treg*tau_t_a*(M1*C_rat)*T0-tau_t_d*Tr;dR=beta_ros*M1-mu_R*R
  dPL=PL*r_L_P*(1-PL/K_L_P)-alpha_AMP*D*PL
  dPLP=PL*gamma_B_P*(1-Eh/K_E)*(1-PLP/K_LP_P)+PLP*r_LP_P*(1-PLP/K_LP_P)-PLP*(alpha_0*M0+alpha_1_e*M1+alpha_1_r*R)
  dCLP=CL*gamma_B_C*(1-Eh/K_E)*(1-CLP/K_LP_C)-CLP*(alpha_0*M0+alpha_1_e*M1+alpha_1_r*R)
  list(c(dPL,dPLP,dCLP,dEh,dEa,dEn,dM0,dM1,dR,dT0,dTr,dD,dS))}) }

mk_seed = function(e0, pathogenic, plp0) {
  y = setNames(numeric(13), STATE_ORDER)
  y["Eh"] = e0; y["En"] = K_E - e0; y["M0"] = M_TOT; y["T0"] = T_TOT
  if (pathogenic) { y["PL"] = 500; y["PLP"] = plp0 }
  y }

snap_fp = function(endpt, fps) {
  if (is.null(fps) || !nrow(fps) || is.na(endpt["Eh"])) return(NA_integer_)
  vars = intersect(FP_MATCH_VARS, names(fps)); ev = endpt[vars]
  dd = vapply(seq_len(nrow(fps)), function(i) {
    fv = unlist(fps[i, vars]); sc = abs(fv) + 1
    sum(((ev - fv)/sc)^2, na.rm = TRUE) }, numeric(1))
  if (!length(dd) || all(!is.finite(dd))) return(NA_integer_)
  which.min(dd) }

## full 13-state vector of a fixed point (M0, T0 slaved by conservation)
fp_full_state = function(f) {
  s = setNames(rep(NA_real_, length(STATE_ORDER)), STATE_ORDER)
  for (v in intersect(names(f), STATE_ORDER)) s[v] = f[[v]]
  s["M0"] = M_TOT - if (is.null(f$M1) || is.na(f$M1)) 0 else f$M1
  s["T0"] = T_TOT - if (is.null(f$Tr) || is.na(f$Tr)) 0 else f$Tr
  s }

## ============================================================================
## grid_info(): the GRID-seeded finder -- what EXISTS, saddles included.
## Runs FIRST and INDEPENDENTLY of the ray, because it answers an existence
## question that needs no trajectories. In the previous version this block sat
## at the END of ray_one, AFTER the early return for "endpoint_finder_empty" --
## so on the 6.6% of combos where endpoint seeding failed, the grid pass never
## ran and n_field_stable came back NA by construction. That made it impossible
## to tell whether those combos genuinely have no fixed points or whether only
## the seeding failed. Now the grid result is attached to EVERY return path.
## ============================================================================
grid_info = function(par, KG, BETA, RLP, theta) {
  out = list(n_field_stable = NA_integer_, n_field_saddle = NA_integer_,
             lo_field = NA_real_, hi_field = NA_real_, sad_field = NA_real_,
             sad_state = NULL, eig_sad = NA_real_, fp_all = NULL)
  if (GRID_FP != 1) return(out)
  fpg = tryCatch(as.data.frame(find_fixed_points(par, Kg=KG, beta_ros=BETA, r_L_P=RLP,
          M_tot=M_TOT, T_tot=T_TOT, theta_treg=theta, grid_seed=TRUE)),
         error = function(e) NULL)
  if (is.null(fpg) || !nrow(fpg)) return(out)
  fpg$cmp = sqrt((K_E - fpg$Eh)^2 + fpg$PLP^2)
  stg = fpg[fpg$type == "stable", , drop = FALSE]
  sdg = fpg[fpg$type != "stable", , drop = FALSE]
  out$n_field_stable = nrow(stg); out$n_field_saddle = nrow(sdg)
  if (nrow(stg)) { out$lo_field = min(stg$cmp); out$hi_field = max(stg$cmp) }
  ## the saddle BETWEEN the two attractors is the separatrix landmark
  if (nrow(sdg) && nrow(stg) > 1) {
    inb = sdg[sdg$cmp > min(stg$cmp) & sdg$cmp < max(stg$cmp), , drop = FALSE]
    if (nrow(inb)) {
      w = which.min(abs(inb$cmp - (min(stg$cmp) + max(stg$cmp))/2))
      out$sad_field = inb$cmp[w]
      out$sad_state = fp_full_state(inb[w, ])
      out$eig_sad   = inb$max_Re_eig[w] } }
  out$fp_all = fpg          # kept for the endpoint-seeding fallback below
  out }

## ============================================================================
## one condition: integrate the whole ray, seed the finder from the endpoints,
## map every Eh0 to the fixed point it reaches.
## ============================================================================
## `gi` is the grid_info() result for THIS Treg condition, computed ONCE per
## combo by the caller. It must not be computed in here: fixed points do not
## depend on initial conditions, so it is identical across PLP0 levels, and
## calling it per level ran the grid finder 8 times per pathogenic combo for one
## distinct answer (~4x the whole per-combo cost).
ray_one = function(P, par, KG, BETA, RLP, theta, pathogenic, plp0, gi) {
  ## ---- integrate the Eh0 ray --------------------------------------------
  ## ONE integration, factored out: the probe below and the bisection further
  ## down must use identical settings, or the threshold would be located under
  ## dynamics the labels were not computed with.
  one_end = function(e0) {
    o = tryCatch(ode(mk_seed(e0, pathogenic, plp0), c(0, TMAX), rhs_app, P,
                     ## maxsteps scales WITH TMAX. At 8000 steps and TMAX = 4000 a
                     ## stiff ray was already close to the cap; at TMAX = 20000 it
                     ## would silently truncate and return a non-converged state,
                     ## which is the exact failure the longer horizon is meant to
                     ## fix. 200000 is what the local timing runs used.
                     ##
                     ## hmax IS NOT OPTIONAL. deSolve defaults it to
                     ## max(diff(times)), which for a two-point output grid is
                     ## TMAX itself -- i.e. no cap at all. lsoda then strides the
                     ## stiff tail, the local error blows up in every state, and
                     ## PLP crosses zero where the logistic term is quadratic and
                     ## negative, so it escapes to -Inf. That is where dev_max =
                     ## 1.3e154 came from. Measured on row_ind 1772 / Kg 30 /
                     ## beta_ros 0.0045 / r_L_P 0.08 / theta 1 / Eh0 15:
                     ##   hmax Inf,200,100,50,20 -> NaN at t = 1381
                     ##   hmax 10               -> PLP = -1.5e156, still garbage
                     ##   hmax 5                -> t = 20000, min PLP over the
                     ##                            whole path +5.5e-21, chronic
                     ##                            at Eh = 77.40
                     ## The failure is NOT always loud: at Eh0 = 35 the unbounded
                     ## run returns finite numbers with Eh = 101.3, i.e. above
                     ## K_E, and at Eh0 = 20 it silently picks a different basin.
                     ## Clamping the state does not fix this -- it removes the
                     ## symptom in PLP while every other state is still being
                     ## integrated with the same oversized steps.
                     ##
                     ## lsodar, NOT lsoda, and the extinction root from
                     ## lib/skeleton_defs.R. Without it a pathogen driven to
                     ## 1e-100 can regrow 100 decades once M1 is suppressed, and
                     ## which way it goes is decided by the sign of the
                     ## integrator's rounding error. The event sets PL and PLP to
                     ## exactly zero when BOTH are below EXT_TOL, from which point
                     ## both derivatives are zero. It adds an output row at the
                     ## root time; nothing here or downstream assumes a two-row
                     ## result -- this reads .N, and fp_settle_stable() reads
                     ## nrow(out) -- but check any new consumer that does.
                     ##
                     ## atol/rtol are explicit now. atol must sit well below
                     ## EXT_TOL or the root fires on noise; see the note at the
                     ## ATOL definition above.
                     method = "lsodar", rootfun = rootfun_ext, events = EVENTS_EXT,
                     maxsteps = MAXSTEPS, hmax = HMAX, atol = ATOL, rtol = RTOL),
                 error = function(e) NULL, warning = function(w) NULL)
    ## A DROPPED WARNING IS NOT A SUCCESS. This handler swallows warnings as well
    ## as errors, so a run that hit maxsteps and was truncated used to be
    ## indistinguishable from a converged one. Re-check the horizon explicitly:
    ## the last output time must be TMAX, or the integration did not finish.
    if (is.null(o) || nrow(o) < 2) return(NULL)
    if (!isTRUE(all.equal(unname(o[nrow(o), 1]), TMAX))) return(NULL)
    unlist(as.data.table(o)[.N, ..STATE_ORDER]) }

  ## ---- STEP 1: THE PROBE -------------------------------------------------
  E = lapply(EH0_SCAN, one_end)
  ok = !vapply(E, is.null, logical(1))
  if (!any(ok)) return(c(list(note = "all_integrations_failed"), gi))

  ## ---- ENDPOINT-seeded finder: the ray-reachable attractors -------------
  sm = do.call(rbind, lapply(E[ok], function(z)
        c(M1 = unname(z["M1"]), PLP = unname(z["PLP"]), L = unname(z["Ea"] + z["En"]))))
  sm = unique(sm[is.finite(sm[,"M1"]) & is.finite(sm[,"PLP"]) & is.finite(sm[,"L"]), , drop = FALSE])
  fpr = if (!nrow(sm)) NULL else
    tryCatch(as.data.frame(find_fixed_points(par, Kg=KG, beta_ros=BETA, r_L_P=RLP,
              M_tot=M_TOT, T_tot=T_TOT, theta_treg=theta,
              extra_seeds=sm, grid_seed=FALSE)), error = function(e) NULL)
  ## ---- FALLBACK when endpoint seeding returns nothing --------------------
  ## Measured on batch 1 (21_empty_probe.R, 300 combos): of the combos where the
  ## endpoint-seeded finder came back empty, 100% of PATHOGENIC ones have >=1
  ## stable fixed point that the grid pass finds, and 0% of STERILE ones do.
  ## So on the pathogenic side the attractors were always there and only the
  ## seeding failed -- 6.5% of rows were being discarded as missing when they
  ## were merely unseeded, and not uniformly: the loss rate runs 15% at
  ## r_L_P = 0.04 down to 2.7% at 0.48, so it biased every pooled figure.
  ## Snapping endpoints against the FIELD set instead is safe: an off-ray fixed
  ## point simply never becomes any endpoint's nearest match, so this recovers
  ## the ray data without letting unreachable states into it. seed_mode records
  ## which path was taken so the two can always be separated downstream.
  seed_mode = "endpoint"
  if (is.null(fpr) || !nrow(fpr)) {
    if (!is.null(gi$fp_all) && nrow(gi$fp_all)) {
      fpr = gi$fp_all; seed_mode = "grid_fallback"
    } else return(c(list(note = "endpoint_finder_empty", seed_mode = "none"),
                    gi[setdiff(names(gi), "fp_all")]))
  }
  fpr$cmp = sqrt((K_E - fpr$Eh)^2 + fpr$PLP^2)
  fpr = fpr[order(fpr$cmp), , drop = FALSE]           # id 1 = lowest burden
  rownames(fpr) = NULL

  ## ---- STEP 3: THREE LABELS ---------------------------------------------
  ## snap_fp returns the nearest fixed point of ANY type, so `stab` is not
  ## optional bookkeeping -- it is what separates "reached an attractor" from
  ## "came to rest on something that is not one".
  idx  = rep(NA_integer_, length(EH0_SCAN))
  for (k in which(ok)) idx[k] = snap_fp(E[[k]], fpr)
  stab = !is.na(idx) & fpr$type[idx] %in% "stable"
  ## the IC -> FP map, one character per probed Eh0. '.' = no stable match.
  map = paste(ifelse(stab, as.character(idx), "."), collapse = "")

  ## CLUSTER FIRST, then label. Two stable fixed points that fp_same() calls one
  ## attractor are one attractor for every purpose here, so labelling by bare
  ## index would let a pair inside one cluster manufacture a crossing. lo_ids is
  ## every reached fixed point in the LOWEST-burden cluster; the label is
  ## membership in it, i.e. attractor identity, never a cutoff on a value.
  ##
  ## The grouping is Chebyshev on (Eh, PLP), not a gap in the composite -- see the
  ## RAY_CLUST_TOL note above for the measured case the composite could not
  ## separate, and for why the tolerance is 1 rather than the composite's 2. fp_cluster() returns cluster ids ordered
  ## by increasing composite, so cluster 1 is still the lowest-burden attractor
  ## and everything downstream that reads lo/hi keeps its meaning.
  reached = sort(unique(idx[stab]))
  n_ray = 0L; lo_ids = integer(0); lo_i = NA_integer_; hi_i = NA_integer_
  if (length(reached)) {
    .cl   = fp_cluster(fpr$Eh[reached], fpr$PLP[reached], fpr$cmp[reached])
    n_ray = length(unique(.cl))
    lo_ids = reached[.cl == 1L]
    ## representative of a cluster = its lowest-composite member, so lo_i and
    ## hi_i remain single fixed points and the landmark states stay well defined.
    .rep  = function(k) { m = reached[.cl == k]; m[which.min(fpr$cmp[m])] }
    lo_i  = .rep(1L)
    hi_i  = if (n_ray >= 2L) .rep(n_ray) else NA_integer_
  }

  gl = rep(NA_integer_, length(EH0_SCAN))
  gl[stab] = as.integer(idx[stab] %in% lo_ids)
  g = EH0_SCAN[!is.na(gl)]; v = gl[!is.na(gl)]

  ## ---- STEP 4: CROSSINGS, THEN BISECT THE ONE THAT EXISTS ----------------
  tr = if (length(v) > 1) which(v[-1] != v[-length(v)]) else integer(0)
  ## coarse crossing locations, kept for continuity: the midpoint of each
  ## bracketing pair. Eh0_star below supersedes them when there is exactly one.
  cross = if (length(tr)) (g[tr] + g[tr + 1]) / 2 else numeric(0)
  note = if (!length(v)) "no_labels"
    else if (all(v == 1L)) "lo_basin_throughout"
    else if (all(v == 0L)) "hi_basin_throughout"
    else if (length(tr) == 1) "ok"
    else sprintf("non_monotonic_%d", length(tr))

  ## The label of a candidate Eh0, on exactly the rule used for the probe.
  lab_of = function(e0) {
    z = one_end(e0); if (is.null(z)) return(NA)
    i = snap_fp(z, fpr); if (is.na(i)) return(NA)
    if (!identical(fpr$type[i], "stable")) return(NA)
    i %in% lo_ids
  }
  Eh0_star = NA_real_; dir = NA_character_; n_bisect = 0L
  if (length(tr) == 1L) {
    lo = g[tr]; hi = g[tr + 1]
    ## which side of the crossing holds the low-burden basin. Read, not assumed:
    ## an arm can cross going up in Eh0 or going down, and both occur.
    up = (v[tr + 1] == 1L)
    ok_bis = TRUE
    while (hi - lo > EH0_TOL) {
      m = (lo + hi) / 2; r = lab_of(m); n_bisect = n_bisect + 1L
      if (is.na(r)) { ok_bis = FALSE; break }
      if (r == up) hi = m else lo = m
    }
    if (ok_bis) {
      Eh0_star = (lo + hi) / 2
      dir = if (up) "increasing" else "decreasing"
    } else note = "bisect_unlabelled"
  }

  ## ---- RELATIVE deviation: is each endpoint actually AT its matched FP? --
  ## snap_fp assigns by nearest neighbour, which always returns something. This
  ## is the check that the assignment means anything. Relative, because PL ~ 1000
  ## and R ~ 1 and an absolute norm would only ever report PL.
  ## fpr carries only the states the finder solves for; fp_full_state() derives
  ## M0/T0 from the conservation laws and pads the rest with NA. Indexing fpr by
  ## all 13 names directly is an "undefined columns selected" error.
  dev_states = if (is.null(DEV_STATES)) STATE_ORDER else DEV_STATES
  devs = vapply(which(ok & stab), function(k) {
    fpv = fp_full_state(fpr[idx[k], ])[dev_states]
    ev  = E[[k]][dev_states]
    keep = is.finite(fpv) & is.finite(ev)
    if (!any(keep)) NA_real_ else max(abs(ev[keep] - fpv[keep]) / (abs(fpv[keep]) + 1))
  }, numeric(1))
  dev_max = if (length(devs)) max(devs, na.rm = TRUE) else NA_real_
  dev_med = if (length(devs)) median(devs, na.rm = TRUE) else NA_real_

  base = list(note = note, map = map, n_ray = n_ray,
       ## the COMPLETE finder output, one vector per state. Length = however many
       ## fixed points were returned; nothing is selected, nothing is dropped.
       fp_n = nrow(fpr), fp_tbl = list(fpr),
       dev_max = dev_max, dev_med = dev_med, dev_n = length(devs),
       n_fail = sum(!ok), n_unmatched = sum(!stab),
       lo_cmp = if (!is.na(lo_i)) fpr$cmp[lo_i] else NA_real_,
       hi_cmp = if (!is.na(hi_i)) fpr$cmp[hi_i] else NA_real_,
       lo_state = if (!is.na(lo_i)) fp_full_state(fpr[lo_i, ]) else NULL,
       hi_state = if (!is.na(hi_i)) fp_full_state(fpr[hi_i, ]) else NULL,
       eig_lo = if (!is.na(lo_i)) fpr$max_Re_eig[lo_i] else NA_real_,
       eig_hi = if (!is.na(hi_i)) fpr$max_Re_eig[hi_i] else NA_real_,
       ## star IS the bisected threshold now, NA unless exactly one crossing was
       ## found. star_coarse keeps the probe-resolution value for comparison; it
       ## is what v5 called star, and the two should agree to within PROBE_BY/2.
       star = Eh0_star, dir = dir, n_bisect = n_bisect, probe_by = PROBE_BY,
       star_coarse = if (!length(cross)) NA_real_ else max(cross),
       n_cross = length(cross), cross = cross, seed_mode = seed_mode)
  c(base, gi[setdiff(names(gi), "fp_all")]) }

## ============================================================================
## one combo, both conditions -> ONE row
## ============================================================================
analyse_combo = function(par, trow, row_ind, KG, BETA, RLP) {
  pathogenic = (RLP > 0)
  mkp = function(th){ p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]]
    p$Kg=KG; p$beta_ros=BETA; p$r_L_P=RLP; p$r_LP_P=RLP; p$theta_treg=th; p }

  ## grid finder ONCE per condition, hoisted out of the PLP0 loop (see ray_one)
  gi0 = grid_info(par, KG, BETA, RLP, 0L)
  gi1 = grid_info(par, KG, BETA, RLP, 1L)

  ## sterile has no pathogen to seed, so its ray does not depend on PLP0 --
  ## run it ONCE and do not multiply the cost by length(PLP0_LEVELS)
  ## v4 looped PLP0 here. One level now -- see the header.
  levels_here = PLP0_MAIN
  out_rows = vector("list", 1L)
  for (li in 1L) {
  plp0 = PLP0_MAIN
  A = ray_one(mkp(0L), par, KG, BETA, RLP, 0L, pathogenic, plp0, gi0)
  B = ray_one(mkp(1L), par, KG, BETA, RLP, 1L, pathogenic, plp0, gi1)
  g = function(z, f, d = NA_real_) { v = z[[f]]; if (is.null(v) || !length(v)) d else v[1] }

  ## ---- rescue_gap -------------------------------------------------------
  ## "How much wider is the set of initial conditions on this line that end in the
  ## LOW-burden basin, with Tregs on?" THE SIGN DEPENDS ON THE DIRECTION, because
  ## the low-basin set is [star, EH0_HI] when the arm runs increasing and
  ## [EH0_LO, star] when it runs decreasing. Using off - on for both reports every
  ## decreasing pair backwards -- measured on the five decreasing combinations in
  ## the v5 pathogenic_help measurable pool, all five read positive under the
  ## unsigned formula and all five are Tregs NARROWING the set.
  ##
  ## Arms running in OPPOSITE directions are not comparable and their difference is
  ## not a rescue; those get NA and their own kind, not a repaired number.
  ## rescue_gap_kind says which case produced the value, so no downstream filter
  ## has to infer it from NA patterns.
  .sa = g(A,"star"); .sb = g(B,"star"); .da = g(A,"dir",""); .db = g(B,"dir","")
  .na = g(A,"note",""); .nb = g(B,"note","")
  gap = NA_real_; gap_kind = "none"
  if (is.finite(.sa) && is.finite(.sb) && nzchar(.da) && identical(.da, .db)) {
    gap = if (identical(.da, "increasing")) .sa - .sb else .sb - .sa
    gap_kind = "measured"
  } else if (is.finite(.sa) && identical(.nb, "lo_basin_throughout")) {
    ## ON reaches the low basin from every probed Eh0: the whole line is rescued.
    ## Exact on the line's domain, not a bound, but a different quantity from a
    ## threshold difference -- hence its own kind.
    gap = if (identical(.da, "increasing")) .sa - EH0_LO else EH0_HI - .sa
    gap_kind = "on_lo_basin_throughout"
  } else if (is.finite(.sa) && is.finite(.sb)) {
    gap_kind = "opposite_directions"
  }

  rec = data.table(
    row_ind = row_ind, Kg = KG, beta_ros = BETA, r_L_P = RLP,
    ## THE GLOBAL KEY, CARRIED NOT RECONSTRUCTED. row_ind restarts at 1 in every
    ## tau grid, so it identifies nothing on its own; batch and tau_id are stamped
    ## into the grid by 00_parameter_grids.R and travel with the parameter row. They
    ## used to be rebuilt in 02_assemble.R and then inverted by arithmetic in
    ## three downstream scripts against a stride that was nominally a variable.
    batch  = as.integer(trow[["batch"]]),
    tau_id = as.integer(trow[["tau_id"]]),
    PLP0 = if (pathogenic) plp0 else 0,
    ## ---- THE SETTINGS THIS ROW WAS COMPUTED UNDER ------------------------
    ## Written per row, not assumed. ray_clust_tol in particular decides how many
    ## attractors this row has, and 02_assemble.R re-clusters the full
    ## fixed-point set at MERGE_TOL -- if the two differ, the labels and the
    ## fixed-point set describe different clusterings. That used to be a comment;
    ## carrying the value makes it an assertion downstream.
    ray_clust_tol = RAY_CLUST_TOL, ext_tol = EXT_TOL,
    atol = ATOL, rtol = RTOL, hmax = HMAX, tmax = TMAX, maxsteps = MAXSTEPS,
    ## ---- ray-reachable (endpoint-seeded finder) --------------------------
    n_ray_off = g(A,"n_ray"), n_ray_on = g(B,"n_ray"),
    lo_off = g(A,"lo_cmp"), hi_off = g(A,"hi_cmp"),
    lo_on  = g(B,"lo_cmp"), hi_on  = g(B,"hi_cmp"),
    eig_lo_off = g(A,"eig_lo"), eig_hi_off = g(A,"eig_hi"),
    eig_lo_on  = g(B,"eig_lo"), eig_hi_on  = g(B,"eig_hi"),
    ## Eh0_off / Eh0_on are the BISECTED thresholds, to EH0_TOL. NA unless that arm
    ## had exactly one crossing. *_coarse is the probe-resolution value for
    ## comparison; dir says which way the arm runs; n_cross is what decides whether
    ## a threshold exists at all; probe_by records the resolution the crossing count
    ## was established at, so a zero never reads as proof on its own.
    Eh0_off = g(A,"star"), Eh0_on = g(B,"star"),
    dir_off = g(A,"dir",""), dir_on = g(B,"dir",""),
    Eh0_off_coarse = g(A,"star_coarse"), Eh0_on_coarse = g(B,"star_coarse"),
    n_bisect_off = g(A,"n_bisect"), n_bisect_on = g(B,"n_bisect"),
    probe_by = PROBE_BY, eh0_tol = EH0_TOL,
    n_cross_off = g(A,"n_cross"), n_cross_on = g(B,"n_cross"),
    cross_off = paste(A$cross, collapse=";"), cross_on = paste(B$cross, collapse=";"),
    ## ---- existence (grid-seeded finder; saddles) --------------------------
    n_field_stable_off = g(A,"n_field_stable"), n_field_stable_on = g(B,"n_field_stable"),
    n_field_saddle_off = g(A,"n_field_saddle"), n_field_saddle_on = g(B,"n_field_saddle"),
    lo_field_off = g(A,"lo_field"), hi_field_off = g(A,"hi_field"), sad_field_off = g(A,"sad_field"),
    lo_field_on  = g(B,"lo_field"), hi_field_on  = g(B,"hi_field"), sad_field_on  = g(B,"sad_field"),
    eig_sad_off = g(A,"eig_sad"), eig_sad_on = g(B,"eig_sad"),
    ## ---- derived ----------------------------------------------------------
    rescue_gap = gap, rescue_gap_kind = gap_kind,
    hi_lo_gap_off = g(A,"hi_cmp") - g(A,"lo_cmp"),
    ## ---- the IC -> FP map, one char per probed Eh0 in EH0_SCAN -------------
    map_off = if (is.null(A$map)) NA_character_ else A$map,
    map_on  = if (is.null(B$map)) NA_character_ else B$map,
    seed_mode_off = g(A,"seed_mode",""), seed_mode_on = g(B,"seed_mode",""),
    n_fail_off = g(A,"n_fail"), n_fail_on = g(B,"n_fail"),
    n_unmatched_off = g(A,"n_unmatched"), n_unmatched_on = g(B,"n_unmatched"),
    note_off = g(A,"note",""), note_on = g(B,"note",""))

  ## full 13-state of every landmark, both conditions (kept for continuity with
  ## v4-era scripts; the authoritative record is the fp_* list-columns below)
  add = function(st, pre) for (v in STATE_ORDER)
    rec[[paste0(pre, "_", v)]] <<- if (is.null(st)) NA_real_ else unname(st[v])
  add(A$lo_state, "lo_off"); add(A$hi_state, "hi_off"); add(A$sad_state, "sad_off")
  add(B$lo_state, "lo_on");  add(B$hi_state, "hi_on");  add(B$sad_state, "sad_on")

  ## ---- THE FIXED POINT SET, as list-columns -------------------------------
  ## One numeric vector per (condition, state), length = number of fixed points
  ## the finder returned for that condition. map_<cond> indexes into these, so
  ## every probed Eh0 resolves to real coordinates -- which is exactly what v4
  ## could not do once RAY_CLUST_TOL nulled a landmark.
  addfp = function(z, pre) {
    tb = if (is.null(z$fp_tbl)) NULL else z$fp_tbl[[1]]
    ## same reason as above: go through fp_full_state() row by row rather than
    ## reading STATE_ORDER off fpr, which does not have all of them.
    FS = if (is.null(tb) || !nrow(tb)) NULL else
         do.call(rbind, lapply(seq_len(nrow(tb)), function(i) fp_full_state(tb[i, ])))
    for (v in STATE_ORDER)
      rec[[paste0("fp_", pre, "_", v)]] <<-
        list(if (is.null(FS)) numeric(0) else as.numeric(FS[, v]))
    for (v in c("cmp", "max_Re_eig"))
      rec[[paste0("fp_", pre, "_", v)]] <<-
        list(if (is.null(tb) || !nrow(tb) || !v %in% names(tb)) numeric(0) else as.numeric(tb[[v]]))
    rec[[paste0("fp_", pre, "_type")]] <<-
      list(if (is.null(tb) || !nrow(tb)) character(0) else as.character(tb$type))
    rec[[paste0("fp_", pre, "_n")]] <<- if (is.null(tb)) 0L else nrow(tb)
  }
  addfp(A, "off"); addfp(B, "on")
  for (nm in c("dev_max","dev_med","dev_n")) {
    rec[[paste0(nm, "_off")]] = g(A, nm); rec[[paste0(nm, "_on")]] = g(B, nm) }
  out_rows[[li]] = rec
  }
  rbindlist(out_rows, fill = TRUE) }

## ============================================================================
## driver
## ============================================================================
args = commandArgs(trailingOnly = TRUE)
n1 = as.integer(ifelse(length(args) >= 1, args[1], 1))
n2 = as.integer(ifelse(length(args) >= 2, args[2], 1))
n3 = as.integer(ifelse(length(args) >= 3, args[3], 1))

tau_grid = readRDS(paste0("./data/tau_grid_v6_", n3, ".rds"))
## THE GRID MUST CARRY THE KEY. 00_parameter_grids.R stamps batch, tau_id and
## tau_stride into every grid it writes, and the sweep passes them straight
## through so that nothing downstream reconstructs them by arithmetic. A grid
## written before that change would produce rows with an NA key that only fails
## much later, in a merge, so it is rejected here.
if (!all(c("batch","tau_id") %in% names(tau_grid)))
  stop("tau_grid_v6_", n3, ".rds predates the stamped key -- rerun 00_parameter_grids.R")
stopifnot(all(tau_grid$batch == as.integer(n3)), !anyDuplicated(tau_grid$tau_id))

## ----------------------------------------------------------------------------
## CAN THESE PARAMETERS EVEN PRODUCE DAMAGE?  Reported before any work is done.
## The sterile branch reduces to one dimensionless group,
##      Lambda  = tau_r * beta_ros * K_M / (mu_R * kappa_E)
##      Eh_min  = K_E / (1 + Lambda)
## Eh_min is the LOWEST epithelium a (row, beta_ros) pair can reach even with
## every macrophage inflamed. A row with Eh_min = 89 will sit at high Eh with
## M1 = 92 and look like a discovery; it is a pinned parameter set. This block
## exists because that went unnoticed until it turned up as a "fully inflamed,
## intact barrier" fixed point in the figures.
## ----------------------------------------------------------------------------
{
  .g = expand.grid(r = seq_len(nrow(tau_grid)), b = BETA_ROS_VALS)
  .lam = tau_grid$tau_r[.g$r] * .g$b * FP_FIXED$K_M /
         (tau_grid$mu_R[.g$r] * tau_grid$kappa_E[.g$r])
  .ehm = FP_FIXED$K_E / (1 + .lam)
  cat(sprintf("[Lambda] tau_r*beta*K_M/(mu_R*kappa_E): p05 %.2f | med %.2f | p95 %.1f
",
              quantile(.lam, .05), median(.lam), quantile(.lam, .95)))
  cat(sprintf("[Lambda] Eh_min: cannot drop below 80 in %.2f%% of (row x beta) | below 50 in %.2f%% | collapses (<5) in %.2f%%
",
              100*mean(.ehm > 80), 100*mean(.ehm > 50), 100*mean(.ehm < 5)))
  if (mean(.ehm > 80) > 0.02)
    cat("[Lambda] *** WARNING: >2% of the grid CANNOT produce damage at any M1.
",
        "          Those rows will return high-Eh high-M1 equilibria that are
",
        "          pinned by parameters, not dynamics. Raise the beta_ros floor.
", sep="")
}
split_equal = function(x, k) if (k == 1) list(`1` = x) else split(x, cut(seq_along(x), k, labels = FALSE))
## ROW_N subsamples the FULL row set BEFORE chunking, so ROW_N is the total
## number of tau rows the whole array covers. (Applying it after the split would
## subsample within a chunk, where ROW_N is larger than the chunk and the call
## would silently do nothing.) tau rows are an unstructured parameter sample, so
## a random subset is unbiased -- every proportion is the same population
## quantity with a wider interval. The seed is fixed so all array tasks draw
## the SAME subset and then partition it.
all_rows = seq_len(nrow(tau_grid))
## CELL-LIST MODE (see the header): only the tau rows of this batch that carry a
## listed cell. ROW_N does not apply here -- the list already says which rows.
CELLS = Sys.getenv("CELLS", "")
.TC = NULL
if (nzchar(CELLS)) {
  .TC = as.data.table(readRDS(CELLS))
  if (!all(c("batch", "row_ind", "Kg", "r_L_P") %in% names(.TC)))
    stop("CELLS file ", CELLS, " needs columns batch, row_ind, Kg, r_L_P")
  all_rows = sort(unique(.TC[batch == n3, row_ind]))
  if (!length(all_rows)) stop("[cells] no listed cell in batch ", n3)
  cat(sprintf("[cells] %s: batch %d, %d tau rows, %d cells x %d beta_ros levels\n", CELLS, n3,
              length(all_rows), .TC[batch == n3, .N], length(BETA_ROS_VALS)))
} else if (ROW_N > 0 && ROW_N < length(all_rows)) {
  set.seed(ROW_SEED); all_rows = sort(sample(all_rows, ROW_N)) }
loop_over_tau = split_equal(all_rows, n1)[[n2]]

combos_grid = CJ(Kg = KG_VALS, beta_ros = BETA_ROS_VALS, r_L_P = RLP_VALS, sorted = FALSE)
## the combinations of one tau row: the full block, or that row's listed cells at
## every BETA_ROS_VALS level
combos_of = function(ri) {
  if (is.null(.TC)) return(combos_grid)
  cr = .TC[batch == n3 & row_ind == ri, .(Kg, r_L_P)]
  CJ(k = seq_len(nrow(cr)), beta_ros = BETA_ROS_VALS)[, .(Kg = cr$Kg[k], beta_ros, r_L_P = cr$r_L_P[k])]
}
combos = combos_grid
dir_name = paste0(OUT_PREFIX, n3)
dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Output: %s | rows: %d | combos/row: %s | probe: %d values | PLP0: %s | grid FP: %d\n",
            dir_name, length(loop_over_tau),
            if (is.null(.TC)) as.character(nrow(combos)) else "listed cells x beta_ros levels", length(EH0_SCAN),
            PLP0_MAIN, GRID_FP))

## ---- TIMEIT: sample a few combos, project the row and batch cost -----------
## A whole row is ~845 combos and hours long, so never time one locally. This
## runs TIME_N combos spread across the grid (so the sample spans sterile vs
## pathogenic and the Kg range, which is where the cost varies) and multiplies.
if (as.integer(Sys.getenv("TIMEIT", "0")) == 1) {
  TIME_N = as.integer(Sys.getenv("TIME_N", "5"))
  row_ind = loop_over_tau[1]
  trow = tau_grid[row_ind, ]
  combos = combos_of(row_ind)
  par = as.list(trow); par$row_ind = NULL; par$treg_eff = NULL; par$phi_instability = NULL
  ## the stamped key travels in `trow` and is written into the output row; it is
  ## not a model parameter and must not reach the RHS closure or the finder.
  par$batch = NULL; par$tau_id = NULL; par$tau_stride = NULL
  ## RANDOM sample, not seq(): the grid has 169 combos per Kg block, so an
  ## evenly-strided index walks the diagonal (Kg, beta and r_L_P all rising
  ## together) and samples one corner of the grid instead of spanning it.
  set.seed(as.integer(Sys.getenv("TIME_SEED", "1")))
  pick = sort(sample(nrow(combos), min(TIME_N, nrow(combos))))
  tt = sapply(pick, function(i) {
    t0 = Sys.time()
    ## print the error instead of returning NULL: a silent FAILED tells you nothing
    d = tryCatch(analyse_combo(par, trow, row_ind, combos$Kg[i], combos$beta_ros[i], combos$r_L_P[i]),
                 error = function(e) { cat("    ERROR: ", conditionMessage(e), "\n", sep=""); NULL })
    dt = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  Kg=%2g beta=%.2f rP=%.2f -> %s | n_ray off/on %s/%s | %.1f s\n",
                combos$Kg[i], combos$beta_ros[i], combos$r_L_P[i],
                if (is.null(d)) "FAILED" else "ok",
                if (is.null(d)) "-" else d$n_ray_off[1], if (is.null(d)) "-" else d$n_ray_on[1], dt))
    dt })
  per = mean(tt)
  cat(sprintf("\nmean %.1f s per combo (both conditions, %d probe values + bisection)\n", per, length(EH0_SCAN)))
  cat(sprintf("one tau row (%d combos): %.1f h\n", nrow(combos), per * nrow(combos) / 3600))
  cat(sprintf("a 1000-task x 24 h batch (24000 core-h) buys ~%d tau rows\n",
              round(24000 / (per * nrow(combos) / 3600))))
  quit(save = "no")
}

options(warn = -1)
acc = list(); flush_counter = 0
for (row_ind in loop_over_tau) {
  t0 = Sys.time()
  trow = tau_grid[row_ind, ]
  par = as.list(trow); par$row_ind = NULL; par$treg_eff = NULL; par$phi_instability = NULL
  ## the stamped key travels in `trow` and is written into the output row; it is
  ## not a model parameter and must not reach the RHS closure or the finder.
  par$batch = NULL; par$tau_id = NULL; par$tau_stride = NULL
  combos = combos_of(row_ind)
  rows = vector("list", nrow(combos))
  for (i in seq_len(nrow(combos)))
    rows[[i]] = tryCatch(analyse_combo(par, trow, row_ind,
                            combos$Kg[i], combos$beta_ros[i], combos$r_L_P[i]),
                         error = function(e) NULL)
  keep = !vapply(rows, is.null, logical(1))
  acc = c(acc, list(rbindlist(rows[keep], fill = TRUE)))
  cat(sprintf("[%s] row %d: %d/%d combos, %.2f min\n", Sys.time(), row_ind,
              sum(keep), nrow(combos), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  ## A tau row now costs hours, so flushing every 4 rows risked losing most of
  ## a day's work on a kill. Default is every row: the write is trivial next to
  ## the integration and a dead task then loses at most one row.
  flush_counter = flush_counter + 1
  if (flush_counter %% FLUSH_EVERY == 0) {
    ## .rds, not .csv: list-columns cannot survive a CSV round trip
    saveRDS(rbindlist(acc, fill = TRUE),
            file.path(dir_name, sprintf("fp_v6_flush_%d_chunk%d.rds", flush_counter, n2)))
    acc = list() } }
if (length(acc)) saveRDS(rbindlist(acc, fill = TRUE),
                         file.path(dir_name, sprintf("fp_v6_final_chunk%d_counter%d.rds", n2, flush_counter)))
options(warn = 0)
cat("\nAll done!\n")
