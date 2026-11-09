#!/bin/bash
#SBATCH --job-name=fp_sweep
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --array=0-399
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs_sweep/fp_%x_%A_%a.out
#SBATCH --error=logs_sweep/fp_%x_%A_%a.err

## ---------------------------------------------------------------------------
## sweep_array.sh -- the parameter sweep (pipeline/01_sweep.R) as a Slurm array,
## ONE file for every tau batch.
##
##   BATCH=1 sbatch pipeline/cluster/sweep_array.sh   (or submit_sweep.sh for all)
##
## v5 kept a separate copy of this script per batch, differing only in BATCH and
## the log directory. Two of those copies had already drifted apart on --time.
## BATCH is an environment variable now; sbatch propagates the submitting
## environment by default, so nothing else is needed.
##
## WHAT ONE TASK DOES
## One array task takes a disjoint slice of the tau rows in tau_grid_v6_<BATCH>.rds
## and, for every parameter combination in that slice, integrates the Eh0 line
## under both Treg conditions, locates the fixed points, labels each probed Eh0,
## counts the crossings, and bisects the one crossing if there is exactly one.
##
## Grid per tau row: 441 = 7 Kg (15,20,25,30,35,40,45)
##                       x 9 beta_ros (0.0015..0.0095 step 0.001)
##                       x 7 r_L_P (0..0.48 step 0.08),  at a single PLP0 = 10.
## Previously 315 = 5 Kg (20,30,40,50,60). See the KG_VALS note in pipeline/01_sweep.R
## for why the levels moved: of 67,486 gated rows the old grid put 22,929 at
## Kg 20, 25,557 at 30, 15,862 at 40, 3,057 at 50 and 81 at 60, so the resolution
## now sits where the population is.
##
## COST, AND WHY THE OLD SIZING NO LONGER APPLIES. Two changes push in the same
## direction and NEITHER has been measured on a Della core:
##   1. 441 combinations against 315 is +40% by itself.
##   2. atol/rtol are now set explicitly at 1e-14 / 1e-8 instead of deSolve's
##      1e-6 default. atol has to sit well below the extinction floor EXT_TOL or
##      the root fires on integrator noise; it is not optional, and a tighter
##      tolerance means more internal steps. The size of that is unknown.
## The old measurement -- 19.3 s per combo x 315 = 1.7 h per tau row -- is now a
## LOWER BOUND, not an estimate. The wall clock is raised from 8 h to 24 h to
## cover it, and the timing run below is mandatory rather than advisory.
##
## RE-MEASURE ON DELLA BEFORE SUBMITTING THE FULL SET. Not optional this time:
##   BATCH=1 TIMEIT=1 TIME_N=5 Rscript pipeline/01_sweep.R 400 1 1
## and read the printed hours per tau row. If it exceeds ~4 h, raise --time again
## or cut N_ROW_CHUNKS; do not submit on the old numbers.
##
## SIZING. 400 tasks, 2000 tau rows per batch => 5 rows per task.
##   At the old 1.7 h per tau row that was 8.5 h per task; at +40% and an unknown
##   tolerance cost it is at least 12 h, which is why --time is 24 h.
## An over-request costs queue priority; an under-request costs data, because
## FLUSH_EVERY=1 writes each completed row as it goes, so a task that hits the
## wall part way through leaves a run that LOOKS successful with rows quietly
## missing. Ask for more time than the estimate.
##
## CRITICAL: the --array range and N_ROW_CHUNKS must match.
##     --array=0-(K-1)   <=>   N_ROW_CHUNKS=K
## ---------------------------------------------------------------------------
## Conda's activation scripts read $PS1, which is unset in a non-interactive
## shell. Under `set -u` that is a fatal "PS1: unbound variable" before anything
## runs. Give PS1 a value and drop -u across the conda block only; strictness
## resumes immediately after. Same fix as FP_night_tzero.sh.
export PS1=${PS1:-}
set +u
module load anaconda3/2023.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_Treg
set -u

## the R scripts set their working directory to SRC and use relative paths
## (./lib/..., ./data/tau_grid_v6_<n>.rds). SRC is the repository root on the cluster.
export SRC=${SRC:-$HOME/tregs_ode_paper}
cd "$SRC" || exit 1
mkdir -p logs_sweep

BATCH=${BATCH:-1}                           # which tau_grid_v6_<BATCH>.rds (1..5)
N_ROW_CHUNKS=${N_ROW_CHUNKS:-400}           # = array size above; keep in sync!
ROW_CHUNK_ID=$(( SLURM_ARRAY_TASK_ID + 1 ))

export ROW_N=${ROW_N:-0}          # 0 = all rows in the tau file (2000 per batch)
export ROW_SEED=${ROW_SEED:-20260808}  # same seed in every task -> same subset, partitioned

## TMAX = 12000, derived rather than guessed. A mode decays as exp(-|Re(lam)|t)
## and the downstream gate is dev_med <= 1e-2, so t = 4.6/|Re(lam)| only ARRIVES
## at the gate. Over the 3,403 eigenvalues in v5 batch 1 the fraction of modes
## still above the gate is 13.811% at 4000, 0.852% at 5000, 0.441% at 8000 and
## 0.000% at 12000. The slowest mode (|Re| = 5.29e-4) sits at a 1.7e-3 residual
## at 12000, ~6x margin; worst deviation vs the TMAX=20000 answer, 2.4e-7.
## With HMAX capped the cost is roughly LINEAR in TMAX, so this is a first-order
## cost driver -- do not raise it casually.
export TMAX=${TMAX:-12000}
export MAXSTEPS=${MAXSTEPS:-200000}   # MUST rise with TMAX or stiff runs truncate
## HMAX. NOT optional. deSolve defaults the max internal step to max(diff(times)),
## which for the two-point output grid used here is TMAX itself, i.e. no cap at
## all. lsoda then strides the stiff tail and returns wrong states, sometimes
## finite and plausible ones: hmax 10 gave PLP = -1.5e156, and 20/50/100/Inf gave
## NaN. 5 was correct on every test case.
export HMAX=${HMAX:-5}

## ERROR TOLERANCES AND THE EXTINCTION FLOOR. These used to be deSolve's defaults
## (atol = rtol = 1e-6) and that is no longer safe. EXT_TOL is one bacterium in
## state units under the single conversion 1 unit = 1e11 organisms, and the
## extinction root tests max(PL, PLP) against it; lsoda accepts a local error of
## rtol*|y| + atol, so at atol = 1e-6 the tolerance would be a hundred thousand
## times coarser than the threshold being tested and the event would fire on
## noise. lib/extinction_floor.R refuses to load unless atol is at least two
## decades below EXT_TOL. rtol stays looser because it binds on the large states
## (E_h ~ 100), where it costs the most and buys the least.
export ATOL=${ATOL:-1e-14}
export RTOL=${RTOL:-1e-8}
export EXT_TOL=${EXT_TOL:-1e-11}

## What makes two fixed points ONE attractor: Chebyshev on (Eh, PLP), both
## coordinates within this tolerance, one rule with no special cases. It was a gap
## in the composite burden, which is a projection and merged attractors 27 apart
## in P_LP. The tolerance is 1, not the composite's 2: Chebyshev is the stricter
## test at equal tolerance, and 1 unit is 1e11 organisms in P_LP, the scale at
## which two loads stop being the same outcome. pipeline/02_assemble.R asserts its own
## MERGE_TOL against the value recorded here.
export RAY_CLUST_TOL=${RAY_CLUST_TOL:-1.0}

## THE PROBE. Eh0 = 0, 10, ... 90, 99 -- eleven values.
## NOT 100: in the sterile system Eh0 = 100 gives Eh = K_E, Ea = En = 0, which IS
## an equilibrium (measured max|x(T) - x(0)| = 0 over TMAX = 12000). It is a
## saddle, leading eigenvalue +0.059 along CLP, so a final state resting on it is
## in no basin and must be labelled unassigned rather than high basin. 99 escapes
## it long before TMAX. Eh0 = 0 is kept: it moves, and it is the damaged end.
## Spacing 10, not v5's 5: compared directly on 20 combinations / 36 arms across
## both branches the two spacings gave the same crossing count every time, and
## the threshold is now bisected rather than read off the spacing.
export EH0_LO=${EH0_LO:-0}
export EH0_HI=${EH0_HI:-99}
export PROBE_BY=${PROBE_BY:-10}
## Bisection tolerance for Eh0_star. A PROBE_BY-wide bracket needs
## ceiling(log2(PROBE_BY / EH0_TOL)) halvings, so 7 at these values. 0.1 is well
## inside what the integrator resolves (dev on a correct run is ~1e-5).
export EH0_TOL=${EH0_TOL:-0.1}

export GRID_FP=${GRID_FP:-1}      # also run the grid-seeded finder -> saddles,
                                  # and n_field_stable: what EXISTS off the line
export PLP0_MAIN=${PLP0_MAIN:-10} # single PLP0 level
export OUT_PREFIX=${OUT_PREFIX:-$HOME/scratch/fp_v6_}
## Cell-list mode of 01_sweep.R (not used by the paper): CELLS restricts each tau row
## to its listed (Kg, r_L_P) cells; BETA_ROS_VALS replaces the beta_ros grid. Both
## empty = the main sweep.
export CELLS=${CELLS:-}
export BETA_ROS_VALS=${BETA_ROS_VALS:-}
export FLUSH_EVERY=${FLUSH_EVERY:-1}   # write output after every tau row

echo "sweep: batch $BATCH, row chunk $ROW_CHUNK_ID / $N_ROW_CHUNKS on $(hostname)"
echo "cells: ${CELLS:-(full grid)} | beta_ros: ${BETA_ROS_VALS:-(main grid)}"
echo "ROW_N=$ROW_N ROW_SEED=$ROW_SEED TMAX=$TMAX HMAX=$HMAX"
echo "atol=$ATOL rtol=$RTOL ext_tol=$EXT_TOL ray_clust_tol=$RAY_CLUST_TOL"
echo "probe Eh0=$EH0_LO..$EH0_HI by $PROBE_BY | bisect to $EH0_TOL | GRID_FP=$GRID_FP"
echo "output: ${OUT_PREFIX}${BATCH}"

Rscript pipeline/01_sweep.R \
    $N_ROW_CHUNKS \
    $ROW_CHUNK_ID \
    $BATCH

echo "Chunk completed"
