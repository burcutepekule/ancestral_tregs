#!/bin/bash
#SBATCH --job-name=cor_v6
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --array=0-399
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs_corridor/cor_%x_%A_%a.out
#SBATCH --error=logs_corridor/cor_%x_%A_%a.err

## ---------------------------------------------------------------------------
## pipeline/cluster/corridor_tzero_array.sh -- the t = 0 corridor stage, one branch per submission.
##
##   BRANCH=sterile sbatch pipeline/cluster/corridor_tzero_array.sh
##
## WHAT ONE TASK DOES
## Takes a disjoint slice of the gated combinations in rescue_<BRANCH>.rds and,
## for each, finds a base state on the Eh0 line and then walks outward along every
## swept coordinate under both Treg conditions, bisecting the OFF wall and the ON
## wall wherever the labels change. Output is long: one row per combination x
## coordinate.
##
## COST, MEASURED ON DELLA (2026-09-14, array cor_v6, progress logs):
##   bistable  11-17 min per combination, mono_B 25-38 min; mean ~15 min.
## The pool is EVERY bistable and mono_B row of rescue_<BRANCH>.rds:
##   sterile          89,328  ->  224 per task at 400 tasks  ->  ~57 h per task
##   pathogenic_help  17,495  ->   44 per task at 400 tasks  ->  ~11 h per task
## so the 12 h / 400-task default below cannot finish sterile. The first array
## timed out at ~47 of 224 per task. What it did finish is recovered by
## pipeline/cluster/corridor_merge_partials.sh, and pipeline/cluster/corridor_resubmit_round2.sh sizes the array
## for the remainder from the measured cost. An earlier version of this header
## budgeted 8 core-min and a 9,559-combination pool; both were wrong.
##
## THIS STAGE IS ~15x THE PER-COMBINATION COST OF THE SWEEP. The sweep integrates
## the Eh0 line, 36 integrations per combination. This walks 7 coordinates x 2
## conditions x (an 80-point scan, doubled up to 10 times, plus 30 bisection
## halvings per wall). Do not run it on anything the gate has not already reduced.
##
## CHUNKING. Same contract as the sweep: every task rebuilds the SAME pool
## (deterministic given BRANCH and the rescue file) and takes a disjoint slice, so
## the union of chunks equals one serial run.
##
## CRITICAL: the --array range and N_CHUNKS must match.
##     --array=0-(K-1)   <=>   N_CHUNKS=K
## ---------------------------------------------------------------------------
set -u

## Conda's activation scripts read $PS1, which is unset in a non-interactive
## shell. Under `set -u` that is a fatal "PS1: unbound variable" before anything
## runs. Give PS1 a value and drop -u across the conda block only.
export PS1=${PS1:-}
set +u
module load anaconda3/2023.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_Treg
set -u

SRC=${SRC:-$HOME/tregs_ode_paper}   # the repository root on the cluster
cd "$SRC" || exit 1
mkdir -p logs_corridor

export BRANCH=${BRANCH:-sterile}
export SRC
N_CHUNKS=${N_CHUNKS:-400}                  # = array size above; keep in sync!
CHUNK_ID=$(( SLURM_ARRAY_TASK_ID + 1 ))

export N_CORES=1              # one core per task; the parallelism IS the array
export N_SAMP=0               # 0 = every gated combination. See the note below.
export N_SAMP_MODE=0          # per-mode sampling is for pilots only

## Integration settings. These MUST match pipeline/01_sweep.R, because the labels here
## and the thresholds there have to describe the same dynamics -- a corridor
## measured under a different hmax is not the boundary the sweep located.
export TMAX=${TMAX:-12000}
export HMAX=${HMAX:-5}        # NOT optional -- see pipeline/cluster/submit_sweep.sh
## MUST MATCH pipeline/cluster/submit_sweep.sh. The clustering rule is shared through
## lib/attractor_cluster.R and the tolerance is the sweep's.
export RAY_CLUST_TOL=${RAY_CLUST_TOL:-1.0}

## THE SOLVER SETTINGS, STATED RATHER THAN INHERITED. pipeline/21_corridor_tzero.R gets
## them from lib/extinction_floor.R, whose defaults are these values, so
## nothing changes by exporting them -- except that a future divergence from
## pipeline/cluster/submit_sweep.sh becomes visible in the two files side by side instead of hiding
## in a default.
export ATOL=${ATOL:-1e-14}
export RTOL=${RTOL:-1e-8}
export EXT_TOL=${EXT_TOL:-1e-11}
## The Eh0 probe, used here ONLY to seed the fixed-point finder from the final
## states. v5's version of this stage called the finder with the grid alone, which
## is the same call that made it report monostable on combinations whose Eh0 line
## plainly reached two attractors.
export EH0_LO=${EH0_LO:-0}
export EH0_HI=${EH0_HI:-99}
export PROBE_BY=${PROBE_BY:-10}

## Scan resolution. NSCAN points per outward pass, the range doubling up to
## EXPAND_MAX times, then BIS_ITERS halvings per wall. Cost is LINEAR in NSCAN and
## in EXPAND_MAX, only logarithmic in the final precision -- so if this needs to be
## made cheaper, cut NSCAN before cutting BIS_ITERS.
export NSCAN=${NSCAN:-80}
export BIS_ITERS=${BIS_ITERS:-30}
export EXPAND_MAX=${EXPAND_MAX:-10}

## SWEEP_ONLY restricts which coordinates are measured. A fix confined to one
## coordinate's scan does not justify recomputing the others, which would come
## back identical. RUN_TAG then keeps the partial run's chunks out of the full
## run's filenames, so the merge cannot mix them:
##   BRANCH=sterile SWEEP_ONLY=M1 RUN_TAG=_m1 sbatch pipeline/cluster/corridor_tzero_array.sh
export SWEEP_ONLY=${SWEEP_ONLY:-}
## EXCLUDE: an .rds whose combinations are dropped from the pool before chunking
## (the merged output of a timed-out array). Set by pipeline/cluster/corridor_resubmit_round2.sh.
export EXCLUDE=${EXCLUDE:-}
## Write the .partial after every combination. See pipeline/21_corridor_tzero.R.
export FLUSH_EVERY=${FLUSH_EVERY:-1}
## Where the base state sits on the resolving side of the OFF threshold. The
## default puts it 2% past. Vary it to measure how much the tolerance numbers
## depend on that convention.
export BASE_FRAC=${BASE_FRAC:-0.02,0.05,0.10,0.20,0.35,0.50,0.75,1.00}

export OUT_TAG="${RUN_TAG:-}_chunk$(printf '%04d' "$CHUNK_ID")"

echo "corridor v6: branch $BRANCH, chunk $CHUNK_ID / $N_CHUNKS on $(hostname)"
echo "coords=${SWEEP_ONLY:-all} run_tag=${RUN_TAG:-none} exclude=${EXCLUDE:-none} flush_every=$FLUSH_EVERY base_frac=$BASE_FRAC"
echo "TMAX=$TMAX HMAX=$HMAX | probe $EH0_LO..$EH0_HI by $PROBE_BY | NSCAN=$NSCAN BIS_ITERS=$BIS_ITERS"

Rscript pipeline/21_corridor_tzero.R "$N_CHUNKS" "$CHUNK_ID"

echo "Chunk completed"
