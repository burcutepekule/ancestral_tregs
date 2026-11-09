#!/bin/bash
#SBATCH --job-name=imperfect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-399
#SBATCH --time=23:59:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs_imperfect/imperfect_%A_%a.out
#SBATCH --error=logs_imperfect/imperfect_%A_%a.err

## ---------------------------------------------------------------------------
## pipeline/cluster/imperfect_sweep.sh -- the imperfect-Treg sweep behind
## Fig. 5 and Extended Data Figs. 6 and 7, as a Slurm array (della).
##
## ONE SUBMISSION PER (BRANCH, WAVE); both are required. BRANCH sets the
## pathogen seeding (lib/skeleton_defs.R) and must match the units file;
## WAVE selects one of the waves stamped by pipeline/60_imperfect_sample.R.
## With the stored draw (977 units, R_REP = 10): waves 1, 2 = sterile (203, 202
## units), waves 3, 4 = pathogenic_help (286, 286 units).
##
##   Rscript pipeline/60_imperfect_sample.R      # once, locally
##   # copy data/fig05_units_*.rds to the cluster clone
##   sbatch --export=ALL,BRANCH=sterile,WAVE=1         pipeline/cluster/imperfect_sweep.sh
##   sbatch --export=ALL,BRANCH=sterile,WAVE=2         pipeline/cluster/imperfect_sweep.sh
##   sbatch --export=ALL,BRANCH=pathogenic_help,WAVE=3 pipeline/cluster/imperfect_sweep.sh
##   sbatch --export=ALL,BRANCH=pathogenic_help,WAVE=4 pipeline/cluster/imperfect_sweep.sh
##   # one wave at a time (each uses all 400 tasks). A task killed at the
##   # walltime is resubmitted with --array=<ids>; it resumes (RESUME=1).
##   # then: Rscript pipeline/62_imperfect_read.R  (on the cluster, reads the raw csv)
##
## The task slices (unit, replicate) items; one item is the full 11 x 11 grid
## for one unit and one noise realisation. Cost is dominated by the eta noise
## (the solver restarts at every unit-time jump of the perceived signal).
## Output goes to $SCRATCH_BASE/fig05_raw/<branch>_w<wave>/fig05_raw_t<task>.csv.
## ---------------------------------------------------------------------------
set -eu

export PS1=${PS1:-}
set +u
module load anaconda3/2023.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_Treg
set -u

export SRC=${SRC:-$HOME/tregs_ode_paper}
cd "$SRC" || exit 1

if [ -z "${BRANCH:-}" ]; then
  echo "BRANCH is not set. Submit as: sbatch --export=ALL,BRANCH=sterile,WAVE=1 pipeline/cluster/imperfect_sweep.sh" >&2
  exit 1
fi
if [ -z "${WAVE:-}" ]; then
  echo "WAVE is not set. Each submission runs ONE wave of the sampler's draw." >&2
  exit 1
fi
export WAVE
UNITS=${SRC}/data/fig05_units_${BRANCH}.rds
if [ ! -f "$UNITS" ]; then
  echo "missing $UNITS -- run pipeline/60_imperfect_sample.R locally and copy it here" >&2
  exit 1
fi

export SCRATCH_BASE=${SCRATCH_BASE:-$HOME/scratch}
if [ ! -d "$SCRATCH_BASE" ]; then
  echo "SCRATCH_BASE=$SCRATCH_BASE does not exist. Set it to your group scratch." >&2
  exit 1
fi
mkdir -p logs_imperfect "${SCRATCH_BASE}/fig05_raw/${BRANCH}_w${WAVE}"

export ARRAY_N=${SLURM_ARRAY_TASK_COUNT:-400}
export ARRAY_ID=${SLURM_ARRAY_TASK_ID:-0}
export IN_UNITS="$UNITS"
export OUTDIR="${SCRATCH_BASE}/fig05_raw/${BRANCH}_w${WAVE}"

## grid, replicates, horizon (as used for the stored data)
export PHI_N=${PHI_N:-11}
export ETA_N=${ETA_N:-11}
export R_REP=${R_REP:-10}
export TMAX=${TMAX:-12000}
export TSTEP=${TSTEP:-1}
export TAIL_FRAC=${TAIL_FRAC:-0.5}

## solver settings of the fixed-point sweep; hmax is overridden per call to TSTEP
export ATOL=${ATOL:-1e-14}
export RTOL=${RTOL:-1e-8}
export EXT_TOL=${EXT_TOL:-1e-11}
export HMAX=${HMAX:-5}

echo "imperfect ${BRANCH} wave ${WAVE} task ${ARRAY_ID}/${ARRAY_N} on $(hostname): ${PHI_N}x${ETA_N}, R=${R_REP}, tmax=${TMAX}"
Rscript "${SRC}/pipeline/61_imperfect_sweep.R"
echo "done"
