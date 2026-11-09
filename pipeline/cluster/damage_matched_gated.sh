#!/bin/bash
#SBATCH --job-name=dm_gated
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-399
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs_dm_gated/dm_gated_%A_%a.out
#SBATCH --error=logs_dm_gated/dm_gated_%A_%a.err

## ---------------------------------------------------------------------------
## pipeline/cluster/damage_matched_gated.sh -- the damage-matched experiment
## (pipeline/40_damage_matched.R, MATCH=conv) on the units of
## pipeline/43_damage_matched_gated_units.R: drawn from the whole Tregs-OFF bistable
## set with no Treg-effect criterion, each at its own gated beta_r, among the units
## in which Tregs OFF does not eliminate the pathogen. A unit whose Tregs-OFF run
## eliminates after all is recorded (status off_eliminates) and not searched.
## A Slurm array on della, one unit slice per task: task t runs the units with
## (i - 1) %% 400 == t (75,000 units: 187-188 per task, ~340 integrations each; at
## ~110 s per unit a task takes ~5.7 h). Every unit is saved as it finishes, and the
## units are in random order, so a task stopped at the walltime leaves a random
## subset behind: resubmit it to finish, or merge with ALLOW_MISSING=1.
##
##   # the cluster clone needs (data/ is not in git):
##   #   data/fig04_gated_units.rds, data/tau_grid_v6_{1..5}.rds
##   sbatch pipeline/cluster/damage_matched_gated.sh
##   # a task killed at the walltime: sbatch --array=<ids> ...; it resumes (RESUME=1)
##   # then, on the cluster:
##   #   Rscript pipeline/cluster/damage_matched_merge.R
##   # and copy data/fig04_damage_matched_gated{,_profile}.rds back
## Output: $SCRATCH_BASE/fig04_gated/dm_gated_t<task>.rds (+ _profile.rds)
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
for f in data/fig04_gated_units.rds data/tau_grid_v6_1.rds; do
  if [ ! -f "$f" ]; then echo "missing $SRC/$f -- copy data/fig04_gated_units.rds and data/tau_grid_v6_*.rds here" >&2; exit 1; fi
done

export SCRATCH_BASE=${SCRATCH_BASE:-$HOME/scratch}
if [ ! -d "$SCRATCH_BASE" ]; then
  echo "SCRATCH_BASE=$SCRATCH_BASE does not exist. Set it to your group scratch." >&2
  exit 1
fi
OUTDIR="${SCRATCH_BASE}/fig04_gated"
mkdir -p logs_dm_gated "$OUTDIR"

## the partition is ALWAYS into 400 slices: SLURM_ARRAY_TASK_COUNT is the number of tasks
## of THIS submission, so a resubmitted subset (--array=5,17) would otherwise re-partition
export ARRAY_N=${ARRAY_N:-400}
export ARRAY_ID=${SLURM_ARRAY_TASK_ID:-0}
export OUT="${OUTDIR}/dm_gated_t${ARRAY_ID}.rds"

## the search
export MATCH=${MATCH:-conv}
export BETA_LO=${BETA_LO:-0.0005}
export BETA_HI=${BETA_HI:-10}
export GRID_STEP=${GRID_STEP:-1.03}
export NEAR=${NEAR:-0.02}
export N_PROBE=${N_PROBE:-6}
export TOL=${TOL:-1e-3}
export TOL_BIS=${TOL_BIS:-1e-4}
export MAX_BIS=${MAX_BIS:-30}
## the units (do not change the file once tasks have run: the units are indexed by its rows)
export UNIT_FILE=${UNIT_FILE:-./data/fig04_gated_units.rds}
export SKIP_OFF_ELIM=1 ELIM_TOL=${ELIM_TOL:-1e-3}
export TMAX=${TMAX:-12000}
## one core per task; save after every unit so a resubmitted task resumes
export N_CORES=1 FLUSH=1 RESUME=1

echo "dm_gated task ${ARRAY_ID}/${ARRAY_N} on $(hostname): beta in [${BETA_LO}, ${BETA_HI}], grid x${GRID_STEP}"
Rscript "${SRC}/pipeline/40_damage_matched.R"
echo "done"
