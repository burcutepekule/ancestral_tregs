#!/bin/bash
## ---------------------------------------------------------------------------
## pipeline/cluster/corridor_resubmit_round2.sh -- submit round 2 of the corridor stage: every
## pool combination NOT in corridor_tzero_<BRANCH>_round1.rds, with the array
## sized from the measured cost so that each task fits the wall clock.
##
##   BRANCH=sterile bash pipeline/cluster/corridor_resubmit_round2.sh          # submits
##   BRANCH=sterile DRY=1 bash pipeline/cluster/corridor_resubmit_round2.sh    # prints the plan only
##
## Knobs (env):
##   WALL_H         wall clock per task in hours            default 24
##   MIN_PER_COMBO  measured minutes per combination        default 15
##   MARGIN         safety factor on that                   default 2
##   ARRAY_MAX      tasks per sbatch call                   default 1000
##   EXCLUDE        round-1 file                            default corridor_tzero_<BRANCH>_round1.rds
##   RUN_TAG        tag on the round-2 chunk files          default _round2
##
## Tasks per sbatch call are capped at ARRAY_MAX because Slurm caps the array
## size and della counts array tasks against MaxSubmitJobsPerUser. Every call
## passes the same N_CHUNKS, so the slices are disjoint and their union is the
## remaining pool; a later call that Slurm rejects is simply re-run once earlier
## tasks have drained. The plan (N_CHUNKS above all) is written to
## corridor_round2_plan_<BRANCH>.txt, which pipeline/cluster/corridor_combine_rounds.sh reads.
## ---------------------------------------------------------------------------
set -u
export PS1=${PS1:-}
set +u
module load anaconda3/2023.3 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_Treg
set -u
command -v Rscript >/dev/null || { echo "Rscript not found after activating env_Treg" >&2; exit 1; }

SRC=${SRC:-$HOME/tregs_ode_paper}   # the repository root on the cluster
cd "$SRC" || exit 1
BRANCH=${BRANCH:?set BRANCH}
WALL_H=${WALL_H:-24}
MIN_PER_COMBO=${MIN_PER_COMBO:-15}
MARGIN=${MARGIN:-2}
ARRAY_MAX=${ARRAY_MAX:-1000}
RUN_TAG=${RUN_TAG:-_round2}
EXCLUDE=${EXCLUDE:-./data/corridor_tzero_${BRANCH}_round1.rds}
DRY=${DRY:-0}
PLAN="data/corridor_round2_plan_${BRANCH}.txt"

[ -f "$EXCLUDE" ] || { echo "round-1 file not found: $EXCLUDE -- run pipeline/cluster/corridor_merge_partials.sh first" >&2; exit 1; }
if ls "data/corridor_tzero_${BRANCH}${RUN_TAG}_chunk"*.rds* >/dev/null 2>&1; then
  echo "round-2 chunk files already exist for tag $RUN_TAG -- move them aside or choose another RUN_TAG" >&2; exit 1
fi

## Remaining pool, computed exactly as pipeline/21_corridor_tzero.R will compute it.
N_REMAIN=$(BRANCH=$BRANCH EXCLUDE=$EXCLUDE Rscript -e '
suppressMessages(library(data.table))
CK = c("tau_id","Kg","beta_ros","r_L_P")
R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", Sys.getenv("BRANCH"))))
U = R[rescue_mode %in% c("bistable", "mono_B")]
X = unique(as.data.table(readRDS(Sys.getenv("EXCLUDE")))[, ..CK])
cat(nrow(U[!X, on = CK]))') || exit 1
[ "$N_REMAIN" -gt 0 ] || { echo "nothing remains: round 1 covers the pool"; exit 0; }

PER_TASK=$(awk -v w="$WALL_H" -v m="$MIN_PER_COMBO" -v s="$MARGIN" 'BEGIN{p=int(w*60/(m*s)); print (p<1?1:p)}')
N_CHUNKS=$(( (N_REMAIN + PER_TASK - 1) / PER_TASK ))
EST_H=$(awk -v p="$PER_TASK" -v m="$MIN_PER_COMBO" 'BEGIN{printf "%.1f", p*m/60}')
CORE_H=$(awk -v n="$N_REMAIN" -v m="$MIN_PER_COMBO" 'BEGIN{printf "%.0f", n*m/60}')

echo "== round 2 plan | branch $BRANCH =="
echo "   remaining combinations : $N_REMAIN"
echo "   per task               : $PER_TASK  (~${EST_H} h at $MIN_PER_COMBO min each; wall ${WALL_H} h, margin ${MARGIN}x)"
echo "   tasks (N_CHUNKS)       : $N_CHUNKS  in sbatch calls of at most $ARRAY_MAX"
echo "   total                  : ~$CORE_H core-hours"
LIM=$(sacctmgr -nP show assoc user="$USER" format=maxsubmitjobs,maxjobs 2>/dev/null | head -1)
[ -n "$LIM" ] && echo "   your assoc limits (maxsubmit|maxjobs): $LIM"
MA=$(scontrol show config 2>/dev/null | awk -F= '/MaxArraySize/{gsub(/ /,"",$2); print $2}')
[ -n "$MA" ] && echo "   MaxArraySize: $MA"

{ echo "BRANCH=$BRANCH"; echo "RUN_TAG=$RUN_TAG"; echo "EXCLUDE=$EXCLUDE"; echo "N_CHUNKS=$N_CHUNKS"
  echo "N_REMAIN=$N_REMAIN"; echo "PER_TASK=$PER_TASK"; echo "WALL_H=$WALL_H"; echo "planned=$(date '+%Y-%m-%dT%H:%M:%S')"; } > "$PLAN"
echo "   plan written to $PLAN"

lo=0
while [ "$lo" -lt "$N_CHUNKS" ]; do
  hi=$(( lo + ARRAY_MAX - 1 )); [ "$hi" -ge "$N_CHUNKS" ] && hi=$(( N_CHUNKS - 1 ))
  CMD="BRANCH=$BRANCH N_CHUNKS=$N_CHUNKS RUN_TAG=$RUN_TAG EXCLUDE=$EXCLUDE FLUSH_EVERY=1 sbatch --array=${lo}-${hi} --time=${WALL_H}:00:00 pipeline/cluster/corridor_tzero_array.sh"
  echo "   $CMD"
  if [ "$DRY" != "1" ]; then
    if ! env BRANCH="$BRANCH" N_CHUNKS="$N_CHUNKS" RUN_TAG="$RUN_TAG" EXCLUDE="$EXCLUDE" FLUSH_EVERY=1 \
         sbatch --array="${lo}-${hi}" --time="${WALL_H}:00:00" pipeline/cluster/corridor_tzero_array.sh; then
      echo "sbatch refused --array=${lo}-${hi}. Submit this and any later slice by hand once tasks drain," >&2
      echo "keeping N_CHUNKS=$N_CHUNKS exactly as above." >&2
      exit 1
    fi
  fi
  lo=$(( hi + 1 ))
done
[ "$DRY" = "1" ] && echo "   (DRY=1: nothing submitted)"
exit 0
