#!/bin/bash
## ---------------------------------------------------------------------------
## pipeline/cluster/corridor_drain_submit.sh -- feed the rest of a round-2 array into the
## queue as room appears under della's submit limit.
##
##   BRANCH=sterile START=1000 nohup bash pipeline/cluster/corridor_drain_submit.sh > drain_sterile.log 2>&1 &
##
## Reads corridor_round2_plan_<BRANCH>.txt for N_CHUNKS, RUN_TAG and EXCLUDE, so
## every slice it submits is part of the SAME plan as the slices already in.
## START is the first array index not yet submitted (the one sbatch refused).
## Progress is kept in corridor_round2_next_<BRANCH>.txt, so the loop can be
## killed and restarted without START.
##
## Every POLL seconds: count this user's queued tasks (squeue -r, one line per
## array task, which is how the QOS limit counts), and if there is room submit
## the next slice of at most that many. Exits when the plan is fully submitted.
## Run it under nohup on the login node; it does no work itself.
## ---------------------------------------------------------------------------
set -u
SRC=${SRC:-$HOME/tregs_ode_paper}   # the repository root on the cluster
cd "$SRC" || exit 1
BRANCH=${BRANCH:?set BRANCH}
LIMIT=${LIMIT:-1000}        # QOSMaxSubmitJobPerUser, in tasks
HEADROOM=${HEADROOM:-20}    # stay this far under it: other jobs, races with the scheduler
POLL=${POLL:-600}           # seconds between checks
PLAN="data/corridor_round2_plan_${BRANCH}.txt"
NEXT="data/corridor_round2_next_${BRANCH}.txt"
[ -f "$PLAN" ] || { echo "no $PLAN" >&2; exit 1; }
N_CHUNKS=$(awk -F= '/^N_CHUNKS=/{print $2}' "$PLAN")
RUN_TAG=$(awk -F= '/^RUN_TAG=/{print $2}' "$PLAN")
EXCLUDE=$(awk -F= '/^EXCLUDE=/{print $2}' "$PLAN")
WALL_H=$(awk -F= '/^WALL_H=/{print $2}' "$PLAN")
if [ -f "$NEXT" ]; then next=$(cat "$NEXT"); else next=${START:?set START (first array index not yet submitted)}; fi
echo "[$(date '+%F %T')] plan N_CHUNKS=$N_CHUNKS tag=$RUN_TAG; next index $next; limit $LIMIT, headroom $HEADROOM, poll ${POLL}s"

while [ "$next" -lt "$N_CHUNKS" ]; do
  queued=$(squeue -u "$USER" -h -r 2>/dev/null | wc -l | tr -d ' ')
  room=$(( LIMIT - HEADROOM - queued ))
  if [ "$room" -gt 0 ]; then
    hi=$(( next + room - 1 )); [ "$hi" -ge "$N_CHUNKS" ] && hi=$(( N_CHUNKS - 1 ))
    if out=$(env BRANCH="$BRANCH" N_CHUNKS="$N_CHUNKS" RUN_TAG="$RUN_TAG" EXCLUDE="$EXCLUDE" FLUSH_EVERY=1 \
             sbatch --array="${next}-${hi}" --time="${WALL_H}:00:00" pipeline/cluster/corridor_tzero_array.sh 2>&1); then
      echo "[$(date '+%F %T')] queued $queued -> submitted --array=${next}-${hi}: $out"
      next=$(( hi + 1 )); echo "$next" > "$NEXT"
      continue          # re-count at once; more room may exist
    else
      echo "[$(date '+%F %T')] queued $queued, tried --array=${next}-${hi}, refused: $out"
    fi
  else
    echo "[$(date '+%F %T')] queued $queued, no room"
  fi
  sleep "$POLL"
done
echo "[$(date '+%F %T')] all $N_CHUNKS tasks submitted; done"
