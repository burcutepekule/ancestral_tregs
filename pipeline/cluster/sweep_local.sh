#!/bin/bash
## ---------------------------------------------------------------------------
## sweep_local.sh -- pipeline/01_sweep.R on this machine, N_CORES processes at once.
##
## Same chunking contract as the cluster array: every process builds the SAME
## tau-row subsample (deterministic given ROW_SEED and ROW_N) and takes a disjoint
## slice, so the union is identical to one serial run and the chunks merge exactly
## as the cluster ones do.
##
## SIZING. Measured head to head against v5 on identical hardware: 19.3 s per
## combination against v5's 21.5, and 315 combinations per tau row against 280 --
## so about 1.7 h per tau row on one core, the same as v5. Then
##       wall clock ~= ROW_N / N_CORES * 1.7 h
## Defaults (ROW_N=60, N_CORES=10) => 6 rows per core => ~10 h. Each row is
## flushed as it completes, so killing it early just gives you fewer rows.
## MEASURE IT HERE FIRST rather than trusting that: TIMEIT=1 TIME_N=5 bash this.
##
## Writes to ./data/sweep_raw/fp_local_batch_<BATCH>/ -- a DIFFERENT directory from
## the cluster output, so the two never mix. They cover different tau-row subsets
## and mixing them would make "how many rows" meaningless.
##
## RUN (from the repository root):
##       bash pipeline/cluster/sweep_local.sh
##       ROW_N=120 N_CORES=6 bash pipeline/cluster/sweep_local.sh
##       caffeinate -i bash pipeline/cluster/sweep_local.sh      # macOS: stop it sleeping
## WATCH: tail -f logs_sweep_local/chunk01.log
## THEN:  DATA_DIR=./data/sweep_raw/fp_local_batch_1 BATCH=1 Rscript pipeline/02_assemble.R
## ---------------------------------------------------------------------------
set -u
N_CORES=${N_CORES:-6}
ROW_N=${ROW_N:-60}
BATCH=${BATCH:-1}

## Every constant below carries its justification in sweep_array.sh; the two files
## MUST agree, because a local run and a cluster run are merged into one dataset.
export TMAX=${TMAX:-12000}
export MAXSTEPS=${MAXSTEPS:-200000}
export HMAX=${HMAX:-5}            # NOT optional -- see sweep_array.sh
export EH0_LO=${EH0_LO:-0}        # probe 0, 10, ... 90, 99
export EH0_HI=${EH0_HI:-99}       # NOT 100: that is the sterile L=0 saddle
export PROBE_BY=${PROBE_BY:-10}
export EH0_TOL=${EH0_TOL:-0.1}    # bisection tolerance for Eh0_star
export PLP0_MAIN=${PLP0_MAIN:-10}
export GRID_FP=${GRID_FP:-1}
export FLUSH_EVERY=${FLUSH_EVERY:-1}
export ROW_N
export ROW_SEED=${ROW_SEED:-20260808}     # identical in every process; see header
export OUT_PREFIX=./data/sweep_raw/fp_local_batch_
export ATOL=${ATOL:-1e-14}
export RTOL=${RTOL:-1e-8}
export EXT_TOL=${EXT_TOL:-1e-11}
export RAY_CLUST_TOL=${RAY_CLUST_TOL:-1.0}
unset CELLS BETA_ROS_VALS        # the main grid

cd "$(dirname "$0")/../.." || exit 1
export SRC="$PWD"
if [ ! -f "./data/tau_grid_v6_${BATCH}.rds" ]; then
  echo "no ./data/tau_grid_v6_${BATCH}.rds -- run pipeline/00_parameter_grids.R first" >&2; exit 1
fi

mkdir -p logs_sweep_local "${OUT_PREFIX}${BATCH}"
echo "local sweep | cores $N_CORES | tau rows $ROW_N | TMAX $TMAX | hmax $HMAX"
echo "probe Eh0 $EH0_LO..$EH0_HI by $PROBE_BY | bisect to $EH0_TOL | batch $BATCH"
echo "estimated wall clock: $(python3 -c "print(f'{$ROW_N/$N_CORES*1.7:.1f}')") h"
echo "output: ${OUT_PREFIX}${BATCH}/   logs: logs_sweep_local/"

## TIMEIT short-circuits inside the R script and quits, so run ONE process only.
if [ "${TIMEIT:-0}" = "1" ]; then
  export TIMEIT TIME_N=${TIME_N:-5}
  Rscript pipeline/01_sweep.R "$N_CORES" 1 "$BATCH"
  exit $?
fi

for i in $(seq 1 "$N_CORES"); do
  log=$(printf "logs_sweep_local/chunk%02d.log" "$i")
  nohup Rscript pipeline/01_sweep.R "$N_CORES" "$i" "$BATCH" > "$log" 2>&1 &
  echo "  started chunk $i/$N_CORES -> $log (pid $!)"
done
wait
