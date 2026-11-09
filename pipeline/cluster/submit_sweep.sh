#!/bin/bash
## ---------------------------------------------------------------------------
## submit_sweep.sh -- submit pipeline/cluster/sweep_array.sh once per tau batch.
##
##   bash pipeline/cluster/submit_sweep.sh                 # every batch with a grid
##   BATCHES="1 2" bash pipeline/cluster/submit_sweep.sh   # just those
##   DRY=1 bash pipeline/cluster/submit_sweep.sh           # print the sbatch lines only
##
## The batch list is DERIVED from the tau grids on disk, not hard-coded.
##
## RUN THE TIMING CHECK FIRST -- once, not per batch:
##   BATCH=1 TIMEIT=1 TIME_N=5 Rscript pipeline/01_sweep.R 400 1 1
## ---------------------------------------------------------------------------
set -u
export SRC=${SRC:-$HOME/tregs_ode_paper}
cd "$SRC" || exit 1

if [ -n "${BATCHES:-}" ]; then
  list="$BATCHES"
else
  list=$(ls -1 data/tau_grid_v6_*.rds 2>/dev/null \
         | sed 's/.*tau_grid_v6_\([0-9]*\)\.rds/\1/' | sort -n)
fi
[ -n "$list" ] || { echo "no batches found -- run pipeline/00_parameter_grids.R first" >&2; exit 1; }

echo "batches to submit: $(echo $list | tr '\n' ' ')"
for b in $list; do
  f="data/tau_grid_v6_${b}.rds"
  [ -f "$f" ] || { echo "  SKIP batch $b: $f missing" >&2; continue; }
  cmd="sbatch --job-name=fp_sweep_${b} --export=ALL,BATCH=${b} pipeline/cluster/sweep_array.sh"
  echo "  $cmd"
  [ "${DRY:-0}" = "1" ] || eval "$cmd"
done
