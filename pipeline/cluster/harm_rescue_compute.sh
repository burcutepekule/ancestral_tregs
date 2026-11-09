#!/bin/bash
#SBATCH --job-name=harm_rescue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --output=logs_harm_rescue/%x_%j.out
#SBATCH --error=logs_harm_rescue/%x_%j.err

## ---------------------------------------------------------------------------
## pipeline/cluster/harm_rescue_compute.sh -- pipeline/50_harm_rescue_compute.R on della.
##
## WHY IT IS ONE TASK WITH MANY CORES, not an array. 50_harm_rescue_compute.R
## parallelises inside R with mclapply over combinations, and it writes one
## file at the end that the three plotting scripts read. Splitting it across
## array tasks would mean merging four partial files afterwards and carrying a
## settings stamp on each, which is exactly the arrangement the restructure
## removed. One task, one output, no reconciliation.
##
## ASK FOR MORE CORES IF YOU WANT IT FASTER. The work is one integration per
## (combination, arm) and one per bisection step, all independent, so it scales
## with cores until the node runs out of them. N_CORES follows
## SLURM_CPUS_PER_TASK, so raising the allocation is the only change needed:
##
##      sbatch --cpus-per-task=64 pipeline/cluster/harm_rescue_compute.sh
##
## Two things bound it. mclapply forks within ONE node, so the ceiling is that
## node's core count -- going wider needs an array over chunks of combinations
## and a merge afterwards, which is the sweep's pattern and not built here. And
## a larger allocation queues longer, so on a busy cluster 32 cores starting now
## can finish before 96 cores that start in four hours.
##
## MEMORY IS A TOTAL, NOT PER CORE, AND 48G IS GENEROUS. Measured on a running
## job: the parent process sits at 84 MB and each forked worker at 86 MB, so
## the workers are sharing almost everything copy-on-write rather than each
## holding a copy. What grows with the problem is the trajectory table -- one
## row per combination per arm per output time, which on the full gated set is
## about 18 million rows and 1.3 GB, transiently double that while the
## per-combination pieces are being bound together. Adding a few hundred MB per
## worker on top, a 64-core run should peak well under 20 GB.
##
## --mem-per-cpu would have charged for this per core, so a 64-core request
## would have asked for 256 GB it never touches and queued behind jobs that
## needed it. Check the real figure with `seff <jobid>` after the first run and
## trim this line to it.
##
## RUN: sbatch pipeline/cluster/harm_rescue_compute.sh
##      sbatch --cpus-per-task=64 pipeline/cluster/harm_rescue_compute.sh              # faster if free
##      THR_TREATS=delay sbatch pipeline/cluster/harm_rescue_compute.sh     # delay thresholds only, quick
##      SAMPLE=200 sbatch pipeline/cluster/harm_rescue_compute.sh           # a test on 200 combinations
## ---------------------------------------------------------------------------
set -eu

## Conda's activation scripts read $PS1, which is unset in a non-interactive
## shell. Under `set -u` that is a fatal "PS1: unbound variable" before anything
## runs, so PS1 gets a value and -u is dropped across the conda block only.
export PS1=${PS1:-}
set +u
module load anaconda3/2023.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_Treg
set -u

## the R script uses relative paths (./lib/..., ./data/...) under SRC
cd ${SRC:-$HOME/tregs_ode_paper} || exit 1
mkdir -p logs_harm_rescue

## N_CORES must match the allocation or the cores are either idle or oversubscribed
export N_CORES=${N_CORES:-${SLURM_CPUS_PER_TASK:-32}}
export SAMPLE=${SAMPLE:-0}                  # 0 = every gated combination
export SRC=${SRC:-$HOME/tregs_ode_paper}

## THESE MUST MATCH THE SWEEP. They are the same variables FP_run_v6.sh exports,
## and a figure integrated under different tolerances from the sweep that
## labelled its combinations is not a figure about those labels.
export ATOL=${ATOL:-1e-14}
export RTOL=${RTOL:-1e-8}
export EXT_TOL=${EXT_TOL:-1e-11}
export HMAX=${HMAX:-5}
export TMAX=${TMAX:-25000}

## THE TREATMENT LADDERS, STATED HERE TOO. They are the defaults in
## 50_harm_rescue_compute.R, repeated so a cluster run and a laptop run can be compared
## line by line rather than by trusting that neither was overridden. The two
## descending ladders run to 1/133 so all three cover the same fold range as
## tau_t_d's 128.
export DELAYS=${DELAYS:-0,10,20,40,80,160,320,960,2400}
export SCALES=${SCALES:-1,0.75,0.5,0.35,0.25,0.15,0.1,0.05,0.03,0.02,0.0125,0.0075}
export DEACTS=${DEACTS:-1,2,4,8,16,32,64,128}
export INDUCE=${INDUCE:-1,0.75,0.5,0.35,0.25,0.15,0.1,0.05,0.03,0.02,0.0125,0.0075}

## the threshold pass
export THR_TREATS=${THR_TREATS:-delay,strength,deact,induce}
export THR_BOUNDS=${THR_BOUNDS:-neutral,help}
export DTOL=${DTOL:-1}                      # absolute, time units, for the delay
export LTOL=${LTOL:-0.01}                   # relative, for the three multipliers
export EFFECT_TOL=${EFFECT_TOL:-1.0}

echo "harm rescue compute on $(hostname): ${N_CORES} cores, SAMPLE=${SAMPLE}"
echo "thresholds: ${THR_TREATS} x ${THR_BOUNDS}"
Rscript ${SRC}/pipeline/50_harm_rescue_compute.R
echo "done"
