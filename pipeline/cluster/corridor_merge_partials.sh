#!/bin/bash
## ---------------------------------------------------------------------------
## pipeline/cluster/corridor_merge_partials.sh -- recover what a timed-out corridor array
## finished, for one branch.
##
##   BRANCH=sterile bash pipeline/cluster/corridor_merge_partials.sh
##
## Reads every corridor_tzero_<BRANCH>_chunk*.rds (tasks that finished) and every
## corridor_tzero_<BRANCH>_chunk*.rds.partial (tasks killed at the wall; the
## .partial holds only completed combinations, flushed by pipeline/21_corridor_tzero.R).
## Where both exist for a chunk the finished file wins.
##
## Writes corridor_tzero_<BRANCH>_round1.rds. The pool is permuted with a fixed
## seed before chunking, so the prefix each task completed is a uniform random
## subsample of the pool -- usable on its own, and the EXCLUDE input for
## pipeline/cluster/corridor_resubmit_round2.sh, which measures the rest as round 2.
##
## Refuses to write if two files cover the same (combination, coord), or if any
## measured combination is not in the pool.
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
export BRANCH=${BRANCH:-sterile}
export RUN_TAG=${RUN_TAG:-}

NFULL=$(ls -1 "data/corridor_tzero_${BRANCH}${RUN_TAG}_chunk"*.rds 2>/dev/null | wc -l | tr -d ' ')
NPART=$(ls -1 "data/corridor_tzero_${BRANCH}${RUN_TAG}_chunk"*.rds.partial 2>/dev/null | wc -l | tr -d ' ')
echo "== round 1 | branch $BRANCH | $NFULL finished chunks, $NPART partial chunks =="
if [ $((NFULL + NPART)) -eq 0 ]; then echo "nothing to merge" >&2; exit 1; fi

Rscript -e '
suppressMessages(library(data.table))
br  = Sys.getenv("BRANCH"); tag = Sys.getenv("RUN_TAG", "")
CK  = c("tau_id","Kg","beta_ros","r_L_P")
full = sort(Sys.glob(sprintf("./data/corridor_tzero_%s%s_chunk*.rds", br, tag)))
part = sort(Sys.glob(sprintf("./data/corridor_tzero_%s%s_chunk*.rds.partial", br, tag)))
part = part[!sub("[.]partial$", "", part) %in% full]     # finished file wins
rd = function(f) { x = as.data.table(readRDS(f)); x[, src := basename(f)]; x }
L = rbindlist(lapply(c(full, part), rd), fill = TRUE)
if (uniqueN(L, by = c(CK, "coord")) < nrow(L))
  stop("duplicate (combination, coord) rows across files -- two files cover the same slice")
per = L[, .(combos = uniqueN(.SD)), by = src, .SDcols = CK]
cat(sprintf("[round1] %d finished + %d partial files -> %s rows, %s combinations\n",
            length(full), length(part), format(nrow(L), big.mark=","),
            format(uniqueN(L[, ..CK]), big.mark=",")))
cat(sprintf("[round1] combinations per file: min %d, median %g, max %d\n",
            min(per$combos), median(per$combos), max(per$combos)))
R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", br)))
U = R[rescue_mode %in% c("bistable", "mono_B")]
K = unique(L[, ..CK])
out = K[!U, on = CK]
if (nrow(out)) stop(nrow(out), " measured combinations are not in the pool -- wrong rescue file?")
cat(sprintf("[round1] %s of %s pool combinations measured (%.1f%%); %s remain for round 2\n",
            format(nrow(K), big.mark=","), format(nrow(U), big.mark=","),
            100 * nrow(K) / nrow(U), format(nrow(U) - nrow(K), big.mark=",")))
L[, src := NULL]
OUT = sprintf("./data/corridor_tzero_%s%s_round1.rds", br, tag)
saveRDS(L, OUT)
cat(sprintf("[round1] wrote %s\n", OUT))
print(L[, .N, by = note][order(-N)])
'
rc=$?
OUT="data/corridor_tzero_${BRANCH}${RUN_TAG}_round1.rds"
if [ $rc -ne 0 ] || [ ! -f "$OUT" ]; then echo "== ROUND-1 MERGE FAILED (exit $rc) ==" >&2; exit 1; fi
echo "== round 1 done | $(ls -lh "$OUT" | awk '{print $5}') =="
