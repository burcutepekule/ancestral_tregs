#!/bin/bash
## ---------------------------------------------------------------------------
## pipeline/cluster/corridor_combine_rounds.sh -- round 1 + round 2 -> corridor_tzero_<BRANCH>.rds
##
##   BRANCH=sterile bash pipeline/cluster/corridor_combine_rounds.sh
##
## Reads corridor_tzero_<BRANCH>_round1.rds and the round-2 chunk files named
## in corridor_round2_plan_<BRANCH>.txt. Refuses if any round-2 chunk is
## missing, if any (combination, coord) is covered twice, or if the union does
## not equal the pool exactly. This replaces pipeline/cluster/corridor_merge_chunks.sh for a run
## that was done in two rounds; the downstream contract is unchanged.
##
## To see which round-2 chunks are missing without merging:
##   BRANCH=sterile CHECK=1 bash pipeline/cluster/corridor_combine_rounds.sh
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
PLAN="data/corridor_round2_plan_${BRANCH}.txt"
[ -f "$PLAN" ] || { echo "no $PLAN -- was round 2 submitted with pipeline/cluster/corridor_resubmit_round2.sh?" >&2; exit 1; }
RUN_TAG=$(awk -F= '/^RUN_TAG=/{print $2}' "$PLAN")
N_CHUNKS=$(awk -F= '/^N_CHUNKS=/{print $2}' "$PLAN")
export BRANCH RUN_TAG N_CHUNKS
R1="data/corridor_tzero_${BRANCH}_round1.rds"
[ -f "$R1" ] || { echo "round-1 file $R1 not found" >&2; exit 1; }

MISSING=$(ls -1 "data/corridor_tzero_${BRANCH}${RUN_TAG}_chunk"*.rds 2>/dev/null | sed 's/.*chunk//; s/\.rds$//' \
  | awk -v n="$N_CHUNKS" '{h[$1+0]=1} END{for(i=1;i<=n;i++) if(!(i in h)) print i}')
NMISS=$(echo -n "$MISSING" | grep -c .)
echo "== combine | branch $BRANCH | round-2 tag $RUN_TAG | $((N_CHUNKS - NMISS)) of $N_CHUNKS chunks present =="
if [ "$NMISS" -gt 0 ]; then
  echo "missing round-2 chunks (array index = chunk - 1):" >&2
  echo "$MISSING" | tr '\n' ' ' >&2; echo >&2
  echo "resubmit one with: BRANCH=$BRANCH N_CHUNKS=$N_CHUNKS RUN_TAG=$RUN_TAG EXCLUDE=./data/corridor_tzero_${BRANCH}_round1.rds sbatch --array=<chunk-1> --time=48:00:00 pipeline/cluster/corridor_tzero_array.sh" >&2
  exit 1
fi
[ "${CHECK:-0}" = "1" ] && { echo "all round-2 chunks present"; exit 0; }

Rscript -e '
suppressMessages(library(data.table))
br = Sys.getenv("BRANCH"); tag = Sys.getenv("RUN_TAG")
CK = c("tau_id","Kg","beta_ros","r_L_P")
f  = sort(Sys.glob(sprintf("./data/corridor_tzero_%s%s_chunk*.rds", br, tag)))
L2 = rbindlist(lapply(f, readRDS), fill = TRUE)
L1 = as.data.table(readRDS(sprintf("./data/corridor_tzero_%s_round1.rds", br)))
L  = rbindlist(list(L1, L2), fill = TRUE)
if (uniqueN(L, by = c(CK, "coord")) < nrow(L))
  stop("duplicate (combination, coord) rows -- round 2 re-measured something in round 1")
R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", br)))
U = unique(R[rescue_mode %in% c("bistable", "mono_B"), ..CK])
K = unique(L[, ..CK])
miss = U[!K, on = CK]; extra = K[!U, on = CK]
cat(sprintf("[combine] round 1: %s combinations | round 2: %s combinations in %d chunks\n",
            format(uniqueN(L1[, ..CK]), big.mark=","), format(uniqueN(L2[, ..CK]), big.mark=","), length(f)))
cat(sprintf("[combine] union %s of %s pool combinations; %d missing, %d not in pool\n",
            format(nrow(K), big.mark=","), format(nrow(U), big.mark=","), nrow(miss), nrow(extra)))
if (nrow(miss) || nrow(extra)) stop("union does not equal the pool -- nothing written")
OUT = sprintf("./data/corridor_tzero_%s.rds", br)
saveRDS(L, OUT)
cat(sprintf("[combine] wrote %s: %s rows\n", OUT, format(nrow(L), big.mark=",")))
print(L[, .N, by = note][order(-N)])
'
rc=$?
OUT="data/corridor_tzero_${BRANCH}.rds"
if [ $rc -ne 0 ] || [ ! -f "$OUT" ]; then echo "== COMBINE FAILED (exit $rc) ==" >&2; exit 1; fi
echo "== combine done | $(ls -lh "$OUT" | awk '{print $5}') =="
