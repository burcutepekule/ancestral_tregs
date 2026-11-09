#!/bin/bash
## ---------------------------------------------------------------------------
## pipeline/cluster/corridor_merge_chunks.sh -- concatenate the corridor chunks for one branch.
##
##   BRANCH=sterile bash pipeline/cluster/corridor_merge_chunks.sh
##
## The chunk count is DERIVED from the files present and CHECKED against N_CHUNKS
## before anything is written. A merge over a partial set is the failure this
## pipeline has already paid for: it produces a file that looks complete and is
## quietly missing combinations, and nothing downstream can tell.
## ---------------------------------------------------------------------------
set -u

## Conda, because Rscript is not on the default path on a login node and this
## script previously died with "Rscript: command not found" and then printed
## "merge done" anyway. $PS1 is unset in a non-interactive shell and conda's
## activation reads it, which is fatal under `set -u` -- so give it a value and
## drop -u across the conda block only.
export PS1=${PS1:-}
set +u
module load anaconda3/2023.3 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_Treg
set -u
command -v Rscript >/dev/null || { echo "Rscript still not found after activating env_Treg" >&2; exit 1; }

SRC=${SRC:-$HOME/tregs_ode_paper}   # the repository root on the cluster
cd "$SRC" || exit 1
BRANCH=${BRANCH:-sterile}
N_CHUNKS=${N_CHUNKS:-400}
## RUN_TAG separates a partial re-run (e.g. SWEEP_ONLY=M1) from the full run, so
## the two sets of chunks cannot be merged into one file by accident. Exported
## because the R block below reads it.
export RUN_TAG=${RUN_TAG:-}

NCH=$(ls -1 "data/corridor_tzero_${BRANCH}${RUN_TAG:-}_chunk"*.rds 2>/dev/null | wc -l | tr -d ' ')
echo "== merge | branch $BRANCH | $NCH of $N_CHUNKS chunks present =="
if [ "$NCH" -ne "$N_CHUNKS" ]; then
  echo "REFUSING: $((N_CHUNKS - NCH)) chunk(s) missing. Re-run those array tasks," >&2
  echo "or set N_CHUNKS to the number you actually intend to merge." >&2
  exit 1
fi

Rscript -e '
suppressMessages(library(data.table))
br  = Sys.getenv("BRANCH", "sterile")
tag = Sys.getenv("RUN_TAG", "")
f   = sort(Sys.glob(sprintf("./data/corridor_tzero_%s%s_chunk*.rds", br, tag)))
L  = rbindlist(lapply(f, readRDS), fill = TRUE)
CK = c("tau_id","Kg","beta_ros","r_L_P")
cat(sprintf("[merge] %d files -> %s rows, %s combinations\n", length(f),
            format(nrow(L), big.mark=","), format(uniqueN(L[, ..CK]), big.mark=",")))
if (uniqueN(L, by = c(CK, "coord")) < nrow(L))
  stop("duplicate (combination, coord) rows -- two chunks covered the same slice")
saveRDS(L, sprintf("./data/corridor_tzero_%s%s.rds", br, tag))
cat(sprintf("[merge] wrote ./data/corridor_tzero_%s%s.rds\n", br, tag))
if (nzchar(tag)) {
  cat(sprintf("[merge] coordinates in this run: %s\n", paste(sort(unique(na.omit(L$coord))), collapse=", ")))
  cat("[merge] this is a PARTIAL run -- join it onto the full file on\n")
  cat("        (tau_id, Kg, beta_ros, r_L_P, coord); do not treat it as complete\n")
}
print(L[, .N, by = note][order(-N)])
'
## The exit status is load-bearing. Without this check the script announced
## "merge done" after Rscript failed to start, and the absence of an output file
## was the only clue.
rc=$?
if [ $rc -ne 0 ]; then
  echo "== MERGE FAILED (exit $rc) -- no output written ==" >&2
  exit $rc
fi
OUT="data/corridor_tzero_${BRANCH}${RUN_TAG:-}.rds"
[ -f "$OUT" ] || { echo "== MERGE FAILED: $OUT was not created ==" >&2; exit 1; }
echo "== merge done | $(ls -lh "$OUT" | awk '{print $5}') =="
