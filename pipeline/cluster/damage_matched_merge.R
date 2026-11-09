## ============================================================================
## pipeline/cluster/damage_matched_merge.R -- merges the task files of
## pipeline/cluster/damage_matched_gated.sh into the experiment's data file.
## Stops (writes nothing) if any of the N_UNITS units is missing, unless ALLOW_MISSING=1
## (tasks that reach their wall time leave a random subset of their units undone; the
## units are in random order, so the merged set is a uniform subsample).
## SKIP_PROFILE=1 leaves the (large) scan profiles out.
## in : $SCRATCH_BASE/fig04_gated/dm_gated_t*.rds (+ _profile.rds)
## out: ./data/fig04_damage_matched_gated.rds (+ _profile.rds)
## RUN: ALLOW_MISSING=1 SKIP_PROFILE=1 Rscript pipeline/cluster/damage_matched_merge.R
## ============================================================================
suppressMessages(library(data.table))
SRC = Sys.getenv("SRC", "."); setwd(SRC)
DIR = file.path(Sys.getenv("SCRATCH_BASE", "~/scratch"), Sys.getenv("DIR_NAME", "fig04_gated"))
STEM = Sys.getenv("STEM", "dm_gated"); OUT = Sys.getenv("OUT", "./data/fig04_damage_matched_gated.rds")
N_UNITS = as.integer(Sys.getenv("N_UNITS", "75000"))
fr = sort(Sys.glob(file.path(DIR, paste0(STEM, "_t*.rds")))); fr = fr[!grepl("_profile\\.rds$", fr)]
fp = sub("\\.rds$", "_profile.rds", fr)
## a task still running rewrites its files after every unit; a file caught mid-write
## is read again after a pause rather than stopping the merge
rd = function(f) { for (k in 1:5) { x = tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(x)) return(x); Sys.sleep(3) }
  warning("could not read ", f, " -- left out"); NULL }
R = rbindlist(lapply(fr, rd), fill = TRUE)
P = if (Sys.getenv("SKIP_PROFILE", "0") == "1") data.table() else rbindlist(lapply(fp[file.exists(fp)], rd), fill = TRUE)
R = unique(R, by = c("i", "match")); if (nrow(P)) P = unique(P, by = c("i", "match", "beta"))
miss = setdiff(seq_len(N_UNITS), R$i)
cat(sprintf("[dm merge] %d task files, %d units, %d profile rows; %d units missing\n",
            length(fr), uniqueN(R$i), nrow(P), length(miss)))
if (length(miss)) {
  tasks = sort(unique((miss - 1L) %% 400L))
  cat("missing units' tasks (resubmit with --array=", paste(tasks, collapse = ","), ")\n", sep = "")
  if (Sys.getenv("ALLOW_MISSING", "0") != "1") stop("units missing; set ALLOW_MISSING=1 to merge anyway")
}
saveRDS(R, OUT)
if (nrow(P)) saveRDS(P, sub("\\.rds$", "_profile.rds", OUT))
print(R[, .N, by = .(match, status)][order(match, -N)])
if ("jump_pts" %in% names(R)) cat(sprintf("matched with jump_pts > 10: %d\n", R[grepl("^matched", status) & jump_pts > 10, .N]))
