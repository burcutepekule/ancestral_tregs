SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## 02_assemble.R -- STAGE 0: ingest the v6 sweep, assign the combination key,
## merge near-degenerate fixed points, and say whether the run is worth analysing.
##
## READS .rds, NOT .csv. 01_sweep.R stores the complete fixed-point set as
## list-columns (fp_off_<state>, fp_on_<state>, one vector per combination),
## which cannot survive a CSV round trip.
##
## in : <DATA_DIR>/*.rds          (every chunk file the sweep wrote for ONE batch)
## out: ./data/<OUT_STEM>_batch_<BATCH>.rds
##      ./data/<OUT_STEM>_batch_<BATCH>_untrusted.rds   (the excluded combinations)
##      OUT_STEM = fp_big (main sweep, default) or another stem for a cell-list re-sweep
## RUN: DATA_DIR=./data/sweep_raw/fp_batch_1 BATCH=1 Rscript pipeline/02_assemble.R
##
## ----------------------------------------------------------------------------
## THE KEY. This is the one thing this script exists to get right.
##
## row_ind is assigned PER TAU GRID -- 00_parameter_grids.R does tau_grid[, row_ind := .I]
## -- so row_ind 1863 in tau_grid_v6_1 and in tau_grid_v6_2 are DIFFERENT parameter
## sets. The sweep records only row_ind, and the combination grid CJ(Kg, beta_ros,
## r_L_P) is identical in every batch, so once two batches are pooled the collision
## is near-total rather than rare.
##
## v5 handled this by stamping `batch` here and keying downstream on the FOUR-part
## tuple (batch, row_ind, Kg, beta_ros, r_L_P). That works and it failed repeatedly,
## always the same way: a merge, unique() or dcast() that omitted `batch` was still
## syntactically valid and silently conflated two parameter sets -- dcast in
## particular applies fun.aggregate = length without saying so.
##
## v6 removes the droppable element instead of documenting it:
##       tau_id = 10000 * batch + row_ind
## one integer, globally unique across grids, batch baked in and impossible to
## forget. batch and row_ind are kept as provenance columns; nothing keys on them.
##
## WHERE IT IS ASSIGNED CHANGED. tau_id used to be built HERE, which left every
## consumer that wanted the parameter row to invert it as tau_id %/% TAU_STRIDE
## against a stride this script owned -- three scripts did exactly that, with the
## stride hard-coded. 00_parameter_grids.R now stamps batch, tau_id and tau_stride
## into the grid itself and the sweep carries them through, so the key is a
## property of the data. This script only VERIFIES the arithmetic still holds;
## the check below is what makes tau_id safe, and a grid of 10000+ rows would be
## caught by the assertion in 00_parameter_grids.R before it could collide here.
## ----------------------------------------------------------------------------
suppressMessages(library(data.table))
source("./lib/plp_clean.R")

DATA_DIR  = Sys.getenv("DATA_DIR", "./data/sweep_raw/fp_batch_1")
BATCH     = as.integer(Sys.getenv("BATCH", "1"))
MERGE_TOL = as.numeric(Sys.getenv("MERGE_TOL", "1"))     # Chebyshev on (Eh, PLP);
                                                        # must equal the sweep's
                                                        # RAY_CLUST_TOL, asserted below
TAU_STRIDE = 10000L                                      # see the key note above
OUT_STEM  = Sys.getenv("OUT_STEM", "fp_big")
OUT       = sprintf("./data/%s_batch_%d.rds", OUT_STEM, BATCH)
## CHECK 1 below flags any beta_ros above BETA_MAX: 0.01 for the main grid (which
## tops out at 0.0095), or the largest level of a cell-list re-sweep.
BETA_MAX  = as.numeric(Sys.getenv("BETA_MAX", "0.01"))
STATES    = c("PL","PLP","CLP","Eh","Ea","En","M0","M1","R","T0","Tr","D","S")

## Glob *.rds, not a name pattern. The directory (fp_v6_<batch>) is what identifies
## the run; the flush filenames inside it are an implementation detail of the sweep
## and have already changed once. Anything else in this directory is a bug anyway.
f = sort(Sys.glob(file.path(DATA_DIR, "*.rds")))
if (!length(f)) stop("no .rds under ", DATA_DIR)
cat(sprintf("[read] %d chunk files under %s\n", length(f), DATA_DIR))

## Tolerate a chunk that is mid-write. The sweep writes a NEW file per flush, so
## completed files are never touched again -- but if this runs while the sweep is
## still going, the file being written right now will fail to deserialise and would
## otherwise abort the whole ingest. Skip it and say so; the next ingest picks it
## up. This is what makes it safe to read partial results during a long run.
.parts = lapply(f, function(x) tryCatch(readRDS(x), error = function(e) NULL))
.bad   = f[vapply(.parts, is.null, logical(1))]
if (length(.bad)) {
  cat(sprintf("[read] SKIPPED %d unreadable chunk file(s) -- being written right now?\n",
              length(.bad)))
  for (b in .bad) cat("        ", b, "\n")
}
S = rbindlist(.parts[!vapply(.parts, is.null, logical(1))], fill = TRUE)
if (!nrow(S)) stop("no readable chunks under ", DATA_DIR)
cat(sprintf("[read] %s rows, %d columns  (from %d of %d files)\n",
            format(nrow(S), big.mark=","), ncol(S), length(f) - length(.bad), length(f)))

## ---- THE KEY, NOW CARRIED RATHER THAN ASSIGNED ------------------------------
## batch and tau_id used to be stamped here, which meant every consumer that
## wanted the parameter row had to invert tau_id by arithmetic against a stride
## this script owned. They are now stamped into the tau grids by 00_parameter_grids.R
## and carried through the sweep, so this script only VERIFIES them. A run whose
## chunks lack the columns predates that change and must not be pooled with one
## that has them, because the two would be keyed on different conventions.
if (!all(c("batch", "tau_id") %in% names(S)))
  stop("this run predates the stamped key (no batch/tau_id column) -- rerun ",
       "00_parameter_grids.R and the sweep, or ingest it with the old 02_assemble.R")
if (any(S$batch != BATCH))
  stop("chunks under ", DATA_DIR, " carry batch ",
       paste(sort(unique(S$batch)), collapse = ","), " but BATCH=", BATCH)
if (!isTRUE(all(S$tau_id == TAU_STRIDE * S$batch + S$row_ind)))
  stop("tau_id does not agree with TAU_STRIDE * batch + row_ind -- the grid and ",
       "this script disagree about the stride; do not pool this batch")
CK = c("tau_id", "Kg", "beta_ros", "r_L_P")

n0 = nrow(S); S = unique(S, by = CK)
if (nrow(S) < n0) cat(sprintf("[read] dropped %d duplicate rows (requeued tasks)\n", n0 - nrow(S)))
S[, patho := r_L_P > 0]
cat(sprintf("[key ] tau_id carried from the grid -> %s..%s over %s tau rows\n",
            min(S$tau_id), max(S$tau_id), format(uniqueN(S$tau_id), big.mark=",")))

## ---- THE SETTINGS THE SWEEP RAN UNDER ---------------------------------------
## The sweep stamps these per row. Two things are checked rather than trusted to
## a comment: that a single batch was produced under ONE set of settings, and
## that MERGE_TOL below equals the sweep's RAY_CLUST_TOL. If they differ, the
## labels (map_*, n_cross_*, Eh0_*) and the fixed-point set this script rewrites
## describe two different clusterings, and nothing downstream can tell.
.stamp = c("ray_clust_tol","ext_tol","atol","rtol","hmax","tmax","maxsteps")
.stamp = .stamp[.stamp %in% names(S)]
if (length(.stamp)) {
  .mix = .stamp[vapply(.stamp, function(v) uniqueN(S[[v]]) > 1L, logical(1))]
  if (length(.mix))
    stop("chunks under ", DATA_DIR, " mix settings across rows: ",
         paste(.mix, collapse = ", "), " -- this batch pools two different sweeps")
  cat("[stmp]", paste(sprintf("%s=%g", .stamp, vapply(.stamp, function(v) as.numeric(S[[v]][1]), numeric(1))),
                      collapse = "  "), "\n")
  if ("ray_clust_tol" %in% .stamp && !isTRUE(all.equal(MERGE_TOL, S$ray_clust_tol[1])))
    stop("MERGE_TOL = ", MERGE_TOL, " but the sweep clustered at ray_clust_tol = ",
         S$ray_clust_tol[1], " -- set MERGE_TOL to match, or the merged fixed-point ",
         "set and the sweep's own labels describe different clusterings")
}

## ============================================================================
## MERGE near-degenerate fixed points
## ----------------------------------------------------------------------------
## The sweep writes every fixed point the finder returned and merges nothing, so
## clustering is a decision made HERE, where it is reversible: raise MERGE_TOL,
## rerun this script, keep the sweep. Members of a run separated by <= MERGE_TOL
## on the composite become one point, positioned at the MEAN of the members in
## every state. Singletons pass through untouched.
##
## NOTE what this does NOT do: the sweep's labels (map_*, n_cross_*, Eh0_*) were
## computed against the UNMERGED set, using its own RAY_CLUST_TOL. Merging here at
## a different tolerance remaps the map indices but does not relabel or re-bisect.
## MERGE_TOL is now asserted equal to the sweep's ray_clust_tol above rather than
## left to a reader to honour.
##
## THE RULE IS THE SWEEP'S RULE. This was a gap in the composite burden, and the
## composite is a projection: two states that trade epithelial damage against
## pathogen load can share a value, so it merged attractors that were 27 apart in
## P_LP. It is now the same Chebyshev test the sweep clusters with -- |dEh| and
## |dPLP| both within tolerance -- at the same tolerance, asserted rather than
## trusted. One rule, no special cases. Single-linkage by union-find, because a
## Chebyshev ball is not transitive and in two dimensions there is no sort order
## to break on. Clusters are numbered by increasing composite so id 1 remains the
## lowest-burden attractor and the map remap below keeps its meaning.
##
## WHAT TO WATCH. With no extinction clause the tolerance alone decides whether a
## pathogen-free fixed point merges with one carrying P_LP <= MERGE_TOL, i.e. up
## to 1e11 organisms under the model's conversion. Whether that ever happens is a
## measurement, not an assumption: n_merge_cleared below counts the merges that
## joined a cleared fixed point to a carrying one, per condition. If it is zero,
## the question was moot and the methods can say so. If it is not, decide
## deliberately whether those are one outcome or two.
## ============================================================================
merge_one = function(cmpv, ehv, plpv) {    # -> integer cluster id per fixed point
  n = length(cmpv); if (!n) return(integer(0))
  if (is.null(ehv) || is.null(plpv) || length(ehv) != n || length(plpv) != n)
    stop("merge_one: need fp_*_Eh and fp_*_PLP of the same length as fp_*_cmp")
  parent = seq_len(n)
  root = function(i) { while (parent[i] != i) i = parent[i]; i }
  if (n > 1L) for (i in seq_len(n - 1L)) for (j in (i + 1L):n)
    if (abs(ehv[i] - ehv[j]) <= MERGE_TOL && abs(plpv[i] - plpv[j]) <= MERGE_TOL) {
      ri = root(i); rj = root(j); if (ri != rj) parent[rj] = ri
    }
  lab = vapply(seq_len(n), root, integer(1))
  u = unique(lab)
  u = u[order(vapply(u, function(k) min(cmpv[lab == k]), numeric(1)))]
  match(lab, u) }

## Does any cluster contain BOTH a cleared and a carrying fixed point?
merge_joins_cleared = function(g, plpv) {
  if (!length(g)) return(0L)
  sum(vapply(split(plpv, g), function(z) any(z <= 0) && any(z > 0), logical(1))) }

merge_cond = function(S, pre) {
  cc = paste0("fp_", pre, "_")
  cmpcol = paste0(cc, "cmp"); if (!cmpcol %in% names(S)) { cat("[merge] no", cmpcol, "\n"); return(S) }
  ehcol = paste0(cc, "Eh"); plpcol = paste0(cc, "PLP")
  if (!all(c(ehcol, plpcol) %in% names(S)))
    stop("[merge] ", pre, ": the Chebyshev rule needs ", ehcol, " and ", plpcol)
  .g1 = function(i) merge_one(S[[cmpcol]][[i]], S[[ehcol]][[i]], S[[plpcol]][[i]])
  ng  = vapply(seq_len(nrow(S)), function(i) { g = .g1(i); if (!length(g)) 0L else max(g) }, integer(1))
  nfp = vapply(S[[cmpcol]], length, integer(1))
  hit = which(nfp > ng & ng > 0L)                     # rows where something merges
  cat(sprintf("[merge] %s: %s of %s rows have >=1 near-degenerate pair (tol %g)\n",
              pre, format(length(hit), big.mark=","), format(nrow(S), big.mark=","), MERGE_TOL))
  ## THE DROPPED EXTINCTION CLAUSE, MADE VISIBLE. Counted over every row, not just
  ## the merging ones, so a zero here is evidence and not an empty subset.
  .jc = sum(vapply(seq_len(nrow(S)),
                   function(i) merge_joins_cleared(.g1(i), S[[plpcol]][[i]]) > 0L, logical(1)))
  cat(sprintf("[merge] %s: %s rows where a cluster joins a CLEARED fixed point to a carrying one\n",
              pre, format(.jc, big.mark=",")))
  if (.jc) cat(sprintf("[merge] %s: at MERGE_TOL = %g that merges P_LP = 0 with up to %g units (~%.0e organisms). Decide deliberately.\n",
                       pre, MERGE_TOL, MERGE_TOL, MERGE_TOL * 1e11))
  if (!length(hit)) return(S)

  vcols = c(STATES, "cmp", "max_Re_eig"); vcols = vcols[paste0(cc, vcols) %in% names(S)]
  for (i in hit) {
    g = .g1(i)
    for (v in vcols) {
      x = S[[paste0(cc, v)]][[i]]
      set(S, i, paste0(cc, v), list(as.numeric(tapply(x, g, mean, na.rm = TRUE))))
    }
    ## type: a cluster is "stable" only if every member is
    tcol = paste0(cc, "type")
    if (tcol %in% names(S)) {
      tp = S[[tcol]][[i]]
      set(S, i, tcol, list(unname(tapply(tp, g, function(z) if (all(z == "stable")) "stable" else "mixed"))))
    }
    ## remap the map so each probed Eh0 points at its cluster, not its old index
    mcol = paste0("map_", pre)
    if (mcol %in% names(S) && !is.na(S[[mcol]][i])) {
      ch  = strsplit(S[[mcol]][i], "")[[1]]
      new = ifelse(ch == ".", ".", as.character(g[suppressWarnings(as.integer(ch))]))
      set(S, i, mcol, paste(new, collapse = ""))
    }
    set(S, i, paste0(cc, "n"), max(g))
  }
  S }

S = merge_cond(S, "off"); S = merge_cond(S, "on")
S[, `:=`(n_fp_off = vapply(fp_off_cmp, length, integer(1)),
         n_fp_on  = vapply(fp_on_cmp,  length, integer(1)))]

## ---- negative P_LP: same rule as everywhere else -----------------------------
for (.cn in setdiff(grep("PLP", names(S), value = TRUE), c("PLP0","fp_off_PLP","fp_on_PLP")))
  S = plp_clean(S, plp = .cn, comp = NULL, mode = "cell", label = .cn)

## ---- map_trusted: is this combination's map usable? --------------------------
## Everything downstream indexes map_* into the fp_* vectors, which is only
## meaningful if the final states really are AT the fixed points they were assigned
## to. dev_med is the typical probed Eh0's relative discrepancy. Gate on this
## rather than dropping rows outright at the flag stage: the combinations are real,
## only their MAP is unreliable.
##
## ITS BLIND SPOT, which matters in v6: a final state resting exactly ON an
## unstable fixed point has dev = 0 and passes this gate. That is why the sweep
## checks stability when labelling and does not rely on dev. Do not treat a small
## dev as evidence that every probed Eh0 landed in a basin -- read n_unmatched_*
## and the '.' characters in map_* for that.
MAP_DEV_TOL = as.numeric(Sys.getenv("MAP_DEV_TOL", "1e-2"))
S[, `:=`(map_trusted_off = is.finite(dev_med_off) & dev_med_off <= MAP_DEV_TOL,
         map_trusted_on  = is.finite(dev_med_on)  & dev_med_on  <= MAP_DEV_TOL)]
S[, map_trusted := map_trusted_off & map_trusted_on]
cat(sprintf("[trust] map_trusted (dev_med <= %g): both %.2f%% | off %.2f%% | on %.2f%%\n",
            MAP_DEV_TOL, 100*mean(S$map_trusted), 100*mean(S$map_trusted_off), 100*mean(S$map_trusted_on)))

## ---- EXCLUDE the untrusted combinations --------------------------------------
## Not merely flagged. When the typical final state is far from any fixed point the
## finder returned, the attractor framing fails for that combination outright: the
## endpoint-seeded finder was SEEDED from those states, so the fixed points it
## returned may be genuine solutions of f(x) = 0 while having nothing to do with
## where the line actually goes -- fp_* is contaminated too, not just map_*.
## They are written to a separate file rather than discarded, so the exclusion is
## auditable and the limit-cycle candidates can be studied on their own terms.
##
## WATCH THE BREAKDOWN. If the excluded combinations cluster in one corner of the
## grid rather than scattering, dropping them biases that corner and the exclusion
## has to be stated in the methods, not just performed.
EXCLUDE_UNTRUSTED = !identical(tolower(Sys.getenv("EXCLUDE_UNTRUSTED", "TRUE")), "false")
BAD = S[map_trusted == FALSE]
if (nrow(BAD)) {
  cat(sprintf("\n================ EXCLUDED: untrusted map (dev_med > %g) ================\n", MAP_DEV_TOL))
  cat(sprintf("%s of %s combinations (%.2f%%)  |  sterile %s, pathogenic %s\n",
              format(nrow(BAD), big.mark=","), format(nrow(S), big.mark=","),
              100*nrow(BAD)/nrow(S),
              format(BAD[patho==FALSE,.N], big.mark=","), format(BAD[patho==TRUE,.N], big.mark=",")))
  ## per-level exclusion RATE, not counts: a count is confounded by how many
  ## combinations each level has, and the rate is what reveals clustering.
  for (v in c("Kg","beta_ros","r_L_P")) {
    tot = S[, .N, by = v]; ex = BAD[, .N, by = v]
    m = merge(tot, ex, by = v, all.x = TRUE, suffixes = c("_all","_excl"))
    m[is.na(N_excl), N_excl := 0L][, pct := round(100*N_excl/N_all, 2)]
    setorderv(m, v)
    cat(sprintf("\n-- by %s --\n", v)); print(m[])
  }
  cat("\nA flat rate across all three means the loss is unbiased. A rate that rises\nmonotonically along any axis means it is not, and must be reported.\n")
  saveRDS(BAD, sub("\\.rds$", "_untrusted.rds", OUT))
  cat(sprintf("[trust] wrote %s for inspection\n", sub("\\.rds$", "_untrusted.rds", OUT)))
}
if (EXCLUDE_UNTRUSTED) {
  S = S[map_trusted == TRUE]
  cat(sprintf("[trust] EXCLUDED -- %s combinations carried forward\n", format(nrow(S), big.mark=",")))
} else {
  cat("[trust] EXCLUDE_UNTRUSTED=false -- untrusted combinations KEPT; filter downstream yourself\n")
}

saveRDS(S, OUT)
cat(sprintf("[read] saved %s\n\n", OUT))

## ============================================================================
## CHECKS
## ============================================================================
cat("================ coverage ================\n")
cat(sprintf("tau rows %d | combos/row %.1f | Kg %s | beta_ros %g..%g | r_L_P %g..%g | PLP0 %s\n",
            uniqueN(S$tau_id), nrow(S)/uniqueN(S$tau_id),
            paste(sort(unique(S$Kg)), collapse=","), min(S$beta_ros), max(S$beta_ros),
            min(S$r_L_P), max(S$r_L_P), paste(unique(S$PLP0), collapse=",")))
cat(sprintf("sterile %s | pathogenic %s\n", format(S[patho==FALSE,.N], big.mark=","),
            format(S[patho==TRUE,.N], big.mark=",")))
cat(sprintf("probe: %s values, spacing %g, bisected to %g\n",
            nchar(S$map_off[which(!is.na(S$map_off))[1]]), S$probe_by[1], S$eh0_tol[1]))

cat("\n================ CHECK 1: beta_ros scale ================\n")
## braces are load-bearing: at top level a bare `else` on a new line is a parse
## error, because the `if` is already a complete expression by then.
if (max(S$beta_ros) > BETA_MAX) {
  cat(sprintf("*** FAIL max %g above BETA_MAX = %g -- wrong grid?\n", max(S$beta_ros), BETA_MAX))
} else {
  cat(sprintf("PASS max %g (at or below BETA_MAX = %g)\n", max(S$beta_ros), BETA_MAX))
}

cat("\n================ CHECK 2: does the epithelium heal? ================\n")
lowEh = vapply(seq_len(nrow(S)), function(i) { v = S$fp_off_Eh[[i]]; c = S$fp_off_cmp[[i]]
  if (!length(v) || !length(c)) NA_real_ else v[which.min(c)] }, numeric(1))
q = quantile(lowEh, c(0,.25,.5,.75,.95,1), na.rm = TRUE); print(round(q,2))
if (is.na(q[["50%"]])) {
  cat("*** all NA -- no fixed point recovered anywhere; stop and investigate.\n")
} else if (q[["50%"]] < 30) {
  cat(sprintf("*** WARN median %.2f -- the low-burden state is damaged tissue, not healed tissue.\n", q[["50%"]]))
} else {
  cat(sprintf("PASS median %.2f -- the low-burden state is genuinely healthier.\n", q[["50%"]]))
}

cat("\n================ CHECK 3: final state vs matched fixed point ================\n")
## dev_max is the worst probed Eh0, so ONE unconverged value flags a whole
## combination. dev_med is the typical one. Read them together:
##   med tiny + max large -> a few had not settled at TMAX; the map is sound
##   med large            -> final states are systematically not at the fixed
##                           points they were assigned to, and the map is fiction
for (cd in c("off","on")) {
  dm = S[[paste0("dev_max_", cd)]]; dd = S[[paste0("dev_med_", cd)]]
  if (is.null(dm)) next
  cat(sprintf("-- %s --\n", cd))
  cat("  dev_max: "); print(signif(quantile(dm, c(.5,.95,1), na.rm=TRUE), 3))
  if (!is.null(dd)) { cat("  dev_med: "); print(signif(quantile(dd, c(.5,.95,1), na.rm=TRUE), 3)) }
  cat(sprintf("  max > 1e-2 in %.2f%% of combinations", 100*mean(dm > 1e-2, na.rm=TRUE)))
  if (!is.null(dd)) cat(sprintf("  |  med > 1e-2 in %.2f%%", 100*mean(dd > 1e-2, na.rm=TRUE)))
  cat("\n")
}

cat("\n================ CHECK 4: fixed point counts ================\n")
cat("after merging, OFF:\n"); print(table(S$n_fp_off, useNA="ifany"))
cat("after merging, ON:\n");  print(table(S$n_fp_on,  useNA="ifany"))
cat("\nreached from the Eh0 line (n_ray) vs existing in the field (n_field_stable):\n")
print(S[, .N, by = .(n_ray_off, n_field_stable_off)][order(n_ray_off, n_field_stable_off)])
cat("\nn_ray < n_field_stable is NOT an error -- it is an attractor that exists and\n",
    "that the Eh0 line never enters the basin of. OFF: ", 
    sprintf("%.2f%%", 100*mean(S$n_ray_off < S$n_field_stable_off, na.rm=TRUE)),
    " | ON: ", sprintf("%.2f%%", 100*mean(S$n_ray_on < S$n_field_stable_on, na.rm=TRUE)), "\n", sep="")

cat("\n================ CHECK 5: notes and directions ================\n")
cat("note_off:\n"); print(sort(table(S$note_off), decreasing=TRUE))
cat("note_on:\n");  print(sort(table(S$note_on),  decreasing=TRUE))
cat("\ndirection of the located threshold (blank = none located):\n")
print(S[, .N, by = .(dir_off, dir_on)][order(-N)])

cat("\n================ CHECK 6: crossings, and the key consistency test ================\n")
print(S[, .N, by = .(n_cross_off, n_cross_on)][order(-N)])
## With endpoint-only seeding these two statements are equivalent by construction:
## one reached attractor means every probed Eh0 is in the same basin, so the label
## cannot change. If this ever fails, the labelling is wrong -- not the data.
.bad1 = S[(n_ray_off == 1L) != (n_cross_off == 0L), .N]
.bad2 = S[(n_ray_on  == 1L) != (n_cross_on  == 0L), .N]
if (.bad1 + .bad2 == 0) {
  cat("PASS  n_ray == 1  <=>  n_cross == 0  holds for every combination, both conditions.\n")
} else {
  cat(sprintf("*** FAIL %d OFF and %d ON combinations violate n_ray==1 <=> n_cross==0.\n", .bad1, .bad2))
  cat("    The labelling in 01_sweep.R is inconsistent with the attractor count.\n")
}
cat(sprintf("\nmore than one crossing: OFF %.3f%% | ON %.3f%% -- these get Eh0_star = NA\n",
            100*mean(S$n_cross_off > 1, na.rm=TRUE), 100*mean(S$n_cross_on > 1, na.rm=TRUE)))

cat("\n================ CHECK 7: the bisection ================\n")
## Eh0_*_coarse is the probe-resolution threshold, the value v5 would have
## reported. The bisected one must sit within half a probe spacing of it; a larger
## disagreement means the bisection walked into a different crossing.
for (cd in c("off","on")) {
  fine = S[[paste0("Eh0_", cd)]]; crs = S[[paste0("Eh0_", cd, "_coarse")]]
  k = is.finite(fine) & is.finite(crs)
  cat(sprintf("-- %s -- located %s of %s (%.1f%%)", cd, format(sum(is.finite(fine)), big.mark=","),
              format(nrow(S), big.mark=","), 100*mean(is.finite(fine))))
  if (any(k)) {
    d = abs(fine[k] - crs[k]); lim = S$probe_by[1] / 2
    cat(sprintf(" | |fine - coarse| med %.3f max %.3f", median(d), max(d)))
    cat(if (max(d) <= lim + 1e-9) sprintf(" | PASS (<= %.1f)\n", lim)
        else sprintf(" | *** FAIL %d exceed %.1f\n", sum(d > lim), lim))
  } else cat("\n")
  ## over the arms that actually bisected. Including the ones with no crossing
  ## would put the median at 0 and say nothing.
  .nb = S[[paste0("n_bisect_", cd)]][is.finite(fine)]
  cat(sprintf("   bisection steps, over the %s located arms: median %s, max %s\n",
              format(length(.nb), big.mark=","),
              if (length(.nb)) median(.nb, na.rm=TRUE) else NA,
              if (length(.nb)) max(.nb, na.rm=TRUE) else NA))
}

cat("\n================ CHECK 8: the rescue gap ================\n")
print(S[, .(n = .N, med = round(median(rescue_gap, na.rm=TRUE), 3),
            pct_pos = round(100*mean(rescue_gap > 0, na.rm=TRUE), 1),
            pct_neg = round(100*mean(rescue_gap < 0, na.rm=TRUE), 1)),
        by = .(rescue_gap_kind, patho)][order(rescue_gap_kind, patho)])
cat("\nkinds: measured = both arms located and running the SAME direction;\n",
    "on_lo_basin_throughout = the ON line reaches the low basin everywhere;\n",
    "opposite_directions = both located but crossing in opposite senses, which is\n",
    "NOT a rescue and is deliberately NA; none = at least one arm located nothing.\n", sep="")

cat("\n[00_read_first] done. MERGE_TOL was", MERGE_TOL,
    "-- rerun with a different value to re-cluster.\n")
