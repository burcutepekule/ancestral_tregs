SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## 03_select.R -- STAGE 1: gate the Tregs-OFF bistable substrate, classify what
## Tregs did to it, split the branches, and emit what the figures read.
##
## in : ./data/<IN_STEM>_batch_<b>.rds  (from 02_assemble.R, one per batch;
##                                       IN_STEM = fp_big by default)
## out: <OUT_SEL> = ./data/selected_big.rds (everything that passed the gate, all columns)
##      ./data/skeleton_<BRANCH>.rds     (Fig. 3, Extended Data Figs. 1 and 2, Supplementary Fig. 2 and others)
##      ./data/rescue_<BRANCH>.rds       (the t=0 corridor stage)
## RUN: BRANCH=sterile         Rscript pipeline/03_select.R
##      BRANCH=pathogenic_help Rscript pipeline/03_select.R
##      BRANCH=pathogenic_harm Rscript pipeline/03_select.R
##   gate only, no branch files (another sweep stem): IN_STEM=<stem> OUT_SEL=<file> BRANCH=none
##
## This script does no integration and no root-finding. Every quantity below
## already exists in the sweep output; this selects, derives two labels, renames
## into the schema the figures expect, and writes. It runs in seconds.
## ============================================================================
suppressMessages(library(data.table))

## ---- THE KEY ----------------------------------------------------------------
## tau_id is assigned by 02_assemble.R as TAU_STRIDE * batch + row_ind and is
## globally unique across tau grids. batch and row_ind ride along as provenance;
## nothing keys on them. This is the whole reason the key is one integer: in v5
## the key was the tuple (batch, row_ind, Kg, beta_ros, r_L_P) and every merge or
## dcast that omitted `batch` was still valid R and silently conflated two
## parameter sets.
CK = c("tau_id", "Kg", "beta_ros", "r_L_P")

if (!exists("CUTOFF")) CUTOFF = as.numeric(Sys.getenv("GAP_CUT", "20"))
if (!exists("BRANCH")) BRANCH = Sys.getenv("BRANCH", "sterile")
EFFECT_TOL = as.numeric(Sys.getenv("EFFECT_TOL", "1.0"))  # composite units; below this a shift is neutral

## BRANCH may carry a Treg-effect suffix: "sterile" | "pathogenic" |
## "pathogenic_help" | "pathogenic_harm". The pathogenic branch mixes two
## opposite populations -- combinations where Tregs LOWER the burden of the
## low-burden attractor and combinations where they RAISE it. Pooled, their
## distributions are bimodal and the medians cancel to ~0, so every summary on
## them is an average of two opposite effects. Splitting is not cosmetic.
BASE   = sub("_(harm|help)$", "", BRANCH)
EFFECT = if (grepl("_(harm|help)$", BRANCH)) sub("^.*_", "", BRANCH) else NA_character_

## ---- read every batch --------------------------------------------------------
## The batch list is DERIVED from what is on disk. A hard-coded count silently
## drops a batch that has not been ingested yet, and the loss is invisible.
IN_STEM = Sys.getenv("IN_STEM", "fp_big")
OUT_SEL = Sys.getenv("OUT_SEL", "./data/selected_big.rds")
bf = sort(Sys.glob(sprintf("./data/%s_batch_*.rds", IN_STEM)))
bf = bf[grepl(sprintf("%s_batch_[0-9]+\\.rds$", IN_STEM), bf)]      # exclude *_untrusted.rds
if (!length(bf)) stop("no ./data/", IN_STEM, "_batch_<n>.rds -- run 02_assemble.R first")
S = rbindlist(lapply(bf, readRDS), fill = TRUE)
cat(sprintf("[select] %d batch file(s): %s\n", length(bf), paste(basename(bf), collapse = ", ")))
cat(sprintf("[select] %s combinations, %s tau rows\n",
            format(nrow(S), big.mark = ","), format(uniqueN(S$tau_id), big.mark = ",")))

if (!"tau_id" %in% names(S))
  stop("no tau_id column -- these batch files predate the v6 key. Re-run 02_assemble.R.")
if (uniqueN(S, by = CK) < nrow(S))
  stop("CK is not unique across the pooled batches -- tau_id has collided. ",
       "Check TAU_STRIDE in 02_assemble.R and re-ingest EVERY batch.")

## The main PLP0 level only, if the sweep ever writes more than one. Levels are
## not comparable: PLP0 shifts the whole line, so the thresholds and the maps do
## too, and a level has to be chosen rather than pooled.
if ("is_main_level" %in% names(S)) S = S[is_main_level == TRUE]
S[, patho := r_L_P > 0]

## ============================================================================
## THE GATE
## ----------------------------------------------------------------------------
## Applied to the Tregs-OFF arm ONLY. This is the substrate: the systems that
## have two separated outcomes and therefore something for Tregs to rescue.
## Three conditions, all required.
##
## (i) TWO SEPARATED ATTRACTORS. hi_lo_gap_off = hi_cmp - lo_cmp > CUTOFF. This
##     IS the bistability criterion -- there is no separate "pick bistable, then
##     apply a gap". CUTOFF = 20 composite units on a ruler spanning 0-132.
##
## (ii) THE LINE CROSSES ONCE. n_cross_off == 1. v5 approximated this with an rle
##     over map_off (`n_runs_off <= 2`), because the sweep did not count crossings
##     -- it inferred them from the label string, which cannot distinguish a
##     genuine second crossing from an unassigned value in the middle. The sweep
##     now counts crossings directly, so the gate tests the quantity rather than a
##     proxy for it. Combinations with two crossings have no single threshold and
##     already carry Eh0_off = NA; this keeps them out of the substrate entirely.
##
## (iii) THE LINE RUNS INCREASING. dir_off == "increasing": the low-burden basin
##     lies ABOVE the threshold, so more healthy epithelium at t=0 means more
##     likely to resolve. Every sterile arm measured does this -- 296 of 296 in
##     the first v6 batch, no exceptions. The decreasing arms are real and are a
##     pathogenic phenomenon whose frequency rises monotonically with r_L_P, but
##     they are a DIFFERENT configuration, not noisy examples of this one:
##       * the corridor's base state is found by walking e0 along the Eh0 line
##         until the OFF system reaches the low basin. For an increasing arm that
##         lands just inside it; for a decreasing arm the same walk leaves it, so
##         every d_off, d_on and w_X measured from that base state would be
##         measured from the wrong region.
##       * the rescue gap's sign is off - on for increasing and on - off for
##         decreasing, so pooling them under one formula reports half of them
##         backwards.
##     They are analysed separately. GATE_DIR="" turns this condition off for
##     that analysis; it must not be turned off for the corridor argument.
## ============================================================================
GATE_DIR = Sys.getenv("GATE_DIR", "increasing")

## ---- THE GATE'S RULER, AND THE ONE IT WILL PROBABLY BECOME -------------------
## hi_lo_gap_off is a difference of COMPOSITES, and the gate sits on it. The sweep
## no longer clusters attractors on the composite -- it uses a Chebyshev test on
## (Eh, PLP), because the composite is a projection and merged attractors 27 apart
## in P_LP -- so the gate is now filtering on a quantity that cannot see the
## separation the clustering uses. The Chebyshev separation is computed alongside
## and both distributions are printed. The gate STILL RUNS ON THE COMPOSITE, so
## this change does not silently move the substrate under the existing figures;
## the point is to make the right CUTOFF for the Chebyshev distance readable off
## the new sweep instead of guessed before it. Switch with GATE_METRIC=cheb once
## that comparison has been made, and record which was used in the methods.
GATE_METRIC = Sys.getenv("GATE_METRIC", "composite")
if (all(c("hi_off_Eh","hi_off_PLP","lo_off_Eh","lo_off_PLP") %in% names(S))) {
  S[, hi_lo_gap_off_cheb := pmax(abs(hi_off_Eh - lo_off_Eh),
                                 abs(hi_off_PLP - lo_off_PLP))]
  cat("\n[select] separation of the two Tregs-OFF attractors, both rulers (finite rows):\n")
  print(S[is.finite(hi_lo_gap_off) & is.finite(hi_lo_gap_off_cheb),
          .(n = .N,
            composite = paste(sprintf("%.1f", quantile(hi_lo_gap_off, c(.1,.5,.9))), collapse = " / "),
            chebyshev = paste(sprintf("%.1f", quantile(hi_lo_gap_off_cheb, c(.1,.5,.9))), collapse = " / "),
            agree_at_cutoff = sum((hi_lo_gap_off > CUTOFF) == (hi_lo_gap_off_cheb > CUTOFF)))])
} else if (identical(GATE_METRIC, "cheb")) {
  stop("GATE_METRIC=cheb needs the landmark components (hi_off_Eh etc.) in the sweep output")
}
.gap_use = if (identical(GATE_METRIC, "cheb")) S$hi_lo_gap_off_cheb else S$hi_lo_gap_off
cat(sprintf("[select] gating on the %s separation\n", GATE_METRIC))
g_ok = is.finite(.gap_use) & .gap_use > CUTOFF
x_ok = is.finite(S$n_cross_off) & S$n_cross_off == 1L
d_ok = if (nzchar(GATE_DIR)) (!is.na(S$dir_off) & S$dir_off == GATE_DIR) else rep(TRUE, nrow(S))

cat(sprintf("\n[select] gate on %s combinations:\n", format(nrow(S), big.mark = ",")))
cat(sprintf("      hi_lo_gap_off > %g          : %s\n", CUTOFF, format(sum(g_ok), big.mark = ",")))
cat(sprintf("      + n_cross_off == 1          : %s\n", format(sum(g_ok & x_ok), big.mark = ",")))
cat(sprintf("      + dir_off == %-12s   : %s\n",
            if (nzchar(GATE_DIR)) sprintf('"%s"', GATE_DIR) else "(any)",
            format(sum(g_ok & x_ok & d_ok), big.mark = ",")))
SEL = S[which(g_ok & x_ok & d_ok)]
if (!nrow(SEL)) stop("nothing passed the gate -- check GAP_CUT and the direction census")
for (br in c(FALSE, TRUE))
  cat(sprintf("      %-11s %s of %s (%.1f%%)\n", if (br) "pathogenic" else "sterile",
              format(SEL[patho == br, .N], big.mark = ","), format(S[patho == br, .N], big.mark = ","),
              100 * SEL[patho == br, .N] / max(S[patho == br, .N], 1L)))

## ---- WHAT THE GATE DISCARDED, and why ----------------------------------------
## Printed rather than left implicit: if one reason dominates, the substrate is
## defined by that reason and the methods have to say so.
.why = data.table(
  reason = c("hi_lo_gap_off <= CUTOFF (not two separated attractors)",
             "n_cross_off != 1 (no threshold, or more than one)",
             sprintf("dir_off != %s", if (nzchar(GATE_DIR)) GATE_DIR else "(not applied)")),
  n = c(sum(!g_ok), sum(g_ok & !x_ok), sum(g_ok & x_ok & !d_ok)))
.why[, pct := round(100 * n / nrow(S), 2)]
cat("\n[select] discarded, in gate order (each count is what the PREVIOUS conditions left):\n")
print(.why[])

## ============================================================================
## WHAT TREGS DID TO IT
## ----------------------------------------------------------------------------
##   bistable      : ON still reaches two attractors -- the boundary moved, both remain
##   mono_B        : ON reaches only the LOW-burden attractor -- every initial
##                   condition on the line now reaches it
##   mono_A        : ON reaches only the HIGH-burden attractor -- Tregs removed the
##                   low-burden basin from the line
##   unresolved_on : no ON attractor recovered at all. NOT one of the three modes;
##                   a missing measurement, labelled so it cannot fall through to
##                   NA and silently vanish from every `rescue_mode == "..."` filter.
##
## note_on == "lo_basin_throughout" is the sweep's own statement that the ON line
## reaches the low basin from every probed Eh0, and is trusted over the composite
## comparison when the two disagree.
## ============================================================================
SEL[, rescue_mode := fifelse(!is.finite(n_ray_on), "unresolved_on",
                     fifelse(n_ray_on >= 2L, "bistable",
                      fifelse(note_on == "lo_basin_throughout", "mono_B",
                       fifelse(is.finite(lo_on) & is.finite(lo_off) & is.finite(hi_off) &
                               abs(lo_on - lo_off) <= abs(lo_on - hi_off), "mono_B", "mono_A"))))]
if (SEL[rescue_mode == "unresolved_on", .N])
  cat(sprintf("\n[select] NOTE %d combinations have no ON attractor (rescue_mode = unresolved_on);\n     carried in selected_big.rds, excluded by the figures' mode filters.\n",
              SEL[rescue_mode == "unresolved_on", .N]))
cat("\n[select] rescue_mode:\n")
print(dcast(SEL[, .N, by = .(patho, rescue_mode)], patho ~ rescue_mode, value.var = "N", fill = 0L))

## ---- NO FABRICATED hi_on -----------------------------------------------------
## v5 had a `merged_on` flag here that set chr_on := res_on when the ON map showed
## ">= 2 symbols" and chr_on was NA -- meant to catch two attractors collapsing to
## within RAY_CLUST_TOL. It never did. Its symbol count included the '.' character
## (an initial condition matching no stable fixed point), so across the whole v5
## file all 89 rows it flagged had exactly ONE attractor index and at least one
## dot, and none had two. It then wrote a HIGH-burden composite equal to the
## low-burden one, which made d_hi = hi_on - hi_off strongly negative by
## construction and those combinations look like a large Treg improvement.
## v6 fabricates nothing: when the ON line reaches one attractor, hi_on is NA and
## stays NA. Genuine near-degenerate pairs are handled in 02_assemble.R, at
## MERGE_TOL, on the full fixed-point set, where it is reversible.

## ---- Treg effect: does Tregs ON improve the reachable outcome? ---------------
## Judged on d_lo = lo_on - lo_off, the shift of the LOW-burden attractor on the
## composite. Negative = Tregs lowered the burden = help. Same sign convention as
## the skeleton figure's shift panel, so the split and the plot cannot disagree.
##
## THIS IS NOT THE SAME QUANTITY AS THE RESCUE GAP, and the two can disagree in
## the same combination. d_lo asks "is the outcome better if you get there";
## rescue_gap asks "from how much more initial damage can you get there". v5
## measured combinations labelled `help` whose rescue gap was negative. Whichever
## one a figure conditions on has to be stated.
## RECOMPUTED FROM THE COMPONENTS, NOT READ OFF A STORED SCALAR. lo_off and lo_on
## are composites the sweep happened to compute; the components (lo_off_Eh,
## lo_off_PLP, lo_on_Eh, lo_on_PLP) are recorded alongside them, so the effect
## metric is a decision made HERE, cheaply and reversibly, rather than baked into
## a nine-hour cluster job. Reweighting help and harm -- or asking what the split
## looks like on E_h alone, or on P_LP alone -- then costs one local rerun.
##
## THE WEIGHTING IS DELIBERATELY FLAT. sqrt((K_E-Eh)^2 + PLP^2) values one unit of
## epithelial damage equally with one unit of pathogen load. That is a modelling
## choice, not a neutral default: it is why Eh 38.7 / P_LP 0 and Eh 43.7 /
## P_LP 26.9 come out 1.10 apart against EFFECT_TOL = 1, i.e. five units of
## epithelium traded for twenty-seven units of pathogen scored as very nearly a
## tie. K_E = K_LP_P = 100 is what makes the flat weighting defensible -- the two
## coordinates are already on a common scale -- and the methods should say so.
##
## This is an ORDERING question and needs a scalar by construction, which is why
## the composite is kept here even though it was removed from the CLUSTERING in
## the sweep. There it was deciding identity -- how many attractors there are --
## and a projection cannot do that, because two states that trade the coordinates
## against each other share a value. Different question, different tool.
.cmpf = function(Eh, PLP) sqrt((100 - Eh)^2 + PLP^2)
if (all(c("lo_on_Eh","lo_on_PLP","lo_off_Eh","lo_off_PLP") %in% names(SEL))) {
  SEL[, lo_on_cmp_chk  := .cmpf(lo_on_Eh,  lo_on_PLP)]
  SEL[, lo_off_cmp_chk := .cmpf(lo_off_Eh, lo_off_PLP)]
  .bad = SEL[is.finite(lo_on) & is.finite(lo_on_cmp_chk) &
             abs(lo_on - lo_on_cmp_chk) > 1e-6, .N]
  if (.bad) cat(sprintf("\n[select] WARNING %d rows where lo_on disagrees with the composite of its own components\n", .bad))
  SEL[, d_lo_sel := lo_on_cmp_chk - lo_off_cmp_chk]
} else {
  cat("\n[select] NOTE landmark components absent; falling back to the stored composites\n")
  SEL[, d_lo_sel := lo_on - lo_off]
}

SEL[, treg_effect := fifelse(!is.finite(d_lo_sel), "unknown",
                      fifelse(d_lo_sel < -EFFECT_TOL, "help",
                       fifelse(d_lo_sel >  EFFECT_TOL, "harm", "neutral")))]
cat(sprintf("\n[select] treg_effect (on the low-burden attractor, tol %g composite units):\n", EFFECT_TOL))
print(dcast(SEL[, .N, by = .(patho, treg_effect)], patho ~ treg_effect, value.var = "N", fill = 0L))
cat("\n[select] treg_effect against the sign of the rescue gap -- these are different questions:\n")
print(SEL[rescue_gap_kind == "measured",
          .(n = .N, gap_pos = sum(rescue_gap > 0), gap_zero = sum(rescue_gap == 0),
            gap_neg = sum(rescue_gap < 0)), by = .(patho, treg_effect)][order(patho, treg_effect)])

saveRDS(SEL, OUT_SEL)
cat(sprintf("\n[select] wrote %s (%s combinations, %d columns)\n",
            OUT_SEL, format(nrow(SEL), big.mark = ","), ncol(SEL)))

## BRANCH=none: gate only. No skeleton or rescue file.
if (identical(BRANCH, "none")) quit(save = "no")

## ============================================================================
## THE BRANCH
## ============================================================================
B = SEL[patho == (BASE == "pathogenic")]
if (!is.na(EFFECT)) B = B[treg_effect == EFFECT]
if (!nrow(B)) stop("no selected combinations for BRANCH=", BRANCH)
cat(sprintf("\n[select] BRANCH=%s -> base=%s effect=%s -> %s combinations\n",
            BRANCH, BASE, if (is.na(EFFECT)) "(pooled)" else EFFECT,
            format(nrow(B), big.mark = ",")))

gk = function(n) if (n %in% names(B)) B[[n]] else NA_real_

## gap_censored: Tregs ON reaches the low basin from every probed Eh0, so Eh0_on*
## was never located and the recorded gap is a LOWER BOUND, not a measurement --
## and usually a large one, so read without the flag these look like the biggest
## rescues in the file. Taken from the sweep's own rescue_gap_kind rather than
## re-derived from a note string.
B[, gap_censored := rescue_gap_kind == "on_lo_basin_throughout"]

SK = data.table(
  B[, ..CK],
  batch = B$batch, row_ind = B$row_ind,          # provenance only; CK is tau_id
  rescue_mode = B$rescue_mode, treg_effect = B$treg_effect,
  rescue_gap = B$rescue_gap, rescue_gap_kind = B$rescue_gap_kind,
  gap_censored = B$gap_censored,
  ## the threshold, per arm: bisected value, the direction it runs, the
  ## probe-resolution value, and the crossing count that decides whether it exists
  Eh0_off = B$Eh0_off, Eh0_on = B$Eh0_on,
  dir_off = B$dir_off, dir_on = B$dir_on,
  Eh0_off_coarse = B$Eh0_off_coarse, Eh0_on_coarse = B$Eh0_on_coarse,
  n_cross_off = B$n_cross_off, n_cross_on = B$n_cross_on,
  note_off = B$note_off, note_on = B$note_on,
  ## reached from the line vs existing in the field -- never conflated
  n_ray_off = B$n_ray_off, n_ray_on = B$n_ray_on,
  n_field_stable_off = gk("n_field_stable_off"), n_field_stable_on = gk("n_field_stable_on"),
  on_has_hi = is.finite(B$hi_on), on_has_sad = is.finite(gk("sad_field_on")),
  ## landmark composites and their shifts
  lo_off = B$lo_off, lo_on = B$lo_on, hi_off = B$hi_off, hi_on = B$hi_on,
  sad_off = gk("sad_field_off"), sad_on = gk("sad_field_on"),
  d_lo = B$lo_on - B$lo_off, d_hi = B$hi_on - B$hi_off,
  d_sad = gk("sad_field_on") - gk("sad_field_off"),
  hi_lo_gap_off = B$hi_lo_gap_off,
  ## leading eigenvalues, and their shifts
  e_lo_off = gk("eig_lo_off"), e_lo_on = gk("eig_lo_on"),
  e_hi_off = gk("eig_hi_off"), e_hi_on = gk("eig_hi_on"),
  e_sad_off = gk("eig_sad_off"), e_sad_on = gk("eig_sad_on"),
  de_lo = gk("eig_lo_on") - gk("eig_lo_off"), de_hi = gk("eig_hi_on") - gk("eig_hi_off"),
  de_sad = gk("eig_sad_on") - gk("eig_sad_off"),
  ## P_LP of every landmark. Without these the shift panels can only show the
  ## COMPOSITE, and in the pathogenic branch the composite mixes damage and
  ## pathogen load, so a panel built on it cannot say which of the two moved.
  lo_PLP_off = gk("lo_off_PLP"), lo_PLP_on = gk("lo_on_PLP"),
  hi_PLP_off = gk("hi_off_PLP"), hi_PLP_on = gk("hi_on_PLP"),
  sad_PLP_off = gk("sad_off_PLP"), sad_PLP_on = gk("sad_on_PLP"),
  d_lo_Eh  = gk("lo_on_Eh")  - gk("lo_off_Eh"),
  d_hi_Eh  = gk("hi_on_Eh")  - gk("hi_off_Eh"),
  d_lo_PLP = gk("lo_on_PLP") - gk("lo_off_PLP"),
  d_hi_PLP = gk("hi_on_PLP") - gk("hi_off_PLP"),
  lo_Eh_off = gk("lo_off_Eh"), lo_M1_off = gk("lo_off_M1"),
  lo_Eh_on  = gk("lo_on_Eh"),  lo_M1_on  = gk("lo_on_M1"),
  hi_Eh_off = gk("hi_off_Eh"), hi_M1_off = gk("hi_off_M1"),
  hi_Eh_on  = gk("hi_on_Eh"),  hi_M1_on  = gk("hi_on_M1"),
  sad_Eh_off = gk("sad_off_Eh"), sad_M1_off = gk("sad_off_M1"),
  sad_Eh_on  = gk("sad_on_Eh"),  sad_M1_on  = gk("sad_on_M1"))
saveRDS(SK, sprintf("./data/skeleton_%s.rds", BRANCH))

## rescue_<BRANCH>.rds -- what the t=0 corridor stage reads. dir_off travels with
## it because the base state is found by walking the Eh0 line toward the low basin,
## and which way that is depends on the direction.
saveRDS(B[, c(CK, "Eh0_off", "Eh0_on", "dir_off", "dir_on", "n_cross_off", "n_cross_on",
              "rescue_mode", "treg_effect", "rescue_gap", "rescue_gap_kind",
              "gap_censored", "hi_lo_gap_off", "n_ray_on", "n_field_stable_on"),
          with = FALSE],
        sprintf("./data/rescue_%s.rds", BRANCH))
cat(sprintf("[select] wrote ./data/skeleton_%s.rds (%s combinations) and ./data/rescue_%s.rds\n",
            BRANCH, format(nrow(SK), big.mark = ","), BRANCH))
cat(sprintf("[select] saddle recovered for %.1f%% of combinations OFF (endpoint seeding could never do this)\n",
            100 * mean(is.finite(SK$sad_off))))
cat(sprintf("[select] gap_censored (ON reaches the low basin throughout): %s of %s (%.1f%%)\n",
            format(sum(SK$gap_censored, na.rm = TRUE), big.mark = ","),
            format(nrow(SK), big.mark = ","), 100 * mean(SK$gap_censored, na.rm = TRUE)))
