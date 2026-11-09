## ============================================================================
## lib/fig_damage_matched.R -- shared code of Fig. 4 and Extended Data Fig. 5: the damage-matched experiment
## (pipeline/40_damage_matched.R on the units of pipeline/43_damage_matched_gated_units.R,
## run by pipeline/cluster/damage_matched_gated.sh).
## ----------------------------------------------------------------------------
##   fig_dm_matched(ratio_min)   the pairs that reach the figures, in this order:
##     1. a matched status (matched_up / matched_same / matched_down) with finite
##        ON(beta*) integrals and |match_err| <= MATCH_TOL (1e-4: the bisection target)
##     2. a robust match: no jump at beta* (jump_pts <= JUMP_MAX = 10 percentage points
##        across beta* x [0.97, 1.03]) and one upward crossing of the cost profile
##        (n_cross_up == 1), so beta* does not depend on the branch of a discontinuous profile
##     3. the brake earns its keep: Tregs ON at beta_r causes at least DAM_MIN (default
##        5%) less cumulative epithelial damage than Tregs OFF at beta_r, accumulated over
##        the pair's own common window win_B = max(tconv_off, tconv_onB)
##        (D_win_onB <= (1 - DAM_MIN) D_win_off, from pipeline/44_window_damage.R). The
##        criterion is therefore read over the same kind of window as the match it selects
##        for, and it is strong enough that no pair of the set falls below beta* = 1.05 beta_r.
##        DAM_BASIS=run applies the same floor to the integrals over the whole run
##        (D_total_onB <= (1 - DAM_MIN) D_total_off), which selects a sample that overlaps
##        this one in 97.9% of the robust pairs (savings correlated at 0.994).
##   STEP 3 REPLACES TWO RETIRED CRITERIA: "Tregs act within the run" (t_div_M1 < 12,000,
##   a 1% difference in M1) and "beta* >= 1.05 beta_r". The first was a timing threshold
##   standing in for an effect size. Of the 15,272 pairs it alone removed, the brake saved
##   a median 0.04% of the cumulative damage (max 0.66%) and 15,271 were neutral at their
##   own beta_r: their Treg channel is inert (S decays 7.5x faster than in the pairs that
##   pass, so the loop gain tau_10 (sigma_S/tau_S) (tau_t_a/tau_t_d) is 81x lower), while
##   their challenge is the same (equal median K_g, r_P, beta_r, E_h(0), and a cumulative
##   damage and peak load within 10% of the passing pairs). A floor on the damage the
##   brake saves removes them in the currency the experiment matches, and subsumes the
##   ratio cut: at 5% the minimum beta*/beta_r is 1.0526.
##   The Treg effect of each unit's combination at its own beta_r (help / neutral / harm,
##   from the unit file) is merged in as `treg_effect` for description only.
## ============================================================================
suppressMessages({library(data.table); library(deSolve); library(ggplot2); library(patchwork)})

FIG_DM_DATA  = Sys.getenv("FIG_DM_DATA", "./data/fig04_damage_matched_gated.rds")
FIG_DM_UNITS = Sys.getenv("UNIT_FILE", "./data/fig04_gated_units.rds")
FIG_DM_WIN   = Sys.getenv("WIN_FILE", "./data/fig06_window_damage.rds")

fig_dm_matched = function(path = FIG_DM_DATA, units = FIG_DM_UNITS, match_mode = "conv",
                          match_tol = as.numeric(Sys.getenv("MATCH_TOL", "1e-4")),
                          jump_max = as.numeric(Sys.getenv("JUMP_MAX", "10")),
                          dam_min = as.numeric(Sys.getenv("DAM_MIN", "0.05")),
                          dam_basis = Sys.getenv("DAM_BASIS", "window"),
                          win_file = FIG_DM_WIN, verbose = TRUE) {
  R = as.data.table(readRDS(path))[match == match_mode]
  U = as.data.table(readRDS(units))[, .(i = .I, treg_effect)]
  R = merge(R[, setdiff(names(R), "treg_effect"), with = FALSE], U, by = "i")
  M0 = R[status %in% c("matched_up", "matched_same", "matched_down")]
  M1 = M0[is.finite(auc_P_onS) & is.finite(P_win_onS) & abs(match_err) <= match_tol]
  M2 = M1[is.finite(jump_pts) & jump_pts <= jump_max & n_cross_up == 1]
  stopifnot(dam_basis %in% c("window", "run"))
  if (dam_basis == "window") {
    if (!file.exists(win_file)) stop(win_file, " not found -- run: N_CORES=9 Rscript pipeline/44_window_damage.R")
    W = as.data.table(readRDS(win_file))[, .(i, win_B, D_win_off, D_win_onB, save_win)]
    M2 = merge(M2, W, by = "i")
    M  = M2[save_win >= dam_min]
  } else {
    M  = M2[D_total_onB <= (1 - dam_min) * D_total_off]
  }
  ## ---- exclusion of combinations whose Tregs-ON threshold inverts ----------------
  ## Author's ruling 2026-11-06. In these combinations the Tregs-OFF arm runs
  ## increasing, as it should, but the Tregs-ON arm runs decreasing: with Tregs the
  ## low-burden attractor is reached from the MORE injured initial conditions. A
  ## regime in which a larger initial injury gives the better outcome is outside the
  ## range of behaviour this work interprets, so it is removed rather than reported.
  ## The runs are numerically clean (one crossing per arm, no solver note), so this is
  ## a modelling judgement, not a repair.
  ## Kept kinds: `measured` (single increasing Tregs-ON crossing) and
  ## `on_lo_basin_throughout` (Tregs-ON low-burden over the whole line). This is the
  ## same filter lib/fig_skeleton.R:54 applies to BI/MB, so the counts reported in the
  ## text agree with the sample sizes shown in every figure.
  ## EXCLUDE_DECREASING=0 restores the pre-2026-11-06 population.
  if (Sys.getenv("EXCLUDE_DECREASING", "1") == "1") {
    KK = c("batch", "row_ind", "Kg", "beta_ros", "r_L_P")
    UD = as.data.table(readRDS(units))[, i := .I][, beta_ros := b_off]
    SEL = unique(as.data.table(readRDS("./data/selected_big.rds"))[
            , c(KK, "rescue_gap_kind"), with = FALSE])
    UD = merge(UD, SEL, by = KK, all.x = TRUE)
    OK = UD[rescue_gap_kind %in% c("measured", "on_lo_basin_throughout"), i]
    n_pre = nrow(M); M = M[i %in% OK]
    if (verbose) cat(sprintf(paste0("[fig04 data] EXCLUDE_DECREASING: %d of %d pairs kept ",
                                    "(%d dropped for an inverted or unlabelled Tregs-ON threshold)\n"),
                             nrow(M), n_pre, n_pre - nrow(M)))
  }

  if (verbose) cat(sprintf(paste0("[fig04 data] %s: %d units; %d matched; %d within |match_err| <= %g; %d robust (no jump, one crossing); ",
                                  "%d in which the brake saves at least %g%% of the cumulative damage over the %s (%d up, %d same, %d down); ",
                                  "beta*/beta_r min %.4f median %.2f\n"),
                           basename(path), nrow(R), nrow(M0), nrow(M1), match_tol, nrow(M2), nrow(M), 100 * dam_min,
                           if (dam_basis == "window") "common window" else "whole run",
                           sum(M$status == "matched_up"), sum(M$status == "matched_same"), sum(M$status == "matched_down"),
                           min(M$b_star / M$b_off), median(M$b_star / M$b_off)))
  M[]
}
