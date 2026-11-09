## ============================================================================
## pipeline/45_short_window.R -- the cumulative lamina-propria pathogen load of the
## Tregs-OFF arm at beta_r and the Tregs-ON arm at beta_r*, accumulated over the SHORTER
## of the two convergence times, for the pairs in which the matched window reverses the
## early pathogen advantage.
## ----------------------------------------------------------------------------
## WHY THIS EXISTS. Fig. 4 compares the two arms over the matching window, the LONGER
## of the two convergence times (pipeline/40_damage_matched.R). That window is the one the
## damage match is entitled to, but it is generous to whichever arm settles later: a small
## residual difference in P_LP, integrated over a long tail, can outweigh a large early
## reduction. In 1,567 of the 11,940 simulations both the cumulative and the steady-state
## pathogen load come out higher with Tregs, and the question is whether that is the tail
## or the acute course. Re-accumulating P_LP to min(tconv_off, tconv_onS) -- the window
## both arms have earned, where neither is integrated past its own settling -- separates
## the two. It is a diagnostic for the Results statement, not a second definition of the
## experiment: the damage match itself is untouched.
##
## in : ./data/fig04_damage_matched_gated.rds, ./data/fig04_gated_units.rds,
##      ./data/fig06_window_damage.rds, ./data/tau_grid_v6_*.rds, ./lib/*
## out: ./data/fig06_short_window.rds  one row per reversal pair: i, win_short,
##      P_short_off, P_short_onS, save_short (= 1 - P_short_onS / P_short_off)
## RUN: N_CORES=9 Rscript pipeline/45_short_window.R   (~1 min for 1,567 pairs)
## ============================================================================
suppressMessages({library(data.table); library(deSolve); library(parallel)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
assign("BRANCH", Sys.getenv("BRANCH", "pathogenic_harm"), envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R")); suppressMessages(source("./lib/extinction_floor.R"))
suppressMessages(source("./lib/fig_damage_matched.R"))
N_CORES = as.integer(Sys.getenv("N_CORES", "9"))
TT_N    = as.integer(Sys.getenv("TT_N", "600"))
OUT     = Sys.getenv("OUT", "./data/fig06_short_window.rds")

## the 11,940 of Fig. 4, then the reversal subset: cumulative AND steady-state both worse
M = as.data.table(fig_dm_matched(verbose = FALSE))
U = as.data.table(readRDS("./data/fig04_gated_units.rds"))[, .(i = .I, batch, row_ind)]
M = merge(M[, setdiff(names(M), c("batch", "row_ind")), with = FALSE], U, by = "i")
S = M[P_win_onS > P_win_off & P_end_onS > P_end_off]
stopifnot(nrow(M) == 11940L)
cat(sprintf("[short] %d of %d simulations reverse (%.1f%%); re-integrating on %d cores\n",
            nrow(S), nrow(M), 100 * nrow(S) / nrow(M), N_CORES))

trap = function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
one = function(k) {
  u = S[k]; W = min(u$tconv_off, u$tconv_onS)
  if (!is.finite(W) || W <= 0) return(NULL)
  trow = tg[.(u$batch, u$row_ind), nomatch = NULL]
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = u$Kg; p$r_L_P = u$r_L_P; p$r_LP_P = u$r_L_P
  TT = c(0, exp(seq(log(min(0.05, W / 1000)), log(W), length.out = TT_N)))
  ## OFF at beta_r, ON at the damage-matched beta_r*
  arm = function(th, br) {
    o = tryCatch(ode_ext(seed(u$Eh0), TT, rhs_app,
                         c(p, list(theta_treg = th, beta_ros = br)), maxsteps = 1e5),
                 error = function(e) NULL)
    if (is.null(o) || nrow(o) < length(TT)) return(NA_real_)
    d = as.data.table(o); trap(d$time, d$PLP)
  }
  data.table(i = u$i, win_short = W,
             P_short_off = arm(0, u$b_off), P_short_onS = arm(1, u$b_star))
}
t0 = Sys.time()
X = rbindlist(mclapply(seq_len(nrow(S)), one, mc.cores = N_CORES))
X = X[is.finite(P_short_off) & is.finite(P_short_onS) & P_short_off > 0]
X[, save_short := 1 - P_short_onS / P_short_off]
saveRDS(X, OUT)

## the Results statement, recomputed from what was just integrated
Y = merge(S[, .(i, tconv_off, tconv_onS, P_end_off, P_end_onS)], X, by = "i")
Y[, `:=`(resid_gap = P_end_onS - P_end_off, tconv_ratio = tconv_onS / tconv_off)]
cat(sprintf("[short] %d pairs integrated in %.1f min -> %s\n", nrow(X),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), OUT))
cat(sprintf("[short] lower with Tregs over the shorter window: %d of %d (%.1f%%), median saving %.1f%%\n",
            Y[save_short > 0, .N], nrow(Y), 100 * mean(Y$save_short > 0), 100 * median(Y$save_short)))
E = Y[save_short <= 0]
cat(sprintf("[short] exceptions: %d (%.1f%%); residual gap median %.2f vs %.2f in the rest; tconv ratio median %.2f vs %.2f\n",
            nrow(E), 100 * nrow(E) / nrow(Y), median(E$resid_gap),
            median(Y[save_short > 0]$resid_gap), median(E$tconv_ratio),
            median(Y[save_short > 0]$tconv_ratio)))
cat(sprintf("[short] median tconv ratio over all reversal pairs: %.1f\n", median(Y$tconv_ratio)))
