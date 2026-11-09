## ============================================================================
## pipeline/44_window_damage.R -- the cumulative epithelial damage of the Tregs-OFF and
## the brake-alone (Tregs-ON at beta_r) arm over THEIR OWN COMMON WINDOW, for every
## robust pair of the damage-matched experiment.
## ----------------------------------------------------------------------------
## WHY THIS EXISTS. The experiment matches cost over a common window: the longer of the
## two convergence times of the pair being compared. The selection criterion -- the brake
## alone lowers cumulative damage by at least DAM_MIN -- was read off D_total_off and
## D_total_onB, the integrals over the whole run (TMAX), so criterion and match used
## different integrals. Here the same two arms are re-integrated and the epithelial
## deficit K_E - E_h is accumulated to win_B = max(tconv_off, tconv_onB), the window the
## brake-alone comparison is entitled to, so the criterion can be applied in the currency
## the experiment matches.
## The window is the pair's own: no arm is integrated beyond the time both have settled.
##
## in : ./data/fig04_damage_matched_gated.rds, ./data/fig04_gated_units.rds,
##      ./data/tau_grid_v6_*.rds, ./lib/*
## out: ./data/fig06_window_damage.rds  one row per robust pair: i, win_B, D_win_off,
##      D_win_onB, save_win (= 1 - D_win_onB / D_win_off)
## RUN: N_CORES=9 Rscript pipeline/44_window_damage.R   (~13 min for 48,842 pairs)
## ============================================================================
suppressMessages({library(data.table); library(deSolve); library(parallel)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
assign("BRANCH", Sys.getenv("BRANCH", "pathogenic_harm"), envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R")); suppressMessages(source("./lib/extinction_floor.R"))
N_CORES = as.integer(Sys.getenv("N_CORES", "9"))
TT_N    = as.integer(Sys.getenv("TT_N", "600"))
KE      = as.numeric(Sys.getenv("K_E", "100"))
OUT     = Sys.getenv("OUT", "./data/fig06_window_damage.rds")

## the robust pairs: steps 1 and 2 of lib/fig_damage_matched.R, before the damage floor
R = as.data.table(readRDS("./data/fig04_damage_matched_gated.rds"))[match == "conv"]
U = as.data.table(readRDS("./data/fig04_gated_units.rds"))[, .(i = .I, batch, row_ind)]
R = merge(R[, setdiff(names(R), c("batch", "row_ind")), with = FALSE], U, by = "i")
S = R[grepl("^matched", status) & is.finite(auc_P_onS) & is.finite(P_win_onS) & abs(match_err) <= 1e-4 &
      is.finite(jump_pts) & jump_pts <= 10 & n_cross_up == 1]
cat(sprintf("[window] %d robust pairs x 2 arms on %d cores\n", nrow(S), N_CORES))

trap = function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
one = function(k) {
  u = S[k]; W = max(u$tconv_off, u$tconv_onB)
  if (!is.finite(W) || W <= 0) return(NULL)
  trow = tg[.(u$batch, u$row_ind), nomatch = NULL]
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = u$Kg; p$r_L_P = u$r_L_P; p$r_LP_P = u$r_L_P; p$beta_ros = u$b_off
  TT = c(0, exp(seq(log(min(0.05, W / 1000)), log(W), length.out = TT_N)))
  arm = function(th) {
    o = tryCatch(ode_ext(seed(u$Eh0), TT, rhs_app, c(p, list(theta_treg = th)), maxsteps = 1e5),
                 error = function(e) NULL)
    if (is.null(o) || nrow(o) < length(TT)) return(NA_real_)
    d = as.data.table(o); trap(d$time, KE - d$Eh)
  }
  data.table(i = u$i, win_B = W, D_win_off = arm(0), D_win_onB = arm(1))
}
t0 = Sys.time()
X = rbindlist(mclapply(seq_len(nrow(S)), one, mc.cores = N_CORES))
X = X[is.finite(D_win_off) & is.finite(D_win_onB) & D_win_off > 0]
X[, save_win := 1 - D_win_onB / D_win_off]
saveRDS(X, OUT)
cat(sprintf("[window] %d pairs integrated in %.1f min -> %s\n", nrow(X),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), OUT))
for (f in c(0.02, 0.05, 0.10)) cat(sprintf("[window] saving >= %2.0f%% over the window: %d pairs (%.1f%% of the robust set)\n",
                                           100 * f, X[save_win >= f, .N], 100 * mean(X$save_win >= f)))
