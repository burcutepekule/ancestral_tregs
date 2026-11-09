## ============================================================================
## pipeline/48_gated_fourth_arm.R -- the fourth arm for every gain case of
## the gated sample: Tregs OFF at beta_r* (the same aggression without the brake).
## ----------------------------------------------------------------------------
## Gain case: Tregs OFF at beta_r persists, Tregs ON at beta_r* eliminates
## (data/fig04_damage_matched_gated.rds; damage matched within MATCH_TOL, beta_r* >=
## RATIO_MIN beta_r). For each, Tregs OFF at beta_r* is integrated and its damage over
## the pair's matching window compared with Tregs OFF at beta_r: cost_noTreg =
## D(OFF, beta_r*) / D(OFF, beta_r). Where it is close to 1 the extra aggression was
## nearly free without Tregs and the gain is not attributable to them; where it is
## well above 1, Tregs are what made the aggression affordable.
## Also recorded: whether OFF at beta_r* eliminates, its final E_h and M1, and the
## damage of the three stored arms re-integrated on the same grid.
## in : ./data/fig04_damage_matched_gated.rds, ./data/fig04_gated_units.rds, ./lib/*
## out: ./data/aux_gated_fourth_arm.rds
## RUN: N_CORES=6 Rscript pipeline/48_gated_fourth_arm.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
assign("BRANCH", Sys.getenv("BRANCH", "pathogenic_harm"), envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R"))
suppressMessages(source("./lib/extinction_floor.R"))
suppressMessages({library(data.table); library(parallel)})
MATCH_TOL = as.numeric(Sys.getenv("MATCH_TOL", "1e-4"))
RATIO_MIN = as.numeric(Sys.getenv("RATIO_MIN", "1.05"))
ELIM_TOL  = as.numeric(Sys.getenv("ELIM_TOL", "1e-3"))
N_CORES   = as.integer(Sys.getenv("N_CORES", "6"))
TMAX      = as.numeric(Sys.getenv("TMAX", "12000"))
KE        = as.numeric(Sys.getenv("K_E", "100"))
OUT       = "./data/aux_gated_fourth_arm.rds"
TT = c(0, exp(seq(log(0.05), log(TMAX), length.out = 400)))   # the experiment's grid
trap = function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

R = as.data.table(readRDS("./data/fig04_damage_matched_gated.rds"))[match == "conv"]
U = as.data.table(readRDS("./data/fig04_gated_units.rds"))[, .(i = .I, treg_effect)]
R = merge(R, U, by = "i")
el = function(tc, pe) is.finite(tc) | (is.finite(pe) & pe < ELIM_TOL)
R[, `:=`(el_off = el(t_clear_off, P_end_off), el_onS = el(t_clear_onS, P_end_onS))]
G = R[grepl("^matched", status) & abs(match_err) <= MATCH_TOL & b_star >= RATIO_MIN * b_off & !el_off & el_onS]
cat(sprintf("[fourth arm] %d gain cases\n", nrow(G)))

one = function(k) {
  u = G[k]; trow = tg[.(u$batch, u$row_ind), nomatch = NULL]
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = u$Kg; p$r_L_P = u$r_L_P; p$r_LP_P = u$r_L_P
  run = function(beta, th) { pp = p; pp$beta_ros = beta
    o = tryCatch(as.data.table(ode_ext(seed(u$Eh0), TT, rhs_app, c(pp, list(theta_treg = th)), maxsteps = 1e5)), error = function(e) NULL)
    if (is.null(o) || nrow(o) < length(TT)) NULL else o }
  sm = function(o, sfx) { if (is.null(o)) return(setNames(as.list(rep(NA_real_, 5)), paste0(c("D_win", "t_clear", "P_end", "Eh_end", "M1_end"), sfx)))
    k = o$time <= u$win
    setNames(list(trap(o$time[k], KE - o$Eh[k]), o$time[which(o$PLP <= EXT_TOL)[1]], o$PLP[nrow(o)], o$Eh[nrow(o)], o$M1[nrow(o)]),
             paste0(c("D_win", "t_clear", "P_end", "Eh_end", "M1_end"), sfx)) }
  as.data.table(c(list(i = u$i), sm(run(u$b_off, 0), "_OFF"), sm(run(u$b_off, 1), "_ON_same"), sm(run(u$b_star, 1), "_ON_star"), sm(run(u$b_star, 0), "_OFF_star")))
}
t0 = Sys.time()
X = rbindlist(mclapply(seq_len(nrow(G)), one, mc.cores = N_CORES), fill = TRUE)
cat(sprintf("[fourth arm] %d cases x 4 arms in %.1f min\n", nrow(X), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
X = merge(G[, .(i, tau_id, batch, row_ind, Kg, r_L_P, Eh0, b_off, b_star, ratio = b_star / b_off, win, t_div = t_div_M1, treg_effect,
                jump_pts, n_cross_up, P_end_off_exp = P_end_off)], X, by = "i")
X[, `:=`(cost_noTreg = D_win_OFF_star / D_win_OFF, match_check = D_win_ON_star / D_win_OFF,
         offstar_elim = is.finite(t_clear_OFF_star) | (is.finite(P_end_OFF_star) & P_end_OFF_star < ELIM_TOL))]
saveRDS(X, OUT)
q = function(x) paste(signif(quantile(x, c(.1, .25, .5, .75, .9), na.rm = TRUE), 3), collapse = " / ")
cat(sprintf("[fourth arm] match re-check D(ON*)/D(OFF) on this grid: %s\n", q(X$match_check)))
cat(sprintf("[fourth arm] cost of the same aggression without Tregs, D(OFF, beta*)/D(OFF, beta): %s\n", q(X$cost_noTreg)))
for (th in c(1.02, 1.05, 1.1, 1.2, 1.5)) cat(sprintf("   cost_noTreg > %.2f: %d (%.0f%%)\n", th, sum(X$cost_noTreg > th, na.rm = TRUE), 100 * mean(X$cost_noTreg > th, na.rm = TRUE)))
cat(sprintf("[fourth arm] OFF at beta* eliminates too: %d of %d; its final E_h minus OFF at beta_r: %s\n", sum(X$offstar_elim, na.rm = TRUE), nrow(X), q(X$Eh_end_OFF_star - X$Eh_end_OFF)))
cat("\nby the Treg effect at the pair's own beta_r:\n")
print(X[, .(cases = .N, cost_med = round(median(cost_noTreg, na.rm = TRUE), 3), cost_gt_1.1 = sum(cost_noTreg > 1.1, na.rm = TRUE), cost_gt_1.2 = sum(cost_noTreg > 1.2, na.rm = TRUE),
            offstar_elim = sum(offstar_elim, na.rm = TRUE), ratio_med = round(median(ratio), 2)), by = treg_effect][order(-cases)])
cat("\ncost of the same aggression without Tregs against the size of the beta_r increase:\n")
X[, rcl := cut(ratio, c(1.05, 1.5, 2, 3, 5, Inf), right = FALSE)]
print(X[, .(cases = .N, cost_med = round(median(cost_noTreg, na.rm = TRUE), 3), cost_gt_1.1 = sum(cost_noTreg > 1.1, na.rm = TRUE)), by = rcl][order(rcl)])
cat("\nthe 10 cases where Tregs matter most (largest cost without them), robust matches:\n")
print(X[jump_pts <= 10 & n_cross_up == 1][order(-cost_noTreg)][1:10, .(i, tau_id, Kg, r_L_P, Eh0, b_off, ratio = round(ratio, 2), cost_noTreg = round(cost_noTreg, 2),
      P_end_OFF = round(P_end_OFF, 1), P_end_ON_same = round(P_end_ON_same, 1), t_clear_ON_star = round(t_clear_ON_star), Eh_OFF = round(Eh_end_OFF, 1),
      Eh_ON_star = round(Eh_end_ON_star, 1), Eh_OFF_star = round(Eh_end_OFF_star, 1), offstar_elim, treg_effect)])
cat(sprintf("[fourth arm] wrote %s\n", OUT))
