## ============================================================================
## pipeline/41_damage_matched_trajectories.R -- the trajectories behind Fig. 4 B and C.
## ----------------------------------------------------------------------------
## The experiment (40_damage_matched.R) stores per-unit summaries only. Panels B and C
## need time courses: they are re-integrated here for N_TRAJ pairs drawn at random
## (seed SEED) among the pairs of the figure (lib/fig_damage_matched.R: fig_dm_matched,
## with DAM_MIN as set -- the one population every panel of Fig. 4 uses),
## with exactly the experiment's code path -- same parameters, initial condition and
## solver -- on a finer output grid (TT_TRAJ_N log-spaced points instead of 400).
## Three runs per pair:
##   OFF       Tregs OFF at beta_r
##   ON_same   Tregs ON  at beta_r
##   ON_star   Tregs ON  at beta* (equal cumulative damage over the window)
## Stored states: time, P_LP, E_h, T_r, M1.
##
## in : ./data/fig04_damage_matched_gated.rds, ./data/fig04_gated_units.rds,
##      ./data/tau_grid_v6_*.rds, ./lib/*
## out: ./data/fig04b_v3_trajectories.rds  list(units = the drawn pairs, traj = time courses)
## RUN: N_TRAJ=1000 SEED=1 N_CORES=6 Rscript pipeline/41_damage_matched_trajectories.R   (the defaults; minutes)
## ============================================================================
suppressMessages({library(data.table); library(parallel); library(deSolve)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
assign("BRANCH", Sys.getenv("BRANCH", "pathogenic_harm"), envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R"))
suppressMessages(source("./lib/extinction_floor.R"))
source("./lib/fig_damage_matched.R")
N_TRAJ    = as.integer(Sys.getenv("N_TRAJ", "1000"))
N_CORES   = as.integer(Sys.getenv("N_CORES", "6"))
SEED      = as.integer(Sys.getenv("SEED", "1"))
TMAX      = as.numeric(Sys.getenv("TMAX", "12000"))
TT_TRAJ_N = as.integer(Sys.getenv("TT_TRAJ_N", "2000"))
OUT       = Sys.getenv("OUT", "./data/fig04b_v3_trajectories.rds")
TT_TRAJ   = c(0, exp(seq(log(0.05), log(TMAX), length.out = TT_TRAJ_N)))

M = fig_dm_matched()
run = function(p, e0, beta, th) {
  pp = p; pp$beta_ros = beta
  o = tryCatch(ode_ext(seed(e0), TT_TRAJ, rhs_app, c(pp, list(theta_treg = th)), maxsteps = 1e5),
               error = function(e) NULL)
  if (is.null(o) || nrow(o) < length(TT_TRAJ)) return(NULL)
  as.data.table(o)[, .(time, PLP, Eh, Tr, M1)]
}
traj_unit = function(k) {
  u = S[k]
  trow = tg[.(u$batch, u$row_ind), nomatch = NULL]
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = u$Kg; p$r_L_P = u$r_L_P; p$r_LP_P = u$r_L_P
  o1 = run(p, u$Eh0, u$b_off, 0); o2 = run(p, u$Eh0, u$b_off, 1); o3 = run(p, u$Eh0, u$b_star, 1)
  if (is.null(o1) || is.null(o2) || is.null(o3)) return(NULL)
  rbind(o1[, run := "OFF"], o2[, run := "ON_same"], o3[, run := "ON_star"])[, i := u$i][]
}
set.seed(SEED)
S = M[sample.int(.N, min(N_TRAJ, .N))]
cat(sprintf("[fig04 traj] re-integrating %d pairs x 3 runs on %d cores\n", nrow(S), N_CORES))
t0 = Sys.time()
TR = rbindlist(mclapply(seq_len(nrow(S)), traj_unit, mc.cores = N_CORES))
saveRDS(list(units = S, traj = TR), OUT)
cat(sprintf("[fig04 traj] %d units with all three runs in %.1f min -> %s\n", uniqueN(TR$i),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), OUT))
