## ============================================================================
## pipeline/42_first_peak.R -- the FIRST PEAK of the pathogen in each arm of the
## damage-matched experiment, for every pair of the figure's set.
## ----------------------------------------------------------------------------
## WHY THIS EXISTS. 40_damage_matched.R stores peak_P as the maximum of P_LP over the
## whole run. For a pair that settles on a high-burden attractor that maximum is not a
## peak at all: the load climbs to the attractor and stays, so the "peak" reported is the
## steady state. Comparing it with the other arm's genuine peak compares two different
## things, and it is what put 78 pairs in the "more at the peak" stratum of Extended Data Fig. 5.
## Here each arm is re-integrated on a fine grid and the FIRST strict local maximum of
## P_LP is recorded -- the first time the load turns over, whatever happens later.
## NO PROMINENCE RULE IS APPLIED, and none is needed: over the 25,568 series of the set,
## the plain first local maximum and one required to fall 5% below the peak before any
## higher value pick the same point in 12,781 of 12,784 pairs (both arms). The three
## exceptions are shoulders on a climb to a higher attractor, and the comparison they
## enter is unchanged either way. A series that never turns over has no first peak (NA);
## one whose first peak sits within 1% of its final load had no transient peak, which
## `no_transient_*` flags rather than hides.
##
## in : ./data/fig04_damage_matched_gated.rds, ./data/fig04_gated_units.rds,
##      ./data/tau_grid_v6_*.rds, ./lib/*
## out: ./data/fig06_first_peak.rds  one row per pair: i, peak1_*, t_peak1_*, gmax_*, end_*
## RUN: N_CORES=8 Rscript pipeline/42_first_peak.R   (~20 min for 12,784 pairs)
## ============================================================================
suppressMessages({library(data.table); library(parallel); library(deSolve)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
assign("BRANCH", Sys.getenv("BRANCH", "pathogenic_harm"), envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R")); suppressMessages(source("./lib/extinction_floor.R"))
source("./lib/fig_damage_matched.R")
N_CORES   = as.integer(Sys.getenv("N_CORES", "8"))
TMAX      = as.numeric(Sys.getenv("TMAX", "12000"))
TT_TRAJ_N = as.integer(Sys.getenv("TT_TRAJ_N", "2000"))
OUT       = Sys.getenv("OUT", "./data/fig06_first_peak.rds")
TT = c(0, exp(seq(log(0.05), log(TMAX), length.out = TT_TRAJ_N)))

first_peak = function(t, y) {
  k = which(diff(sign(diff(y))) < 0) + 1L
  if (!length(k)) return(c(NA_real_, NA_real_))
  c(y[k[1]], t[k[1]])
}
S = fig_dm_matched()
arm = function(p, e0, beta, th) {
  pp = p; pp$beta_ros = beta
  o = tryCatch(ode_ext(seed(e0), TT, rhs_app, c(pp, list(theta_treg = th)), maxsteps = 1e5), error = function(e) NULL)
  if (is.null(o) || nrow(o) < length(TT)) return(NULL)
  as.data.table(o)[time > 0, .(time, PLP)]
}
one = function(k) {
  u = S[k]; trow = tg[.(u$batch, u$row_ind), nomatch = NULL]
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = u$Kg; p$r_L_P = u$r_L_P; p$r_LP_P = u$r_L_P
  a = arm(p, u$Eh0, u$b_off, 0); b = arm(p, u$Eh0, u$b_star, 1)       # OFF at beta_r, ON at beta*
  if (is.null(a) || is.null(b)) return(NULL)
  fa = first_peak(a$time, a$PLP); fb = first_peak(b$time, b$PLP)
  data.table(i = u$i, peak1_off = fa[1], t_peak1_off = fa[2], gmax_off = max(a$PLP), end_off = tail(a$PLP, 1),
             peak1_onS = fb[1], t_peak1_onS = fb[2], gmax_onS = max(b$PLP), end_onS = tail(b$PLP, 1))
}
cat(sprintf("[first peak] %d pairs x 2 arms on %d cores\n", nrow(S), N_CORES))
t0 = Sys.time()
R = rbindlist(mclapply(seq_len(nrow(S)), one, mc.cores = N_CORES))
R[, `:=`(no_transient_off = is.na(peak1_off) | peak1_off <= 1.01 * end_off,
         no_transient_onS = is.na(peak1_onS) | peak1_onS <= 1.01 * end_onS)]
saveRDS(R, OUT)
cat(sprintf("[first peak] %d pairs in %.1f min -> %s\n", nrow(R), as.numeric(difftime(Sys.time(), t0, units = "mins")), OUT))
cat(sprintf("[first peak] no first peak: OFF %d, ON* %d | no transient peak: OFF %d, ON* %d\n",
            R[is.na(peak1_off), .N], R[is.na(peak1_onS), .N], R[no_transient_off == TRUE, .N], R[no_transient_onS == TRUE, .N]))
B = R[!is.na(peak1_off) & !is.na(peak1_onS)]
cat(sprintf("[first peak] lower with Tregs: %d of %d (%.2f%%) | median reduction %.1f%%\n",
            B[peak1_onS < peak1_off, .N], nrow(B), 100 * B[, mean(peak1_onS < peak1_off)],
            100 * median(1 - B$peak1_onS / B$peak1_off)))
