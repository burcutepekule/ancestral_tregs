## ============================================================================
## pipeline/20_load_balance.R -- the load-balance trajectories behind Fig. 3C and
## Extended Data Fig. 2C (bistable class) and Extended Data Fig. 3 (ON_mono class).
##
## For each sampled combination, one injury that the Tregs-OFF tissue does not
## survive and the Tregs-ON tissue does -- the midpoint of the rescued interval on
## the initial-condition line -- is integrated under Tregs OFF (theta_treg = 0)
## and Tregs ON (theta_treg = 1) to the sweep horizon. With L = Ea + En,
##     dL/dt = Eh (tau_r R + tau_P_B PLP / K_LP_P) - g(L) (M0 + M1) L
##           = killing - clearance            (exact; M0 + M1 is conserved)
## and the overwhelmed-efferocytosis regime is L > Kg.
##
## POPULATION (gate inputs from the core pipeline, copied into ./data/):
##   rescue_<BRANCH>.rds, skeleton_<BRANCH>.rds, selected_big.rds
##   1. every combination of rescue_<BRANCH>.rds, joined to its (batch, row_ind);
##   2. EXCLUDE_DEV (default on): worst per-probe deviation from the assigned fixed
##      point <= MAP_DEV_TOL (1e-2) on BOTH rays (max over the probes, not the median);
##   3. bistable: both threshold directions "increasing", gap not censored,
##      rescue_gap >= GAP_MIN (5 sterile, 2 pathogenic);
##      ON_mono (rescue_mode "mono_B"): Tregs-OFF direction increasing, Eh0_off finite;
##   4. N_PER drawn at random per class (seed SEED), all of a class if fewer.
## Initial condition: seed(e_mid), e_mid = (Eh0_on + Eh0_off)/2 (bistable) or
## Eh0_off/2 (ON_mono). Output every 0.5 time units to TMAX = 12000.
##
## PER-TRAJECTORY SUMMARIES are computed on the FULL 0.5 output grid (not on the
## thinned `pop` series), including the endpoint quantities of the load-balance
## table: minimum of L and its time, load removed before it, M1 and L at the time
## of the Tregs-OFF minimum, minimum of dL/dt, peak M1 while L > Kg, final and
## peak pathogen load.
##
## in : ./data/rescue_<BRANCH>.rds, ./data/skeleton_<BRANCH>.rds, ./data/selected_big.rds,
##      ./data/tau_grid_v6_<n>.rds
## out: ./data/load_balance_<BRANCH><OUT_TAG>.rds
## RUN: BRANCH=sterile         N_PER=500 N_CORES=6 Rscript pipeline/20_load_balance.R
##      BRANCH=pathogenic_help N_PER=500 N_CORES=6 Rscript pipeline/20_load_balance.R
##      GATE_ONLY=1 BRANCH=sterile Rscript pipeline/20_load_balance.R   (counts only, no integration)
## ============================================================================
suppressMessages({library(data.table); library(deSolve); library(parallel)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
BRANCH = Sys.getenv("BRANCH", "sterile"); N_PER = as.integer(Sys.getenv("N_PER", "500"))
N_CORES = as.integer(Sys.getenv("N_CORES", "6")); SEED = as.integer(Sys.getenv("SEED", "1"))
TMAX = as.numeric(Sys.getenv("TMAX", "12000")); DT = as.numeric(Sys.getenv("DT", "0.5"))
GAP_MIN = as.numeric(Sys.getenv("GAP_MIN", if (BRANCH == "sterile") "5" else "2"))
OUT_TAG = Sys.getenv("OUT_TAG", ""); GATE_ONLY = identical(Sys.getenv("GATE_ONLY", "0"), "1")
## STORAGE OF THE DRAWN SERIES IS THINNED BY MOTION, NOT BY TIME: a point is kept when
## log10(L/Kg) has moved by more than MOVE_X or dL/dt by more than MOVE_Y since the last
## kept point, and in any case every STORE_BY time units. Summaries do not use it.
T_STORE = as.numeric(Sys.getenv("T_STORE", "12000")); STORE_BY = as.numeric(Sys.getenv("STORE_BY", "10"))
MOVE_X = as.numeric(Sys.getenv("MOVE_X", "0.004")); MOVE_Y = as.numeric(Sys.getenv("MOVE_Y", "0.02"))
keep_pt = function(tt, xk, yk) {
  lx = log10(pmax(xk, 1e-12)); n = length(tt); keep = logical(n); keep[1] = TRUE; lastx = lx[1]; lasty = yk[1]; lastt = tt[1]
  for (i in 2:n) { if (tt[i] > T_STORE) break
    if (abs(lx[i] - lastx) > MOVE_X || abs(yk[i] - lasty) > MOVE_Y || tt[i] - lastt >= STORE_BY) { keep[i] = TRUE; lastx = lx[i]; lasty = yk[i]; lastt = tt[i] } }
  keep[n] = keep[n] | tt[n] <= T_STORE; keep }
## The sweep matches each probe endpoint to the nearest fixed point with no distance limit
## and its trust flag gates on the MEDIAN probe deviation; a row with one probe far from
## every listed fixed point passes that gate. EXCLUDE_DEV drops rows whose WORST probe
## deviation on either ray exceeds MAP_DEV_TOL.
EXCLUDE_DEV = !identical(Sys.getenv("EXCLUDE_DEV", "1"), "0"); MAP_DEV_TOL = as.numeric(Sys.getenv("MAP_DEV_TOL", "1e-2"))
source("./lib/skeleton_defs.R"); CK = c("tau_id","Kg","beta_ros","r_L_P")
R = as.data.table(readRDS(sprintf("./data/rescue_%s.rds", BRANCH))); S = as.data.table(readRDS(sprintf("./data/skeleton_%s.rds", BRANCH)))
n_rescue = nrow(R)
R = merge(R, unique(S[, .(tau_id, Kg, beta_ros, r_L_P, batch, row_ind)]), by = CK)
cat(sprintf("[load_balance] %s | rescue file %d rows, %d joined to the skeleton | by mode: %s\n", BRANCH, n_rescue, nrow(R),
            paste(R[, .N, by = rescue_mode][, sprintf("%s %d", rescue_mode, N)], collapse = ", ")))
if (EXCLUDE_DEV) {
  Bg = as.data.table(readRDS("./data/selected_big.rds"))[, .(dev_max_on = max(dev_max_on, na.rm = TRUE), dev_max_off = max(dev_max_off, na.rm = TRUE)), by = CK]
  n0 = nrow(R); R = merge(R, Bg, by = CK, all.x = TRUE)
  ok = is.finite(R$dev_max_on) & R$dev_max_on <= MAP_DEV_TOL & is.finite(R$dev_max_off) & R$dev_max_off <= MAP_DEV_TOL
  cat(sprintf("[load_balance] EXCLUDE_DEV: %d of %d rows kept (worst probe deviation <= %g on both rays) | dropped by mode: %s\n",
              sum(ok), n0, MAP_DEV_TOL, paste(R[!ok, .N, by = rescue_mode][, sprintf("%s %d", rescue_mode, N)], collapse = ", ")))
  R = R[ok] }
set.seed(SEED)
BI = R[rescue_mode == "bistable" & rescue_gap >= GAP_MIN & dir_off == "increasing" & dir_on == "increasing" & !gap_censored]
MB = R[rescue_mode == "mono_B" & dir_off == "increasing" & is.finite(Eh0_off)]
cat(sprintf("[load_balance] gate: bistable %d (of %d after the deviation gate) | ON_mono %d (of %d)\n",
            nrow(BI), R[rescue_mode == "bistable", .N], nrow(MB), R[rescue_mode == "mono_B", .N]))
U = rbindlist(list(BI[sample(.N, min(N_PER, .N))], MB[sample(.N, min(N_PER, .N))])); U[, cid := .I]
U[, e_mid := ifelse(rescue_mode == "bistable", (Eh0_on + Eh0_off) / 2, Eh0_off / 2)]
cat(sprintf("[load_balance] %s | %d bistable + %d ON_mono drawn (seed %d) | injury = midpoint of the rescued interval | TMAX %g\n",
            BRANCH, sum(U$rescue_mode == "bistable"), sum(U$rescue_mode == "mono_B"), SEED, TMAX))
if (GATE_ONLY) { saveRDS(U, sprintf("./data/load_balance_sample_%s%s.rds", BRANCH, OUT_TAG)); quit(save = "no") }

integrate_one = function(cc, p0, th) {
  p = p0; p$theta_treg = th
  o = tryCatch(as.data.table(ode_ext(seed(cc$e_mid), seq(0, TMAX, by = DT), rhs_app, p, maxsteps = 4e5)), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  o[, L := Ea + En]; o[, g := p$a_floor + (p$a0 - p$a_floor) * p$Kg^p$n_g / (p$Kg^p$n_g + L^p$n_g)]
  ## killing: ROS (tau_r R) and, on the pathogenic branches, the pathogen (tau_P_B PLP / K_LP_P)
  o[, kill := Eh * (p$tau_r * R + p$tau_P_B * PLP / p$K_LP_P)]; o[, clear := g * (M0 + M1) * L]; o[, dL := kill - clear]; o[, lowg := L > p$Kg]
  ## pathogen balance (identically zero on the sterile branch)
  o[, p_supply := PL * p$gamma_B_P * (1 - Eh / p$K_E) * (1 - PLP / p$K_LP_P) + PLP * p$r_LP_P * (1 - PLP / p$K_LP_P)]
  o[, p_removal := PLP * (p$alpha_0 * M0 + p$alpha_1_e * M1 + p$alpha_1_r * R)]
  list(o = o, p = p) }

one = function(i) {
  cc = U[i]; trow = tg[.(cc$batch, cc$row_ind), nomatch = NULL][1]
  p0 = c(STRUCTURAL, EFFERO); for (nm in KIN) p0[[nm]] = trow[[nm]]; p0$Kg = cc$Kg; p0$beta_ros = cc$beta_ros; p0$r_L_P = cc$r_L_P; p0$r_LP_P = cc$r_L_P
  RUNS = lapply(c(`0` = 0, `1` = 1), function(th) integrate_one(cc, p0, th))
  ## time of the Tregs-OFF global minimum of L; both conditions are read there
  t_off_min = if (!is.null(RUNS[["0"]])) { oo = RUNS[["0"]]$o; oo$time[which.min(oo$L)] } else NA_real_
  rbindlist(lapply(names(RUNS), function(k) {
    r = RUNS[[k]]; if (is.null(r)) return(NULL); o = r$o; p = r$p; th = as.numeric(k)
    k95 = which(o$g / p$a0 >= 0.95)[1]; low = o[lowg == TRUE]; jmin = which.min(o$L); jo = if (is.finite(t_off_min)) which.min(abs(o$time - t_off_min)) else NA_integer_
    data.table(cid = i, cls = cc$rescue_mode, batch = cc$batch, row_ind = cc$row_ind, Kg = cc$Kg, beta_ros = cc$beta_ros, r_L_P = cc$r_L_P, rescue_gap = cc$rescue_gap,
               Eh0_off = cc$Eh0_off, Eh0_on = cc$Eh0_on, e_mid = cc$e_mid, cond = ifelse(th == 1, "Tregs ON", "Tregs OFF"),
               L0 = o$L[1], starts_lowg = o$lowg[1], resolves = !is.na(k95), t95 = if (is.na(k95)) NA_real_ else o$time[k95],
               L_end = tail(o$L, 1), L_min = min(o$L), ends_lowg = tail(o$lowg, 1),
               ## at rest at the horizon: |dL/dt| < 1e-3 at every output time of the last 100 time units
               steady_end = all(abs(tail(o$dL, 201)) < 1e-3),
               dL_turns_pos_in_lowg = any(low$dL > 0), max_dL_lowg = if (nrow(low)) max(low$dL) else NA_real_,
               t_first_pos = if (any(low$dL > 0)) low$time[which(low$dL > 0)[1]] else NA_real_,
               ## exit = the load stays below Kg from here on
               t_exit_lowg = if (!tail(o$lowg, 1)) o$time[max(which(o$lowg)) + 1] else NA_real_,
               M1_max_lowg = if (nrow(low)) max(low$M1) else NA_real_, kill_over_clear_lowg = if (nrow(low)) sum(low$kill) / sum(low$clear) else NA_real_,
               ## ---- full-grid endpoint quantities (load-balance table) ----
               t_end = tail(o$time, 1), t_Lmin = o$time[jmin], frac_removed = (o$L[1] - o$L[jmin]) / o$L[1],
               L_ever_le_Kg = any(o$L <= p$Kg), g_end_over_a0 = tail(o$g, 1) / p$a0,
               t_off_min = t_off_min, L_at_off_min = if (is.na(jo)) NA_real_ else o$L[jo], M1_at_off_min = if (is.na(jo)) NA_real_ else o$M1[jo],
               dL_min = min(o$dL), M1_max = max(o$M1), M1_end = tail(o$M1, 1),
               PLP_end_over_K = tail(o$PLP, 1) / p$K_LP_P, PLP_max_over_K = max(o$PLP) / p$K_LP_P, t_PLP_max = o$time[which.max(o$PLP)],
               ## the drawn series, thinned by motion (see keep_pt)
               pop = list(o[keep_pt(time, L / p$Kg, dL), .(time, L_over_Kg = L / p$Kg, kill_over_clear = kill / pmax(clear, 1e-12), kill_minus_clear = dL, M1, Tr,
                                                                   PLP_over_K = PLP / p$K_LP_P, p_removal_over_supply = p_removal / pmax(p_supply, 1e-12),
                                                                   p_removal_minus_supply = p_removal - p_supply)])) })) }
t0 = Sys.time(); D = rbindlist(mclapply(seq_len(nrow(U)), one, mc.cores = N_CORES), fill = TRUE)
cat(sprintf("[load_balance] %d trajectories (of %d) in %.1f min\n", nrow(D), 2L * nrow(U), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
OUT = sprintf("./data/load_balance_%s%s.rds", BRANCH, OUT_TAG); saveRDS(D, OUT); cat("[load_balance] wrote", OUT, "\n")

## ---- console checks --------------------------------------------------------
cat("\n[load_balance] every injury starts in the low-g regime (L0 > Kg)?", all(D$starts_lowg), "\n")
print(D[, .(n = .N, resolves = sum(resolves), ends_in_lowg = sum(ends_lowg), steady_at_horizon = sum(steady_end), reached_TMAX = sum(t_end >= TMAX)), by = .(cls, cond)][order(cls, cond)])
print(D[, .(n = .N, dL_turns_positive = sum(dL_turns_pos_in_lowg), leaves_lowg = sum(!is.na(t_exit_lowg)),
            M1max_lowg_med = round(median(M1_max_lowg, na.rm = TRUE), 1)), by = .(cls, cond)][order(cls, cond)])
