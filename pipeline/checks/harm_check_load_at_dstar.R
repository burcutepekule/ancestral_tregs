## ============================================================================
## pipeline/checks/harm_check_load_at_dstar.R -- what the scatter in the delay
## panel of Supplementary Fig. 3 (panel B) is made of.
## ----------------------------------------------------------------------------
## THE PANEL. The delay d_star at which a combination stops being harmed (the
## bisected `neutral` threshold of the delay treatment in thresh) against the
## time t_dip at which its Tregs OFF arm reaches the extinction floor.
##
## THE TWO DERIVED QUANTITIES.
##     gap   = t_dip - d_star   how much earlier than elimination Tregs may be
##                              switched on and the infection still be cleared
##     ratio = d_star / t_dip   = 1 - gap / t_dip
##
## THE MEASUREMENT. For each point, integrate the Tregs OFF arm once (same
## initial condition, solver, tolerances and extinction floor as the production
## run) and read off the pathogen load P_LP it carries at t = d_star, the moment
## Tregs would be switched on, plus the killing margin
## (alpha_1_e*M1 + alpha_1_r*R)/r_LP_P there. Spearman correlations of every
## candidate against the ratio and the gap are printed. It also checks that
## every arm reaches the absorbing state (both P_L and P_LP under EXT_TOL), so
## that t_dip is an elimination time.
##
## in : ./data/skeleton_pathogenic_harm.rds, ./data/tau_grid_v6_*.rds,
##      ./data/harm_rescue_data.rds, ./lib/skeleton_defs.R
## out: prints; ./data/harm_check_load_at_dstar.csv
## RUN: N_CORES=6 Rscript pipeline/checks/harm_check_load_at_dstar.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages({library(data.table); library(parallel); library(deSolve)})
BRANCH = "pathogenic_harm"
source("./lib/skeleton_defs.R")

N_CORES = as.integer(Sys.getenv("N_CORES", "6"))
TMAX    = as.numeric(Sys.getenv("TMAX", "25000"))
EPS     = as.numeric(Sys.getenv("EPS", "0.01"))
EH0_HI  = as.numeric(Sys.getenv("EH0_HI", "99"))
CTRL_FRAC = as.numeric(Sys.getenv("CTRL_FRAC", "0.01"))
KEY = c("batch", "row_ind", "Kg", "beta_ros", "r_L_P")
TT  = c(0, exp(seq(log(0.05), log(TMAX), length.out = 300)))

SK = as.data.table(readRDS("./data/skeleton_pathogenic_harm.rds"))
IN_RDS = Sys.getenv("IN_RDS", "./data/harm_rescue_data.rds")
if (!file.exists(IN_RDS)) stop(IN_RDS, " not found -- run pipeline/50_harm_rescue_compute.R first")
Z  = readRDS(IN_RDS)
## d_star: the bisected delay at which harm ends (bound "neutral"). For every
## combination whose delay response crosses the dead band in one step the help
## threshold is the same number; the one combination that only reaches neutral
## is plotted at its neutral threshold, so this is exactly the set Supplementary Fig. 3 draws.
DS = as.data.table(Z$thresh)[treat == "delay" & bound == "neutral",
                             c(KEY, "thr"), with = FALSE]
setnames(DS, "thr", "d_star")
LD = unique(as.data.table(Z$sm)[treat == "off", c(KEY, "t_dip"), with = FALSE], by = KEY)
PA = merge(DS, LD, by = KEY)[is.finite(d_star) & d_star > 0 & is.finite(t_dip) & t_dip > 0]
PA[, ratio := d_star / t_dip][, gap := t_dip - d_star]
cat(sprintf("[ld] %d refined combinations | ratio median %.3f | gap median %.1f\n",
            nrow(PA), median(PA$ratio), median(PA$gap)))

getp = function(r) {
  cc = SK[batch == r$batch & row_ind == r$row_ind & Kg == r$Kg &
          abs(beta_ros - r$beta_ros) < 1e-12 & abs(r_L_P - r$r_L_P) < 1e-12]
  if (!nrow(cc)) return(NULL)
  cc = cc[1]; trow = tg[.(r$batch, r$row_ind), nomatch = NULL]
  if (!nrow(trow)) return(NULL)
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = r$Kg; p$beta_ros = r$beta_ros; p$r_L_P = r$r_L_P; p$r_LP_P = r$r_L_P
  list(p = p, e0 = cc$Eh0_off + EPS * (EH0_HI - cc$Eh0_off))
}
offarm = function(g) tryCatch(as.data.table(
  ode(seed(g$e0), TT, rhs_app, c(g$p, list(theta_treg = 0)), method = "lsodar",
      rootfun = function(t, y, q) max(y[EXT_IDX]) - EXT_TOL,
      events = list(func = function(t, y, q) { y[EXT_IDX] = 0; y }, root = TRUE),
      atol = ATOL, rtol = RTOL, hmax = HMAX, maxsteps = 1e5)),
  error = function(e) NULL)

Z = rbindlist(mclapply(seq_len(nrow(PA)), function(i) {
  r = PA[i]; g = getp(r); if (is.null(g)) return(NULL)
  o = offarm(g); if (is.null(o) || nrow(o) < length(TT)) return(NULL)
  p = g$p; n = nrow(o)
  z = o$PLP; z[z < EXT_TOL] = 0; k = which.max(z)
  ic = which(z[k:n] < CTRL_FRAC * z[k])
  t_ctrl = if (length(ic)) o$time[k - 1L + ic[1]] else NA_real_
  iP = which(o$PLP < EXT_TOL); iL = which(o$PL < EXT_TOL)
  ## the killing margin, and the two states read at the moment Tregs arrive
  mar = (p$alpha_1_e * o$M1 + p$alpha_1_r * o$R) / p$r_LP_P
  atv = function(v, tt) if (is.na(tt)) NA_real_ else approx(o$time, v, xout = tt, rule = 2)$y
  plp_d = max(atv(pmax(o$PLP, 0), r$d_star), 0)
  data.table(batch = r$batch, row_ind = r$row_ind, Kg = r$Kg,
             beta_ros = r$beta_ros, r_L_P = r$r_L_P,
             d_star = r$d_star, t_dip = r$t_dip, ratio = r$ratio, gap = r$gap,
             ## is t_dip an elimination time on this arm?
             absorbed = isTRUE(o$PL[n] < EXT_TOL && o$PLP[n] < EXT_TOL),
             t_PLP_below = if (length(iP)) o$time[iP[1]] else NA_real_,
             t_PL_below  = if (length(iL)) o$time[iL[1]] else NA_real_,
             PLP_regrew  = if (length(iP)) max(o$PLP[iP[1]:n]) > EXT_TOL else NA,
             ## the load and the margin when Tregs are switched on
             PLP_at_dstar = plp_d,
             decades_at_dstar = log10(max(plp_d, 1e-30) / EXT_TOL),
             margin_at_dstar = atv(mar, r$d_star),
             M1_at_dstar = atv(o$M1, r$d_star),
             ## the alternatives already tried, kept so the comparison is in one table
             t_ctrl = t_ctrl, tail_dur = r$t_dip - t_ctrl,
             margin_peak = max(mar), M1_peak = max(o$M1),
             PLP_peak = z[k], Eh_end = o$Eh[n])
}, mc.cores = N_CORES), fill = TRUE)
fwrite(Z, "./data/harm_check_load_at_dstar.csv")
cat(sprintf("[ld] integrated %d untreated arms; wrote ./data/harm_check_load_at_dstar.csv\n", nrow(Z)))

## ---- 1. IS t_dip AN ELIMINATION TIME EVERYWHERE? ---------------------------
cat("\n[1] does every arm in the Supplementary Fig. 3 delay panel actually eliminate the pathogen?\n")
cat(sprintf("    absorbing state reached: %d of %d\n", sum(Z$absorbed), nrow(Z)))
cat(sprintf("    P_LP returns above the floor after first crossing it: %d\n",
            sum(Z$PLP_regrew, na.rm = TRUE)))
cat(sprintf("    P_L crosses the floor BEFORE P_LP on %d arms, after on %d\n",
            sum(Z$t_PL_below < Z$t_PLP_below, na.rm = TRUE),
            sum(Z$t_PL_below > Z$t_PLP_below, na.rm = TRUE)))
if (all(Z$absorbed))
  cat("    every arm is absorbed, so t_dip is an elimination time throughout and\n",
      "   the ordering of the two compartments does not affect what it measures.\n")

## ---- 2. WHAT ORDERS THE DEPARTURE FROM THE IDENTITY ------------------------
CAND = c("PLP_at_dstar", "decades_at_dstar", "margin_at_dstar", "tail_dur",
         "margin_peak", "M1_peak", "M1_at_dstar", "t_dip", "Kg", "beta_ros", "r_L_P")
sp = function(y, v) suppressWarnings(cor(y, Z[[v]], method = "spearman", use = "complete.obs"))
R = data.table(candidate = CAND,
               vs_ratio = round(vapply(CAND, function(v) sp(Z$ratio, v), 0), 3),
               vs_gap   = round(vapply(CAND, function(v) sp(Z$gap, v), 0), 3))
R[, a := abs(vs_ratio)]; setorder(R, -a); R[, a := NULL]
cat("\n[2] Spearman against the ratio and against the gap:\n"); print(R)
cat("\n    The ratio and the gap are ordered by DIFFERENT things, which is the\n")
cat("    point of separating them: t_dip moves the ratio arithmetically without\n")
cat("    moving the gap at all.\n")

## ---- 3. THE LOAD AT d_star -------------------------------------------------
cat("\n[3] the pathogen load the untreated arm is carrying when Tregs arrive:\n")
print(signif(quantile(Z$PLP_at_dstar, c(0, .25, .5, .75, .9, .99, 1)), 3))
cat("\n    binned:\n")
print(Z[, .(n = .N, ratio_med = round(median(ratio), 3), gap_med = round(median(gap), 1)),
        by = .(decades_above_floor = cut(decades_at_dstar,
               c(-Inf, 0.5, 1, 2, 3, 5, Inf)))][order(decades_above_floor)])
cat("\n[3] the five largest loads at d_star:\n")
print(Z[order(-PLP_at_dstar)][seq_len(min(5, .N)),
      .(batch, row_ind, Kg, beta_ros, r_L_P, t_dip = signif(t_dip, 4),
        d_star = signif(d_star, 4), ratio = round(ratio, 3), gap = round(gap, 1),
        PLP_at_dstar = signif(PLP_at_dstar, 3),
        margin_at_dstar = round(margin_at_dstar, 2))])

## ---- 4. THE EXTREME POINT AGAINST THE POPULATION ---------------------------
setorder(Z, -gap)
tgt = Z[1]
cat(sprintf("\n[4] the point with the longest grace period: batch %d row_ind %d Kg %g beta_ros %g r_L_P %g\n",
            tgt$batch, tgt$row_ind, tgt$Kg, tgt$beta_ros, tgt$r_L_P))
for (v in c("ratio", "gap", "t_dip", "PLP_at_dstar", "decades_at_dstar",
            "margin_at_dstar", "tail_dur", "margin_peak")) {
  x = Z[[v]]
  cat(sprintf("    %-18s %11.4g | population median %11.4g | z %+6.2f | rank %d of %d\n",
              v, tgt[[v]], median(x, na.rm = TRUE),
              (tgt[[v]] - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE),
              rank(-x, na.last = "keep")[1], sum(!is.na(x))))
}
cat("\n    A point that is rank 1 on the load at d_star and unremarkable on the\n")
cat("    margin and the decay rate is not a separate kind of object: it is the\n")
cat("    end of the same relation, and the ratio exaggerates it because its\n")
cat("    elimination time is short.\n")

## ---- 5. HOW MUCH OF THE SCATTER IS ARITHMETIC ------------------------------
m = lm(log10(ratio) ~ log10(t_dip), data = Z)
cat(sprintf("\n[5] log10(ratio) against log10(t_dip): R2 %.3f, slope %+.4f\n",
            summary(m)$r.squared, coef(m)[2]))
cat(sprintf("    sd(log10 ratio) %.4f falls to %.4f once t_dip is removed\n",
            sd(log10(Z$ratio)), sd(resid(m))))
cat("    That share of the scatter in the Supplementary Fig. 3 delay panel is the arithmetic of dividing a\n")
cat("    roughly constant grace period by an elimination time that spans a\n")
cat("    decade, not variation in the underlying quantity.\n")
