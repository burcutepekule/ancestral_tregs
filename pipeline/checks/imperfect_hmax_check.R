## pipeline/checks/imperfect_hmax_check.R -- open item 10: does hmax = 1 (jump times not given as tcrit)
## change the outcome? 3 units with a split cell at eta = 0.5, 2 noise realisations, TMAX = 3000,
## hmax 1 vs 0.1. Result 2026-10-01: max relative difference of tail L 9.1e-10, same outcome. ~4 min.
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages({library(data.table); library(deSolve)})
BRANCH = "pathogenic_help"
suppressMessages(source("./lib/skeleton_defs.R")); suppressMessages(source("./lib/ODE_SYSTEM_v6.R"))
U = as.data.table(readRDS("./data/fig05_units_pathogenic_help.rds"))
PU = as.data.table(readRDS("./data/fig05_switch_prob_units.rds"))
PC = as.data.table(readRDS("./data/fig05_switch_prob_cells.rds"))
## pick units with split cells at eta 0.5 (most sensitive)
sp = PC[abs(eta-0.5)<1e-9 & p_open>0 & p_open<1 & unit %in% U$unit][order(unit)]
set.seed(7); pick = sp[sample(.N, 3)]
print(pick)
TMAX = 3000; TT = seq(0, TMAX, 1)
res = list()
for (i in seq_len(nrow(pick))) {
  cc = U[unit == pick$unit[i]]; trow = tg[.(cc$batch, cc$row_ind)]
  p0 = c(STRUCTURAL, EFFERO); for (nm in KIN) p0[[nm]] = trow[[nm]][1]
  p0$Kg = cc$Kg; p0$beta_ros = cc$beta_ros; p0$r_L_P = cc$r_L_P; p0$r_LP_P = cc$r_L_P
  p = c(p0, list(theta_treg = 1, phi_instability = pick$phi[i], treg_eff = pick$eta[i]))
  for (s in 1:2) {
    set.seed(1000 + s); nf = make_noise_fun(TMAX, pick$eta[i])
    for (h in c(1, 0.1)) {
      t0 = Sys.time()
      o = ode_ext(c(seed(cc$Eh0), CumC = 0), TT, make_ode_system(nf), p, hmax = h)
      tl = o[o[, 1] > TMAX / 2, ]
      L = mean(tl[, "Ea"] + tl[, "En"])
      res[[length(res) + 1]] = data.table(unit = cc$unit, phi = pick$phi[i], eta = pick$eta[i], noise = s, hmax = h,
        L_tail = L, Kg = cc$Kg, Eh_tail = mean(tl[, "Eh"]), secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
      print(res[[length(res)]])
    }
  }
}
R = rbindlist(res)
print(dcast(R, unit + noise ~ hmax, value.var = "L_tail"))
W = dcast(R, unit + noise ~ hmax, value.var = "L_tail"); setnames(W, c("unit", "noise", "h01", "h1"))
cat(sprintf("max relative difference %.2g\n", max(abs(W$h01 - W$h1) / W$h1)))
