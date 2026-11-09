## ============================================================================
## pipeline/25_corridor_persistence.R -- does the corridor of Supplementary Fig. 1A stay open
## in combinations other than the worked example?
##
## Supplementary Fig. 1A re-locates both basin boundaries along En at 60 times for ONE sterile
## combination. This script asks the same question of a random sample of the
## combinations the figure's panels C and D are drawn on, at three times only.
##
## PER COMBINATION. The t = 0 walls are already measured (x_off, x_on in the
## corridor file), so only the later times need work: the "below" seed -- the
## initial condition just inside the resolving side, the one panel A takes its
## slice along -- is integrated under each condition to t, and the En wall is
## located from the state it has reached, by the scan and bisection of
## lib/corridor_defs.R (find_wall). Four wall searches per combination.
##
## A wall that the scan cannot resolve (no label change in range, or more than
## one) is recorded as NA with a status, NOT as a closed corridor: the two are
## different answers and collapsing them would manufacture the result.
##
## in : ./data/corridor_q_<BRANCH><CORR_TAG>.rds, ./data/selected_big.rds, tau grids
## out: ./data/corridor_persist_<BRANCH>.csv
## RUN: BRANCH=sterile PERSIST_N=50 PERSIST_CORES=8 Rscript pipeline/25_corridor_persistence.R
## ============================================================================
suppressMessages({library(data.table); library(parallel)})
SRC = Sys.getenv("SRC", "."); setwd(SRC)
Sys.setenv(COMPOSITE = "1")
if (!nzchar(Sys.getenv("RAY_CLUST_TOL"))) Sys.setenv(RAY_CLUST_TOL = "1")
BR      = Sys.getenv("BRANCH", "sterile")
## N_PICK, not NSAMP: lib/corridor_defs.R defines a variable of that name (its own
## time-sample count, default 40) and is sourced below, so the assignment would be
## silently overwritten and every run would draw 40 combinations whatever was asked.
N_PICK  = as.integer(Sys.getenv("PERSIST_N", "500"))
NC      = as.integer(Sys.getenv("PERSIST_CORES", "6"))
TIMES   = as.numeric(strsplit(Sys.getenv("TIMES", "1000,10000"), ",")[[1]])
WIN     = Sys.getenv("WIN", "En")
SEED    = as.integer(Sys.getenv("SEED", "1"))
if (BR == "sterile" && !nzchar(Sys.getenv("CORR_TAG"))) Sys.setenv(CORR_TAG = "_round1")
assign("BRANCH", BR, envir = globalenv())
suppressMessages(source("./lib/skeleton_defs.R")); suppressMessages(source("./lib/extinction_floor.R"))
suppressMessages(source("./lib/corridor_defs.R"))

CK   = c("tau_id","Kg","beta_ros","r_L_P")
STEP = 10 / 2^7
EFF  = unique(as.data.table(readRDS("./data/selected_big.rds"))[, c(CK, "treg_effect"), with = FALSE])
L    = merge(load_corridor(BR), EFF, by = CK, all.x = TRUE)
L[, gain := 100 * (x_on - x_off) / (x_off - base)]
## the population of panels C and D, restricted to the coordinate panel A uses
P = L[treg_effect == "help" & coord == WIN & q_kind == "measured" & rescue_mode == "bistable" &
      is.finite(gain) & gain > 0 & is.finite(rescue_gap) & !gap_censored & rescue_gap > STEP * 1.01]
cat(sprintf("[persist] %s | %s combinations in the plotted population along %s\n",
            BR, format(nrow(P), big.mark = ","), WIN))
set.seed(SEED); S = P[sample(.N, min(N_PICK, .N))]
cat(sprintf("[persist] CONFIG asked for %d, population has %d, sampled %d | times %s | %d cores | NSCAN_D %d BIS_D %d\n",
            N_PICK, nrow(P), nrow(S), paste(TIMES, collapse = ", "), NC, NSCAN_D, BIS_D))

one = function(i) {
  pick = as.list(S[i])
  out = data.table(tau_id = pick$tau_id, Kg = pick$Kg, beta_ros = pick$beta_ros, r_L_P = pick$r_L_P,
                   rescue_gap = pick$rescue_gap, base = pick$base,
                   off_t0 = pick$x_off, on_t0 = pick$x_on, w_t0 = pick$x_on - pick$x_off)
  res = tryCatch({
    eng = combo_engine(pick); sds = make_seeds(pick)
    y0  = eng$SETC(eng$ybase, WIN, sds[["below"]])
    ST  = names(eng$ybase)
    for (tt in TIMES) {
      w = sapply(c(0L, 1L), function(th) {
        o = eng$run(y0, th, c(0, tt))
        if (is.null(o) || nrow(o) < 2) return(NA_real_)
        yt = setNames(as.numeric(o[nrow(o), ST]), ST)
        unname(find_wall(eng, WIN, yt, th)[["b"]]) })
      set(out, j = sprintf("off_t%g", tt), value = w[1])
      set(out, j = sprintf("on_t%g",  tt), value = w[2])
      set(out, j = sprintf("w_t%g",   tt), value = w[2] - w[1])
    }
    out[, status := "ok"]; out
  }, error = function(e) { out[, status := paste("error:", conditionMessage(e))]; out })
  cat(sprintf("[persist] %s %3d/%d tau_id %d%s\n", BR, i, nrow(S), pick$tau_id,
              if (!is.null(res$status)) paste0(" ", res$status[1]) else ""))
  res }

R = rbindlist(mclapply(seq_len(nrow(S)), one, mc.cores = NC), fill = TRUE)
OUT = if (WIN == "En") sprintf("./data/corridor_persist_%s.csv", BR) else sprintf("./data/corridor_persist_%s_%s.csv", BR, WIN)
fwrite(R, OUT)
cat(sprintf("\n[persist] wrote %s (%d rows)\n", OUT, nrow(R)))
for (tt in TIMES) {
  wc = R[[sprintf("w_t%g", tt)]]
  cat(sprintf("[persist] t = %6g | resolved %d/%d | open (w>0) %d | median width %.2f (t0 median %.2f)\n",
              tt, sum(is.finite(wc)), nrow(R), sum(wc > 0, na.rm = TRUE),
              median(wc, na.rm = TRUE), median(R$w_t0, na.rm = TRUE)))
}
