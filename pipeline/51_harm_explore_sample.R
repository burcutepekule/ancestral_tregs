## ============================================================================
## pipeline/51_harm_explore_sample.R -- the help-against-harm trajectory sample
## that panel A of Supplementary Fig. 3 draws.
## ----------------------------------------------------------------------------
## ONE PANEL'S DATA. The median P_LP trajectory with its interquartile band, for
## four groups: combinations Tregs help and combinations Tregs harm, each under
## Tregs OFF and under Tregs ON.
##
## WHY IT CANNOT COME FROM harm_rescue_data.rds. That file holds the HARMED
## combinations only. This panel is a contrast between the harmed and the
## helped, so half its content is not in there. It therefore has its own source,
## selected_big.rds (written by the selection stage), its own random sample and
## its own cache.
##
## THE INITIAL CONDITION IS MATCHED ACROSS THE TWO CONDITIONS. Both arms start
## just inside the Tregs OFF threshold, on the resolving side,
## e0 = Eh0_off + EPS*(EH0_HI - Eh0_off), so the only difference between them
## is theta_treg.
##
## THE DISPLAY FLOOR IS 1 UNIT, AND IT IS A DISPLAY CHOICE, NOT THE MODEL'S.
## The integration's own floor is EXT_TOL = 1e-11. The panel clamps P_LP at
## PLP_FLOOR = 1 when drawing (Supplementary Fig. 3 does the same with PLP_DISP), so clearance is
## not visible in that panel; it is not the panel's job.
##
## THE BANDS ARE QUANTILES ACROSS COMBINATIONS at each output time, not a
## confidence interval on the median.
##
## in : ./data/selected_big.rds, ./data/tau_grid_v6_<n>.rds, ./lib/skeleton_defs.R
## out: ./data/harm_explore_cache.rds   (the integrated sample, reused when valid)
##      optionally a preview PNG at OUT_FIG (default: none)
## RUN: N_PER=250 N_CORES=6 Rscript pipeline/51_harm_explore_sample.R
##      RECOMPUTE=1 N_PER=250 N_CORES=6 Rscript pipeline/51_harm_explore_sample.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages({library(data.table); library(deSolve); library(parallel)
                  library(ggplot2)})
BRANCH = "pathogenic_harm"
source("./lib/skeleton_defs.R")

N_PER   = as.integer(Sys.getenv("N_PER", "250"))
N_CORES = as.integer(Sys.getenv("N_CORES", "4"))
TMAX    = as.numeric(Sys.getenv("TMAX", "12000"))
EPS     = as.numeric(Sys.getenv("EPS", "0.01"))
EH0_HI  = as.numeric(Sys.getenv("EH0_HI", "99"))
SEED    = as.integer(Sys.getenv("SEED", "3"))
GAP_MIN = as.numeric(Sys.getenv("GAP_MIN", "1e-9"))
FLIP_KEEP = as.integer(Sys.getenv("FLIP_KEEP", "0")) == 1L
PLP_FLOOR = as.numeric(Sys.getenv("PLP_FLOOR", "1"))
BANDS   = Sys.getenv("BANDS", "both")
if (!BANDS %in% c("both", "on", "none")) stop("BANDS must be both, on or none")
CACHE   = Sys.getenv("CACHE", "./data/harm_explore_cache.rds")
RECOMP  = identical(Sys.getenv("RECOMPUTE", "0"), "1")
OUT_FIG = Sys.getenv("OUT_FIG", "")   # empty: write the cache only
FIG_W   = as.numeric(Sys.getenv("FIG_W", "2600"))
FIG_H   = as.numeric(Sys.getenv("FIG_H", "2100"))
FIG_RES = as.numeric(Sys.getenv("FIG_RES", "200"))
BASE_SZ = as.numeric(Sys.getenv("BASE_SIZE", "17"))

## WHAT MAKES A CACHE STALE: everything that changes which combinations are
## drawn or how they are integrated. The display floor is deliberately NOT in
## the key -- it is applied when drawing, so changing it must not force an hour
## of re-integration.
KEY = list(N_PER = N_PER, SEED = SEED, TMAX = TMAX, HMAX = HMAX, EPS = EPS,
           EH0_HI = EH0_HI, GAP_MIN = GAP_MIN, FLIP_KEEP = FLIP_KEEP,
           ATOL = ATOL, RTOL = RTOL, EXT_TOL = EXT_TOL, schema = 1L)

## ---- the gate, the same one every population figure applies ----------------
## FLIP_KEEP removes dir_on == "decreasing": the Tregs ON resolving side is
## inverted there, so no boundary shift can be read and rescue_gap is undefined.
## ON_mono carries dir_on = NA and is not affected, since the filter is on
## "decreasing" rather than on "not increasing". GAP_MIN then applies to the
## bistable combinations only, where the gap is a distance between two
## crossings; for ON_mono the recorded gap is the whole Eh0 range and filtering
## on it would mean something else.
gate = function(D) {
  n0 = nrow(D)
  if (!FLIP_KEEP) D = D[is.na(dir_on) | dir_on != "decreasing"]
  n1 = nrow(D)
  D = D[rescue_mode != "bistable" | (is.finite(rescue_gap) & rescue_gap >= GAP_MIN)]
  cat(sprintf("[ex] gate: %s in, %s after dropping dir_on = decreasing, %s after GAP_MIN = %g on the bistable ones\n",
              format(n0, big.mark = ","), format(n1, big.mark = ","),
              format(nrow(D), big.mark = ","), GAP_MIN))
  D
}

compute = function() {
  S = as.data.table(readRDS("./data/selected_big.rds"))
  P = S[patho == TRUE & treg_effect %in% c("harm", "help") &
        is.finite(Eh0_off) & Eh0_off < 95]
  if (!nrow(P)) stop("no pathogenic help/harm combinations in selected_big.rds")
  P = gate(P)
  if (!nrow(P)) stop("the gate removed every combination; lower GAP_MIN")
  set.seed(SEED)
  P = P[, .SD[sample(.N, min(N_PER, .N))], by = treg_effect]
  cat("[ex] sampled:\n"); print(P[, .N, by = treg_effect])

  ## A common output grid, dense early because M1 peaks and the pathogen turns
  ## around within the first ~100 time units, then logarithmic to TMAX.
  tt = sort(unique(c(seq(0, 300, length.out = 300),
                     exp(seq(log(300), log(TMAX), length.out = 120)))))

  one = function(i) {
    cc = P[i]
    trow = tg[.(cc$batch, cc$row_ind), nomatch = NULL]
    if (!nrow(trow)) return(NULL)
    p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
    p$Kg = cc$Kg; p$beta_ros = cc$beta_ros
    p$r_L_P = cc$r_L_P; p$r_LP_P = cc$r_L_P
    e0 = cc$Eh0_off + EPS * (EH0_HI - cc$Eh0_off)
    ## ode_ext: lsodar with the extinction root and the sweep's atol/rtol,
    ## trimmed back to tt (lib/extinction_floor.R, via skeleton_defs). Without
    ## the floor a cleared OFF arm can regrow out of numerical noise, which is
    ## exactly the difference between the groups this panel is drawn to show.
    arm = function(th) {
      o = tryCatch(ode_ext(seed(e0), tt, rhs_app, c(p, list(theta_treg = th)),
                           maxsteps = 1e5), error = function(e) NULL)
      if (is.null(o) || nrow(o) < length(tt) || !all(is.finite(o[nrow(o), ]))) return(NULL)
      as.data.table(o)[, theta := th][]
    }
    a = arm(0L); b = arm(1L)
    if (is.null(a) || is.null(b)) return(NULL)
    rbind(a[, .(id = i, grp = cc$treg_effect, theta, time, PLP)],
          b[, .(id = i, grp = cc$treg_effect, theta, time, PLP)])
  }
  RES = mclapply(seq_len(nrow(P)), function(i) tryCatch(one(i), error = function(e) NULL),
                 mc.cores = N_CORES)
  RES = RES[!vapply(RES, is.null, logical(1))]
  if (!length(RES)) stop("every integration failed")
  out = list(key = KEY, tr = rbindlist(RES, fill = TRUE))
  cat(sprintf("[ex] %d of %d combinations integrated under both conditions\n",
              uniqueN(out$tr$id), nrow(P)))
  saveRDS(out, CACHE); out
}

## Reuse only a cache built with the same settings, and say why not.
D = NULL
if (!RECOMP && file.exists(CACHE)) {
  cand = readRDS(CACHE)
  k = cand$key
  bad = if (is.null(k)) "it carries no settings key" else {
    d = names(KEY)[!vapply(names(KEY), function(z)
          isTRUE(all.equal(KEY[[z]], k[[z]])), logical(1))]
    if (length(d)) paste("these differ:", paste(d, collapse = ", ")) else ""
  }
  if (!nzchar(bad)) { cat("[ex] reusing", CACHE, "\n"); D = cand
  } else cat("[ex]", CACHE, "not reusable --", bad, "-- re-integrating\n")
}
if (is.null(D)) D = compute()

TR = as.data.table(D$tr)
GLEV = c("help", "harm"); GLAB = c("Help", "Harm")
TR[, grp := factor(grp, levels = GLEV, labels = GLAB)]
TR[, cond := factor(ifelse(theta == 0L, "Tregs OFF", "Tregs ON"),
                    levels = c("Tregs OFF", "Tregs ON"))]
TR[, PLP_disp := pmax(PLP, PLP_FLOOR)]
cat(sprintf("[ex] display floor %g | at the last output time %.1f%% of Tregs OFF arms and %.1f%% of Tregs ON arms are at or below it\n",
            PLP_FLOOR,
            100 * TR[cond == "Tregs OFF" & time == max(time), mean(PLP <= PLP_FLOOR)],
            100 * TR[cond == "Tregs ON"  & time == max(time), mean(PLP <= PLP_FLOOR)]))

## ---- the panel -------------------------------------------------------------
COL = setNames(c("#2a78d6", "#eb6834"), GLAB)
LTY = c(`Tregs OFF` = "22", `Tregs ON` = "solid")
## Each guide is one COLUMN, so Help sits above Harm and Tregs OFF above Tregs
## ON and the two columns stand side by side. legend.box = "vertical" would
## instead stack the guides and lay each one out along a row.
th_fig = theme_minimal(base_size = BASE_SZ) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.box = "horizontal", legend.box.just = "top",
        legend.margin = margin(b = 0), legend.spacing.x = unit(BASE_SZ * 0.05, "cm"),
        legend.text  = element_text(size = BASE_SZ * 0.95),
        legend.key.width = unit(BASE_SZ * 0.12, "cm"),
        axis.text    = element_text(size = BASE_SZ * 0.90),
        axis.title   = element_text(size = BASE_SZ * 1.00),
        panel.grid.minor = element_blank())
GUIDE1 = guides(colour   = guide_legend(ncol = 1, order = 1),
                fill     = guide_legend(ncol = 1, order = 1),
                linetype = guide_legend(ncol = 1, order = 2))
## Both conditions get a band; the OFF band is drawn fainter so the four regions
## stay separable where they overlap.
BP = TR[is.finite(PLP_disp), .(med = median(PLP_disp), lo = quantile(PLP_disp, .25),
                               hi = quantile(PLP_disp, .75)), by = .(grp, cond, time)]
keep = switch(BANDS, both = c("Tregs OFF", "Tregs ON"), on = "Tregs ON", none = character(0))
ribbons = lapply(keep, function(cd) geom_ribbon(data = BP[time > 0 & cond == cd],
  aes(ymin = lo, ymax = hi), colour = NA,
  alpha = if (cd == "Tregs OFF") 0.10 else 0.16, show.legend = FALSE))

pA = ggplot(BP[time > 0], aes(time, med, colour = grp, fill = grp, linetype = cond)) +
  ribbons +
  geom_line(linewidth = 0.9) +
  scale_x_log10() +
  ## THE AXIS IS NOT FORCED DOWN TO THE FLOOR. The floor clamps the data; making
  ## it the axis minimum as well would reserve thirteen decades of height for a
  ## region nothing is drawn in whenever no quartile reaches extinction, and
  ## compress the whole episode into the top of the panel. The axis is left to
  ## the quantiles actually plotted, so it reaches the floor exactly when some
  ## band does.
  scale_y_log10() +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
  scale_linetype_manual(values = LTY) +
  ## formatC("g") so a floor of 1 reads "1" rather than "1e+00" and a floor of
  ## 1e-11 still reads "1e-11"
  labs(x = "Time", y = bquote(P[LP] ~ "  (floored at " * .(formatC(PLP_FLOOR, format = "g")) * ")")) +
  th_fig + GUIDE1

if (nzchar(OUT_FIG)) {
  dir.create(dirname(OUT_FIG), showWarnings = FALSE, recursive = TRUE)
  png(OUT_FIG, width = FIG_W, height = FIG_H, res = FIG_RES)
  print(pA)
  invisible(dev.off())
  cat(sprintf("[ex] wrote %s\n", OUT_FIG))
}
cat(sprintf("   FOR THE CAPTION: the median P_LP over %d combinations per group at each\n",
            uniqueN(TR$id) %/% 2L))
cat("   output time, with the interquartile range as a band; blue is the group\n")
cat("   Tregs help, orange the group Tregs harm; dashed is Tregs OFF and solid\n")
cat("   Tregs ON; both arms of a combination start from the same initial\n")
cat(sprintf("   condition, and P_LP is floored at %g for display.\n", PLP_FLOOR))
