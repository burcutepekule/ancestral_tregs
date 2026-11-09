## ============================================================================
## Extended_Data_Fig_5.R -- Extended Data Fig. 5: the three classifications of
## Tregs-ON at beta_r* against Tregs-OFF at beta_r (first peak, integrated over the
## matching window, steady state) for every damage-matched pair of Fig. 4, as a
## Sankey of routes with the effect sizes above each axis. Same data and filter as the
## figure (lib/fig_damage_matched.R).
## THE PEAK AXIS IS THE FIRST PEAK of each arm -- the first time the pathogen turns over
## (pipeline/42_first_peak.R) -- and not peak_P, the maximum over the whole run. The two
## differ for a pair that settles on a high-burden attractor: its load climbs and stays,
## so its maximum IS the attractor, and comparing that with the other arm's genuine peak
## compared a steady state with a peak. On peak_P, 78 pairs were "more at the peak"; on
## the first peak there are none, and those 78 are pairs that end on a worse attractor,
## which the steady-state axis already carries. Fig. 4 measures the same quantity: its
## A and C clocks are each pair's own first peak.
## Drawn without ggalluvial: strata are stacked rectangles, each combination of
## classes is a ribbon between neighbouring axes with a cosine profile, stacked
## within its strata in a fixed order (eliminating pairs first) and drawn largest
## first so that thin routes stay visible.
## in : ./data/fig04_damage_matched_gated.rds, ./data/fig04_gated_units.rds
## out: ./figures/extended_data/Extended_Data_Fig_5.png
## RUN: Rscript figure_scripts/Extended_Data_Fig_5.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages({library(data.table); library(ggplot2)})
source("./lib/fig_damage_matched.R")
ELIM_TOL = as.numeric(Sys.getenv("ELIM_TOL", "1e-3"))
BASE_SZ  = as.numeric(Sys.getenv("BASE_SIZE", "34"))
OUT_FIG  = Sys.getenv("OUT_FIG", "./figures/extended_data/Extended_Data_Fig_5.png")
ORDER    = strsplit(Sys.getenv("ORDER", "peak,int,ss"), ",")[[1]]     # the axes, left to right
GAP      = 0.03      # gap between strata, as a share of all pairs
W        = 0.12      # half-width of a stratum bar
fmt = function(z) format(z, big.mark = ",", scientific = FALSE, trim = TRUE)

M = fig_dm_matched()                                   # the brake saves >= DAM_MIN of the damage
## the first peak of each arm (pipeline/42_first_peak.R), merged on the pair index
FP = Sys.getenv("F_PEAK", "./data/fig06_first_peak.rds")
if (!file.exists(FP)) stop(FP, " not found -- run: N_CORES=8 Rscript pipeline/42_first_peak.R")
P1 = as.data.table(readRDS(FP))[, .(i, peak1_off, peak1_onS, no_transient_off, no_transient_onS)]
M = merge(M, P1, by = "i")
if (M[is.na(peak1_off) | is.na(peak1_onS), .N])
  cat(sprintf("[sankey] %d pairs without a first peak in one arm are dropped\n", M[is.na(peak1_off) | is.na(peak1_onS), .N]))
M = M[!is.na(peak1_off) & !is.na(peak1_onS)]
cat(sprintf("[sankey] no transient peak (the first peak is the final load): OFF %d, ON* %d\n",
            M[no_transient_off == TRUE, .N], M[no_transient_onS == TRUE, .N]))
M[, `:=`(el_off = is.finite(t_clear_off) | P_end_off < ELIM_TOL, el_on = is.finite(t_clear_onS) | P_end_onS < ELIM_TOL)]
M[, fate := fifelse(el_off & el_on, "both_elim", fifelse(!el_off & !el_on, "both_persist", fifelse(el_on, "only_on", "only_off")))]
M = M[fate %in% c("both_persist", "only_on")]
M[, `:=`(peak = peak1_onS < peak1_off, int = P_win_onS < P_win_off, ss = P_end_onS < P_end_off)]   # TRUE = less with Tregs
F = M[, .(N = .N), by = c("fate", ORDER)]
cat(sprintf("[sankey] %d pairs: %d neither eliminates, %d only Tregs-ON eliminates\n", nrow(M), M[fate == "both_persist", .N], M[fate == "only_on", .N]))
print(F[order(fate, -N)])
for (a in ORDER) cat(sprintf("  %-5s less with Tregs in %.1f%%\n", a, 100 * M[, mean(get(a))]))

## ---- layout ---------------------------------------------------------------------
## Within every stratum ("less" on top, "more" below), the ribbons are stacked in one
## fixed order: by fate (eliminating pairs first), then by the classes at the other
## axes. The same order on both sides of a gap keeps ribbons from crossing needlessly.
nA = length(ORDER); tot = sum(F$N)
setorderv(F, c("fate", ORDER), order = c(-1L, rep(-1L, nA)))      # only_on first; TRUE (less) first
F[, id := .I]
## stratum bottoms per axis: "less" (TRUE) above "more" (FALSE), separated by GAP
pos = list()
for (k in seq_len(nA)) {
  a = ORDER[k]
  nL = F[get(a) == TRUE, sum(N)]; nM = tot - nL
  y0 = c(`TRUE` = nM + GAP * tot, `FALSE` = 0)
  P = copy(F)[, cls := get(a)]
  setorderv(P, c("cls", "fate", ORDER), order = c(-1L, -1L, rep(-1L, nA)))
  P[, `:=`(y_lo = y0[as.character(cls)] + cumsum(N) - N), by = cls][, y_hi := y_lo + N]
  pos[[k]] = P[, .(id, axis = k, cls, y_lo, y_hi)]
}
POS = rbindlist(pos)
STR = POS[, .(y_lo = min(y_lo), y_hi = max(y_hi), N = sum(y_hi - y_lo)), by = .(axis, cls)]
STR[, `:=`(lab = sprintf("%s\n%.1f%%", fifelse(cls, "Less", "More"), 100 * N / tot), small = N / tot < 0.04)]
## ribbons: between axis k and k+1, a cosine profile from the left interval to the right
tseq = seq(0, 1, length.out = 60); ease = (1 - cos(pi * tseq)) / 2
RIB = rbindlist(lapply(seq_len(nA - 1), function(k) {
  Lk = POS[axis == k]; Rk = POS[axis == k + 1]
  J = merge(Lk[, .(id, y_lo_l = y_lo, y_hi_l = y_hi)], Rk[, .(id, y_lo_r = y_lo, y_hi_r = y_hi)], by = "id")
  J[, .(seg = k, x = k + W + (1 - 2 * W) * tseq,
        ymin = y_lo_l + (y_lo_r - y_lo_l) * ease, ymax = y_hi_l + (y_hi_r - y_hi_l) * ease), by = id]
}))
## one colour per path: the eliminations in magenta, the persistent pairs by their
## route through the three classes (where "less" = L and "more" = M), in a cold-to-warm
## order from all-less to all-more; routes below MIN_SHARE of pairs are pooled as "other"
## every route is its own entry, however small (no pooling)
F[, route := paste(ifelse(get(ORDER[1]), "L", "M"), ifelse(get(ORDER[2]), "L", "M"), ifelse(get(ORDER[3]), "L", "M"), sep = "")]
F[, key := fifelse(fate == "only_on", "elim", route)]
ROUTES = c("LLL", "LLM", "LML", "LMM", "MLL", "MLM", "MML", "MMM")
RLAB = function(r) paste(ifelse(strsplit(r, "")[[1]] == "L", "less", "more"), collapse = " / ")
## legend order: the eliminations, then the routes from all-less to all-more (as ROUTES lists them)
KEYS = c("elim", intersect(c("LLL", "LLM", "LML", "LMM", "MLL", "MLM", "MML", "MMM"), F[key != "elim", unique(key)]))
KL = sapply(KEYS, function(k) if (k == "elim") sprintf("Only Tregs-ON eliminates (n = %s)", fmt(F[key == "elim", sum(N)]))
                              else sprintf("Neither eliminates, %s (n = %s)", RLAB(k), fmt(F[key == k, sum(N)])))
## the eliminations in navy; the persistent routes, in the order of ROUTES, from dark
## green (less everywhere) through the greens to light and dark orange (more everywhere)
## one colour per route, keyed by the route itself so a colour never changes with the sample
PAL_R = c(LLL = "#9ECAE1", LLM = "#5AAE61", LML = "#1B7837", LMM = "#FDB863", MMM = "#E8220E",
          MLL = "#9E9AC8", MLM = "#6A51A3", MML = "#E7298A")
RK = KEYS[!KEYS %in% c("elim", "other")]
COL_F = setNames(c("#0A3BFF", PAL_R[RK]), KL[c("elim", RK)])
RIB = merge(RIB, F[, .(id, key, N)], by = "id")
RIB[, fate := factor(KL[key], levels = KL[KEYS])]
## drawing order: the largest ribbons first (at the back), the thin ones on top
setorder(RIB, -N, id, seg, x)
RIB[, grp := factor(paste(id, seg), levels = unique(paste(id, seg)))]
AXL = c(peak = "First local\nmaximum", int = "Integrated\n(matching window)", ss = "Steady state")

g = ggplot() +
  geom_ribbon(data = RIB, aes(x = x, ymin = ymin, ymax = ymax, group = grp, fill = fate), alpha = 0.8,
              colour = "black", linewidth = 0.35) +
  geom_rect(data = STR, aes(xmin = axis - W, xmax = axis + W, ymin = y_lo, ymax = y_hi), fill = "white", colour = "grey25", linewidth = 0.6) +
  geom_text(data = STR[small == FALSE], aes(x = axis, y = (y_lo + y_hi) / 2, label = lab), size = BASE_SZ * 0.3, lineheight = 0.9) +
  geom_text(data = STR[small == TRUE], aes(x = axis + W + 0.03, y = (y_lo + y_hi) / 2, label = sub("\n", " ", lab)),
            size = BASE_SZ * 0.26, hjust = 0) +
  scale_fill_manual(values = COL_F, name = NULL) +
  scale_x_continuous(breaks = seq_len(nA), labels = AXL[ORDER], limits = c(0.5, nA + 0.5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = BASE_SZ) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(size = BASE_SZ * 0.9, face = "bold", colour = "grey20"),
        legend.position = "right", legend.text = element_text(size = BASE_SZ * 0.7), legend.key.size = unit(BASE_SZ * 0.05, "cm"),
        legend.title = element_text(size = BASE_SZ * 0.75)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 1),
                             title = sprintf("Route:\n%s", paste(sub("First local", "First peak", sub("\n.*", "", AXL[ORDER])), collapse = " / "))))

## ---- the effect sizes: one histogram per axis above the Sankey ---------------------------
## Tregs-OFF at beta_r minus Tregs-ON at beta_r*, positive = less pathogen with Tregs: the
## peak and the steady state in P_LP units (what matters there is the level itself), the
## integral as % of Tregs-OFF. Histograms, not densities: a bin edge sits at 0, so nothing
## is smeared across the sign. Symmetric log axes; one neutral fill, so the colours of the
## routes below are not reused; the median of each side is labelled.
tr = scales::pseudo_log_trans(sigma = 1, base = 10)
## every pair carries the key of its route, so the histograms can be stacked by route in
## the same colours as the ribbons below
M[, route := paste(ifelse(get(ORDER[1]), "L", "M"), ifelse(get(ORDER[2]), "L", "M"), ifelse(get(ORDER[3]), "L", "M"), sep = "")]
M = merge(M, unique(F[, .(fate, route, key)]), by = c("fate", "route"))
M[, klab := factor(KL[key], levels = KL[KEYS])]
EFFS = list(peak = M[, .(v = peak1_off - peak1_onS, klab)], int = M[, .(v = 100 * (1 - P_win_onS / P_win_off), klab)], ss = M[, .(v = P_end_off - P_end_onS, klab)])
EXL = c(peak = "atop('First local maximum of '*P[LP], 'Tregs-OFF '*minus*' Tregs-ON (units)')", int = "atop('Integrated '*P[LP], 'reduction with Tregs (%)')", ss = "atop('Steady-state '*P[LP], 'Tregs-OFF '*minus*' Tregs-ON (units)')")
hist1 = function(a) {
  D = EFFS[[a]]; r = range(D$v)
  ## symmetric about 0 when both signs occur; from 0 upward when every difference is positive
  ## (the first peak), so no empty negative half is drawn
  bt = if (min(D$v) >= 0) seq(0, max(tr$transform(r)), length.out = 31)
       else seq(-max(abs(tr$transform(r))), max(abs(tr$transform(r))), length.out = 61)   # 0 is a bin edge
  brk = tr$inverse(bt)
  xb = if (a == "int") c(-100, -10, 0, 10, 100) else c(-30, -10, -3, -1, 0, 1, 3, 10, 30)
  xb = xb[xb >= min(brk) & xb <= max(brk)]
  lb = if (a == "int") function(z) paste0(z, "%") else function(z) z
  MD = D[, .(med = median(v), n = .N), by = .(side = v > 0)]
  ## one label for the whole distribution: its median and its range (signed; + = less with Tregs)
  u = if (a == "int") "%" else ""
  f2 = function(z) paste0(ifelse(z < 0, "\u2212", ""), format(signif(abs(z), 2)), u)
  MD = D[, .(side = TRUE, lab = sprintf("median %s, range %s to %s", f2(median(v)), f2(min(v)), f2(max(v))))]
  MD[, `:=`(y = Inf, vj = 1.2)]
  ggplot(D, aes(x = v)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    geom_histogram(aes(fill = klab), breaks = brk, colour = "white", linewidth = 0.1, position = "stack", show.legend = FALSE) +
    scale_fill_manual(values = COL_F) +
    geom_text(data = MD, aes(x = -Inf, y = y, label = lab, vjust = vj), hjust = -0.05, size = BASE_SZ * 0.26, colour = "grey15") +
    ## oob_keep: the outermost bin edges equal the limits up to rounding after the
    ## transform/inverse round trip, and censoring would drop those bins
    scale_x_continuous(trans = tr, breaks = xb, labels = lb, limits = range(brk), oob = scales::oob_keep) +
    scale_y_continuous(breaks = NULL, expand = expansion(mult = c(0, 0.22))) +      # headroom for the median label
    labs(x = parse(text = EXL[[a]]), y = NULL) +
    theme_minimal(base_size = BASE_SZ) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = BASE_SZ * 0.7), axis.title.x = element_text(size = BASE_SZ * 0.72))
}
pd = Reduce(`|`, lapply(ORDER, hist1))
suppressMessages(library(patchwork))
## the legend is collected at the composite level and sits to the right of BOTH rows,
## so the Sankey's panel spans the same width as the three histograms above it
## the legend is pulled left into the Sankey's empty right margin (the strip right of the
## steady-state bar that aligns it with the histograms above)
g = (pd / g) + plot_layout(heights = c(0.42, 1), guides = "collect") &
  theme(legend.position = "right", legend.box.spacing = unit(-BASE_SZ * 7.5, "pt"), legend.justification = c(0, 0.5))
dir.create(dirname(OUT_FIG), showWarnings = FALSE, recursive = TRUE)
png(OUT_FIG, width = 6600, height = 3800, res = 200); print(g); invisible(dev.off())
cat(sprintf("[sankey] wrote %s\n", OUT_FIG))
