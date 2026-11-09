## ============================================================================
## Supplementary_Fig_3.R -- SUPPLEMENTARY FIG. 3: the
## smallest intervention that changes the outcome, against when the pathogen was
## brought under control.
##
## WHAT IS PLOTTED, AND WHY THE SMALLEST. For every harmed combination and every
## intervention, the levels are taken in order from the untreated arm outward and
## the FIRST level at which the combination reaches each outcome class is
## recorded. The figure then reports the SMALLEST such level within each group of
## combinations. That is an existence statement: the least intervention for which
## the outcome is possible at all. It is not a summary of how many combinations
## reach it, and the panels carry no counts; the counts are in the output table.
##
## TWO DISJOINT OUTCOME CLASSES.
##   harm -> neutral   |final burden difference| <= EFFECT_TOL, the Tregs ON arm
##                     ends where the Tregs OFF arm ends
##   harm -> help      final burden difference <  -EFFECT_TOL, the Tregs ON arm
##                     ends BELOW the Tregs OFF arm
## A combination is in one class or the other at a given level, never both, so
## the two are plotted as separate series and one may be absent where the other
## is not.
##
## NOTHING IS BINNED AND NOTHING IS MINIMISED OVER A GROUP. Every panel is one
## point per combination: its own least sufficient level against its own t_dip.
## Taking a minimum over a group of combinations breaks the pairing that the
## figure is about -- a single combination with a short dip time and a small
## sufficient level sets the value for every combination grouped with it -- so
## the pairing is kept and the reader sees the relation directly.
##
## THE X AXIS IS t_dip: P_LP rises to a peak and falls, and t_dip is the time of
## the minimum that follows. DIP_ARM selects the arm it is read from: "off"
## (default) makes it a property of the combination, unaffected by any
## intervention; "on" reads it from the untreated Tregs ON arm. Where the arm
## eliminates the pathogen the fall is monotone and t_dip is the end of the
## integration, which is reported so those combinations can be recognised.
##
## PANEL A CARRIES THE IDENTITY LINE because there the two axes are both times:
## the delay and the dip time of the same combination. A point above the line is
## a combination for which the delay that first produces help is longer than its
## own dip time.
##
## PANEL A USES THE BISECTED THRESHOLD. The delays that were integrated are a
## coarse sequence, so the first one that produces help overstates the threshold
## by up to the gap to the previous level. pipeline/50_harm_rescue_compute.R bisects that interval
## for every combination that has one, and the result travels in the same file
## as the curves; this script plots it on a continuous axis, where the identity
## line is exact. If the file carries no thresholds, panel A falls back to the
## grid level and says so at run time.
##
## in : ./data/harm_rescue_data.rds, written by pipeline/50_harm_rescue_compute.R.
##      ./data/harm_explore_cache.rds, written by pipeline/51_harm_explore_sample.R.
##      ./data/skeleton_pathogenic_harm.rds and the tau grids (row 2 only). This script integrates
##      nothing. Produce the file over the whole gated class first:
##        SAMPLE=0 N_CORES=12 Rscript pipeline/50_harm_rescue_compute.R
## out: ./figures/supplementary_information/Supplementary_Fig_3.png
##      ./source_data/Supplementary_Fig_3_dose_table.csv
##      ./source_data/Supplementary_Fig_3_delay_fit.csv
## RUN: Rscript figure_scripts/Supplementary_Fig_3.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

IN_RDS  = Sys.getenv("IN_RDS", "./data/harm_rescue_data.rds")
TOL     = as.numeric(Sys.getenv("EFFECT_TOL", "1.0"))
ZTOL    = as.numeric(Sys.getenv("ZTOL", "1e-3"))     # below this the load is zero
DIP_ARM = Sys.getenv("DIP_ARM", "off")               # "off" or "on"
NBIN    = as.integer(Sys.getenv("NBIN", "4"))        # bins of the P_LP endpoint
NDIP    = as.integer(Sys.getenv("NDIP", "5"))        # bins of t_dip
OUT_FIG = Sys.getenv("OUT_FIG", "./figures/supplementary_information/Supplementary_Fig_3.png")
OUT_TAB = Sys.getenv("OUT_TAB", "./source_data/Supplementary_Fig_3_dose_table.csv")
OUT_FIT = Sys.getenv("OUT_FIT", "./source_data/Supplementary_Fig_3_delay_fit.csv")
FIG_W   = as.numeric(Sys.getenv("FIG_W", "5400"))
FIG_H   = as.numeric(Sys.getenv("FIG_H", "1900"))
FIG_RES = as.numeric(Sys.getenv("FIG_RES", "200"))
## the two-row canvas is its own size: five panels across the bottom row and
## three across the top, so it is wider and twice as tall as the old one
FIG_W2  = as.numeric(Sys.getenv("FIG_W2", "7200"))
FIG_H2  = as.numeric(Sys.getenv("FIG_H2", "4000"))
## EVERY PIECE OF TEXT IN THE FIGURE IS EXPRESSED AGAINST BASE_SZ -- the theme
## base size, the panel A fit block, the direct labels in C, the tag letters and
## the legend key width -- so this one number sets the whole figure's type scale
## and nothing has to be adjusted alongside it. Raised from 19: at the canvas
## size below, 19pt rendered at roughly 7pt on a printed column and the legends
## were the worst of it.
BASE_SZ = as.numeric(Sys.getenv("BASE_SIZE", "34"))
## ONE MULTIPLIER ON EVERY DRAWN CURVE in the ECDF panel and the second row.
## The base widths differ between them on purpose -- the untreated arm in the
## second row is drawn heavier than the four treated ones so it reads as the
## reference -- so this scales them together instead of replacing them, and the
## relative weighting survives any change. It does not touch the identity line
## in panel B or the explore bands in A, which are references rather than data.
LW_MULT = as.numeric(Sys.getenv("LW_MULT", "1.5"))
SIZE_U  = as.numeric(Sys.getenv("SIZE_U", "3.4"))    # mark size on the three rate panels
## Panel A carries a dozen points where the rate panels carry hundreds, so the
## same mark size reads as sparse there. Larger, and set separately rather than
## scaled from SIZE_U, so changing one does not silently move the other.
SIZE_A  = as.numeric(Sys.getenv("SIZE_A", "5.4"))    # mark size on panel A
## CRIT_EXT=1 colours panel A by extinction instead of by outcome class, which
## is what it did while every point on it was a harm-to-help switch. Off by
## default: d_star is now the end of harm, so the class is informative there.
CRIT_EXT = as.integer(Sys.getenv("CRIT_EXT", "0")) == 1L
## CONT_RATE=1 draws B, C and D on the bisected multiplier instead of the coarse
## grid level. The thresholds for all four treatments are in the data file
## either way; this only decides whether the three rate panels use theirs.
CONT_RATE = as.integer(Sys.getenv("CONT_RATE", "0")) == 1L
## ECDF=0 restores the load-at-dip scatter in B, C and D. Defined here rather
## than beside the panel that uses it because panel() reads it too, to decide
## whether to contribute a legend.
USE_ECDF = as.integer(Sys.getenv("ECDF", "1")) == 1L

if (!file.exists(IN_RDS))
  stop(sprintf("%s not found; run pipeline/50_harm_rescue_compute.R first, or point IN_RDS at an existing data file", IN_RDS))
IN = readRDS(IN_RDS)
SM = as.data.table(IN$sm)
## ONE FILE, SO NOTHING TO RECONCILE. This used to read three -- the curves, the
## dip statistics and the delay thresholds -- each carrying a settings stamp
## that was compared against the others. The comparison could not do the job
## asked of it: the stamp holds only numerical settings, none of which change
## between one sweep and the next, so a companion file left over from a previous
## dataset passed it cleanly. The three are now one list written in a single
## pass, and the stamp is reported rather than checked.
STAMP0 = IN$settings
if (!is.null(STAMP0)) {
  cat("[Supp. Fig. 3] settings: ")
  cat(paste(sprintf("%s=%s", c("ATOL", "RTOL", "TMAX", "HMAX", "EXT_TOL", "EFFECT_TOL"),
                    c(STAMP0$ATOL, STAMP0$RTOL, STAMP0$TMAX, STAMP0$HMAX,
                      STAMP0$EXT_TOL, STAMP0$EFFECT_TOL)), collapse = "  "), "\n")
}
if (!"t_dip" %in% names(SM))
  stop(sprintf("%s predates the t_dip columns; re-run pipeline/50_harm_rescue_compute.R", IN_RDS))
KK = c("batch", "row_ind", "Kg", "beta_ros", "r_L_P")
SM[, effect := fifelse(dcomp < -TOL, "help", fifelse(dcomp > TOL, "harm", "neutral"))]

## THE LEVELS OF EACH TREATMENT, ordered from the untreated arm outward. sgn
## gives the direction of increasing intervention, so that step index 1 is the
## untreated arm in all four.
TRT = list(
## THE VALUE REPORTED IS THE LEAST SUFFICIENT INTERVENTION in every panel, which
## is the FIRST level in the ordering from the untreated arm outward. For tau_10
## and tau_t,a the intervention is a reduction, so the least sufficient one is
## the LARGEST multiplier at which the class is reached -- the most suppression
## the combination still tolerates. For tau_t,d it is the smallest multiplier and
## for the delay the shortest delay. The axis labels say which, because "smallest
## multiplier" and "least intervention" are opposite statements for two of the
## four.
  delay    = list(sgn =  1, name = "Shortest delay that suffices"),
  strength = list(sgn = -1, name = expression("Largest " * tau[10] * " multiplier that suffices")),
  deact    = list(sgn =  1, name = expression("Smallest " * tau[td] * " multiplier that suffices")),
  induce   = list(sgn = -1, name = expression("Largest " * tau[ta] * " multiplier that suffices")))
for (tt in names(TRT))
  TRT[[tt]]$lv = sort(unique(SM[treat == tt, level]), decreasing = TRT[[tt]]$sgn < 0)
if (any(vapply(TRT, function(z) length(z$lv) == 0L, logical(1))))
  stop(IN_RDS, " does not carry all four treatments")

## ONLY COMBINATIONS THAT ARE ACTUALLY HARMED UNTREATED. The class label comes
## from the parameter sweep; a few combinations fall on the other side of the
## tolerance once integrated, and a smallest level is not defined for those.
UNT = SM[treat == "delay" & level == TRT$delay$lv[1]]
HARMED = UNT[effect == "harm", ..KK]
cat(sprintf("[Supp. Fig. 3] %d combinations in the file, %d harmed when integrated untreated, %d dropped\n",
            nrow(UNT), nrow(HARMED), nrow(UNT) - nrow(HARMED)))

## THE PER-COMBINATION COVARIATES, read from a single arm so that they describe
## the combination and not the intervention.
##
## THE DIP STATISTICS COME FROM THE SUMMARY ROWS, and t_zero with them. There
## used to be a second file for these, dip_stats.rds, on the grounds that the
## statistics stored here had been computed on the raw trajectory, where a
## monotone decline through values with no meaning puts the minimum wherever
## the integration stops. That stopped being true when dipstat() started
## flooring P_LP before taking the minimum: run() calls it on every arm, so
## these columns ARE the floored statistics and the second file recomputed them
## from a fresh integration of the same two arms under the same settings.
## t_zero was the one column the fallback path did not carry, which quietly set
## it to NA and left panel A's extinction colouring resting on the dip test
## alone. It is selected here.
ARM = if (DIP_ARM == "off") SM[treat == "off"] else UNT
COV = merge(HARMED, ARM[, c(KK, "t_dip", "PLP_dip", "t_ctrl", "t_zero",
                            "PLP_end", "PL_end"), with = FALSE], by = KK)
setnames(COV, "PLP_end", "PLP_arm")
cat(sprintf("[Supp. Fig. 3] %d of %d harmed combinations reach the floor on that arm and have a\n     defined elimination time\n",
            sum(!is.na(COV$t_zero)), nrow(COV)))
cat(sprintf("[Supp. Fig. 3] t_dip read from the Tregs %s arm\n", toupper(DIP_ARM)))
cat(sprintf("[Supp. Fig. 3] of the harmed combinations, %d end that arm with P_LP below %g and %d with PL below %g\n",
            sum(COV$PLP_arm < ZTOL), ZTOL, sum(COV$PL_end < ZTOL), ZTOL))

SM = merge(SM, HARMED, by = KK)

## THE SMALLEST LEVEL AT WHICH EACH CLASS IS REACHED, per combination.
dose = function(tt) {
  lv = TRT[[tt]]$lv
  D = SM[treat == tt]
  D[, step := match(level, lv)]
  setorderv(D, c(KK, "step"))
  D[, .(step_neutral = { i = which(effect == "neutral"); if (length(i)) step[i[1]] else NA_integer_ },
        step_help    = { i = which(effect == "help");    if (length(i)) step[i[1]] else NA_integer_ },
        nm_neutral = { i = which(effect == "neutral"); length(i) > 0 && length(i) != length(step) - i[1] + 1 },
        nm_help    = { i = which(effect == "help");    length(i) > 0 && length(i) != length(step) - i[1] + 1 }),
    by = KK]
}
DOSE = merge(rbindlist(lapply(names(TRT), function(tt) cbind(treat = tt, dose(tt)))),
             COV, by = KK)
cat("\n[Supp. Fig. 3] reached at any level, and reached then left again at a larger level:\n")
print(DOSE[, .(n = .N, reach_neutral = sum(!is.na(step_neutral)),
               reach_help = sum(!is.na(step_help)),
               left_neutral_again = sum(nm_neutral), left_help_again = sum(nm_help)),
           by = treat])

## ---- THE REFINED DELAY THRESHOLD, if it has been computed -----------------
DSTAR = NULL
DELAY_FIT = NULL      # filled by panel A when the refined threshold is available
if (!is.null(IN$thresh) && nrow(IN$thresh)) {
  ## one row per combination PER TREATMENT PER BOUNDARY. `bound` is "neutral"
  ## where harm ends and "help" where the arm beats no Tregs at all, which is
  ## exactly the distinction `crit` carries, so the join below is on the class
  ## the row is about and each point gets the threshold for its own class.
  DSTAR = as.data.table(IN$thresh)[, c(KK, "treat", "bound", "thr", "nonmono"), with = FALSE]
  cat(sprintf("\n[Supp. Fig. 3] bisected thresholds: %d rows over %d treatments, %d non-monotone\n",
              nrow(DSTAR), uniqueN(DSTAR$treat), sum(DSTAR$nonmono)))
  print(DSTAR[, .N, by = .(treat, bound)][order(treat, bound)])
} else {
  cat("\n[Supp. Fig. 3] the data file carries no bisected thresholds, so panel A shows the\n")
  cat("     coarse grid level and the identity line means nothing on it. That\n")
  cat("     happens when no combination reached help at any ladder delay, which\n")
  cat("     on a small SAMPLE is ordinary and on the full set would be a result.\n")
}

## ---- THE LONG FORM: one row per combination, treatment and class -----------
## the depth of the dip, used as the mark size. lsoda can return a small
## negative value at the floor, so it is clamped at zero, and the scale is
## shared across panels so that a mark of a given size means the same thing in
## all four.
DOSE[, dip := pmax(PLP_dip, 0)]
DIPMAX = max(DOSE$dip, na.rm = TRUE)
cat(sprintf("[Supp. Fig. 3] P_LP at the dip on the Tregs %s arm: median %.3g, max %.3g, %d of %d below %g\n",
            toupper(DIP_ARM), median(DOSE[treat == "delay", dip]), DIPMAX,
            sum(DOSE[treat == "delay", dip] < ZTOL), nrow(DOSE[treat == "delay"]), ZTOL))
L = melt(DOSE, id.vars = c(KK, "treat", "t_dip", "PLP_dip", "dip", "t_ctrl", "t_zero", "PLP_arm"),
         measure.vars = c("step_neutral", "step_help"),
         variable.name = "crit", value.name = "step")
L[, crit := fifelse(crit == "step_neutral", "Harm to neutral", "Harm to help")]
## EXTINCTION IS ITS OWN CATEGORY, AND IT IS SHOWN IN PANEL A ONLY. A combination
## whose arm drives the pathogen below the floor is a different object from one
## that merely dips: the floor is absorbing, so the load cannot return, and that
## is the condition under which a delay can help at all. t_zero is defined exactly
## when the floor was reached, and the dip value is the fallback for files written
## before it was recorded.
##
## IT IS NO LONGER USED FOR PANEL A, and the reason it was is worth keeping.
## While d_star was bisected on the onset of help, every point in panel A was a
## harm-to-help switch by construction, so the outcome class carried no
## information there and the colour was free to say something else -- whether
## the pathogen was eliminated. d_star is now bisected on the END OF HARM, so
## panel A also holds the combinations a delay lifts out of harm without
## carrying them past the help threshold, and the class is informative again.
## Panel A therefore colours by `crit` like the other three, and `crit_ext` is
## kept for the extinct counts in the table and for CRIT_EXT=1, which restores
## the old colouring. Every point in panel A is still extinct -- that is the
## all-or-nothing result -- so nothing is lost by not colouring it: it is
## checked below and stated in the caption.
L[, extinct := (!is.na(t_zero)) | dip < ZTOL]
L[, crit := factor(crit, levels = c("Harm to neutral", "Harm to help"))]
L[, crit_ext := factor(fifelse(extinct, "Pathogen extinct", as.character(crit)),
                       levels = c("Harm to neutral", "Harm to help", "Pathogen extinct"))]
## the boundary each row's class corresponds to, so a row can be joined to the
## threshold that was bisected for THAT class rather than to a single one
L[, bound := fifelse(crit == "Harm to neutral", "neutral", "help")]
L = L[!is.na(step) & is.finite(t_dip)]
cat(sprintf("\n[Supp. Fig. 3] %d of %d rows are combinations whose Tregs %s arm reaches the floor\n",
            sum(L$extinct), nrow(L), toupper(DIP_ARM)))
## The caption says panel A's points are all extinct harm-to-help switches. That
## is a claim about the data, not a design decision, so it is checked here rather
## than assumed -- on a different dataset it may simply stop being true.
##
## IT IS CHECKED ON THE POINTS PANEL A DRAWS, WHICH ARE NOT ALL THE DELAY ROWS.
## Panel A plots d_star, and d_star exists only where the bisection had a help
## level to bracket against a harmful one. A combination whose best outcome at
## any delay is NEUTRAL has no such level, so it has no d_star and is not drawn
## -- the full sweep contains one of these. Checking every delay row instead
## would report two classes in a panel that shows one, which is how this check
## first went wrong.
.dl = if (is.null(DSTAR)) L[treat == "delay"] else
        merge(L[treat == "delay"], DSTAR[treat == "delay", -"treat"],
              by = c(KK, "bound"))[!is.na(thr)]
if (nrow(.dl) && (any(!.dl$extinct) || uniqueN(.dl$crit) > 1L))
  cat(sprintf("[Supp. Fig. 3] NOTE panel A is no longer uniform: %d of %d extinct, classes %s -- the caption needs rewording\n",
              sum(.dl$extinct), nrow(.dl), paste(unique(as.character(.dl$crit)), collapse = ", ")))
cat(sprintf("[Supp. Fig. 3] panel A draws %d points; %d delay rows have no bisected threshold and are not drawn\n",
            nrow(.dl), nrow(L[treat == "delay"]) - nrow(.dl)))

## Reported rather than drawn: whether the least sufficient level tracks the
## combination's own dip time. This is the quantity panel A is about.
cat("\n[Supp. Fig. 3] least sufficient level against t_dip, per combination:\n")
print(L[, .(n = .N, n_extinct = sum(extinct),
            spearman = round(suppressWarnings(cor(step, t_dip, method = "spearman")), 3)),
        by = .(treat, crit)][order(treat, crit)])

## ---- COLOURS AND THEME ----------------------------------------------------
## Muted blue-grey and pale green, checked with the dataviz validator: lightness
## in band and adjacent CVD separation 15.5 deutan, 12.3 tritan. Both sit below
## 3:1 against the page, which the dark outline on every mark supplies instead.
## Muted blue-grey, pale green and a red for extinction, checked with the dataviz
## validator: all three inside the lightness band, worst adjacent CVD separation
## 15.5 deutan and 12.3 tritan. All three sit below 3:1 against the page, which
## the dark outline on every mark supplies instead.
## The blue-grey and the pale green sat too close for a reader scanning four
## panels, so the neutral class is now a magenta. Chosen with the dataviz
## validator rather than by eye, on the three colours together:
##   lightness band  PASS, all three inside L 0.43-0.77
##   chroma floor    PASS, all three >= 0.1 -- the old #7FC08D was 0.099 and
##                   read as grey, which is why the green is deepened too
##   CVD separation  PASS, magenta vs green deltaE 15.2 deutan / 23.3 tritan
##   normal vision   PASS, magenta vs green 29.7
##   contrast        WARN, the green is 2.4:1 against the page. Every mark
##                   carries a dark outline, which is the relief that warning
##                   asks for; the legend is also labelled, so identity is never
##                   colour alone.
##
## The extinct colour is near-black rather than the red it started as. It is the
## ONLY class drawn in panel A, so it is not a categorical slot competing with
## the other two and the CVD checks above do not constrain it; what it has to be
## is legible on a sparse panel, and #1C1C1C against the page is the strongest
## contrast available. It never appears beside the magenta or the green.
CRIT_COL = c("Harm to neutral"  = "#2BBBD7",
             "Harm to help"     = "#7DCCAD",#A5D6A7
             "Pathogen extinct" = "#1C1C1C")
## THE LEGEND TEXT IS SET EXPLICITLY. ggplot2's themes give legend.text
## rel(0.8) of the base size and the legend title rel(1.0), so raising
## BASE_SZ alone leaves the keys reading smaller than the axis titles beside
## them -- and in this figure the legend is where the treatment names live, so
## it is the text that has to be readable. LEG_FRAC is a fraction of BASE_SZ
## rather than an absolute size, so it follows any change to the type scale.
LEG_FRAC = as.numeric(Sys.getenv("LEG_FRAC", "1"))
th = theme_minimal(base_size = BASE_SZ) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        legend.text  = element_text(size = BASE_SZ * LEG_FRAC),
        legend.title = element_text(size = BASE_SZ * LEG_FRAC),
        legend.key.height = unit(BASE_SZ * 0.04, "cm"),   # one row pitch in every legend
        ## panel letters: here and not in plot_annotation(), whose theme does not
        ## reach the tags of nested patchworks. 1.3 gives them the share of the page
        ## width that Supplementary Fig. 2's letters have.
        plot.tag = element_text(size = BASE_SZ * 1.3, face = "bold"))

## the vertical jitter is a fraction of the spacing between two levels, so a
## point cannot be read as lying at a level it does not occupy
JIT = 0.13
panel = function(tt, identity_line = FALSE) {
  lv = TRT[[tt]]$lv
  D = L[treat == tt]
  ## panel A on a continuous axis when the bisected threshold is available: the
  ## quantity is then a delay in time units and not an index into the sequence
  ## CONT_RATE puts the bisected multiplier on B, C and D too. Off by default:
  ## those panels have always shown the coarse level and switching them is a
  ## change to three figures, not a bug fix. The thresholds are computed either
  ## way, so it costs nothing to leave off and nothing to turn on.
  refined = !is.null(DSTAR) && (identity_line || CONT_RATE) &&
            nrow(DSTAR[treat == tt])
  if (refined) {
    D = merge(D, DSTAR[treat == tt, -"treat"], by = c(KK, "bound"))[!is.na(thr)]
    D[, yy := thr]
  } else {
    D[, yy := step]
  }
  ## ---- WHAT GOES ON THE X AXIS, AND WHY IT DIFFERS BY PANEL ----------------
  ## Panel A puts the TIME of the dip on x, because its claim is an identity
  ## between two times: the delay at which help begins equals the time the Tregs
  ## OFF arm clears the pathogen. Both axes have to be times for that line to
  ## mean anything.
  ##
  ## The three rate panels put the DEPTH of the dip on x instead. Their question
  ## is not when the pathogen was at its lowest but how low it got -- how much
  ## suppression a combination can absorb depends on the load that survives, not
  ## on the clock. Time of the dip carried no signal in those panels: the points
  ## spread across four decades of x at every level of y.
  ##
  ## That also retires the size aesthetic. Size encoded P_LP at the dip, which is
  ## now the x position, so keeping it would encode one quantity twice and would
  ## draw the extinct combinations -- the ones at P_LP = 0 -- as the least visible
  ## marks on the panel. One uniform size, no size guide.
  D[, xx := if (identity_line) t_dip else dip]
  ## Panel A draws the three-level scale so extinction is visible; the rate panels
  ## draw the two-level outcome scale. Only the rate panels contribute a fill
  ## guide. That is also what stops the fill legend collecting twice: patchwork
  ## compares guides including the layer that produced them, and panel A's layer
  ## has no size aesthetic, so its fill guide was never identical to the others'
  ## and the two were collected side by side.
  D[, .fill := if (identity_line && CRIT_EXT) crit_ext else crit]
  jit = position_jitter(width = 0, height = if (refined) 0 else JIT)
  g = ggplot(D, aes(xx, yy, fill = .fill)) +
    ## one alpha and one size for every panel: patchwork compares guides
    ## including the layer parameters, so panels drawn differently collect twice
    geom_point(shape = 21, colour = "grey15", stroke = 0.35,
               alpha = if (identity_line) 0.85 else 0.6,
               size = if (identity_line) SIZE_A else SIZE_U,
               ## panel A contributes the guide too now that it carries the
               ## same two-level fill as the rate panels; with the three-level
               ## extinction fill it did not, and collecting both produced two
               ## legends, which is why it was suppressed there
               ## ONE LEGEND FOR THE FIGURE. patchwork collects guides by
               ## comparing them, and a point key from this panel is not the
               ## same object as a line key from the ECDF panels, so both would
               ## be collected and the figure would carry the legend twice.
               ## Panel A yields its guide when B to D are drawing ECDFs: the
               ## two classes and their colours are identical, so one key set
               ## serves all four.
               position = jit,
               show.legend = !(identity_line && (CRIT_EXT || USE_ECDF))) +
    ## drop = FALSE so every rate panel contributes a guide with the SAME two
    ## keys even where one class is absent from that panel; a guide with a
    ## different key set would be collected as a second legend.
    scale_fill_manual(values = CRIT_COL, name = NULL, drop = FALSE) +
    ## the key is fully overridden so every rate panel contributes an identical
    ## guide and the collected legend appears once
    guides(fill = guide_legend(order = 1, override.aes = list(size = 3.2, alpha = 1))) +
    labs(x = if (identity_line) expression("Time of the " * P[LP] ~ "dip, Tregs-OFF arm")
             else expression(P[LP] ~ "at the dip, Tregs-OFF arm"),
         y = if (refined && identity_line) "Delay"
             else if (refined) paste(TRT[[tt]]$name, "(bisected)")
             else TRT[[tt]]$name) + th
  ## PSEUDO-LOG, NOT LOG, ON THE RATE PANELS. P_LP at the dip is exactly zero for
  ## every combination that reaches the extinction floor, and log10(0) is -Inf, so
  ## a log axis would silently DROP precisely the cases the floor was added to
  ## represent. pseudo_log is linear near zero and logarithmic above it, so the
  ## zeros are drawn, at zero, and the small values are still spread out. sigma
  ## sets where it turns over; 0.05 is below the smallest nonzero dip observed.
  g = g + (if (identity_line) scale_x_log10()
           else scale_x_continuous(
                  trans = scales::pseudo_log_trans(sigma = 0.05, base = 10),
                  breaks = c(0, 0.1, 1, 10, 30)))
  g = g + (if (refined) scale_y_log10()
           else scale_y_continuous(breaks = seq_along(lv),
                                   labels = format(lv, trim = TRUE),
                                   limits = c(0.5, length(lv) + 0.5)))
  ## ---- THE FIT, ON THE PANEL AND IN THE TABLE ------------------------------
  ## Panel A is the quantitative result of this figure: the delay at which help
  ## begins against the time the Tregs OFF arm clears the pathogen. Eyeballing a
  ## line is not measuring it, so the panel carries the numbers.
  ##
  ## THREE STATISTICS, EACH ANSWERING A DIFFERENT QUESTION. Spearman says the two
  ## are monotonically related and is robust to the range spanning a decade. The
  ## RATIO thr / t_dip says whether the relation is the identity rather than
  ## merely increasing -- a median near 1 is the claim; its IQR is the spread, and
  ## it is reported instead of a mean because the sample is small. The log-log
  ## SLOPE says whether it is proportional across the range rather than only near
  ## the middle: a slope of 1 with a CI containing 1 is proportionality, a slope
  ## reliably away from 1 would mean the delay scales as a power of the dip time
  ## and the identity line is the wrong reference.
  ##
  ## Computed on the REFINED values only. The coarse grid quantises upward -- the
  ## smallest sufficient level on a grid of 0,10,...,960,2400 always overshoots --
  ## so a ratio taken from it would be biased above 1 by construction and would say
  ## nothing about the system. Where DSTAR is absent the panel carries no fit.
  if (identity_line && refined && nrow(D) >= 3L) {
    .r  = D$thr / D$t_dip
    .sp = suppressWarnings(cor(D$thr, D$t_dip, method = "spearman"))
    .fm = try(lm(log10(thr) ~ log10(t_dip), data = D[thr > 0 & t_dip > 0]),
              silent = TRUE)
    .sl = if (inherits(.fm, "try-error")) NA_real_ else unname(coef(.fm)[2])
    .ci = if (inherits(.fm, "try-error")) c(NA_real_, NA_real_) else
            suppressWarnings(confint(.fm)[2, ])
    DELAY_FIT <<- data.table(
      n = nrow(D), spearman = round(.sp, 3),
      ratio_med = round(median(.r), 3),
      ratio_q25 = round(quantile(.r, 0.25), 3),
      ratio_q75 = round(quantile(.r, 0.75), 3),
      loglog_slope = round(.sl, 3),
      slope_lo = round(.ci[1], 3), slope_hi = round(.ci[2], 3),
      pct_off_extinct = round(100 * mean(D$extinct), 2))
    cat("\n[Supp. Fig. 3] panel A fit, delay at which help begins vs time of the dip:\n")
    print(DELAY_FIT)
    ## Placed in the panel's own data space rather than at a fixed device
    ## position, so it stays put if the axis limits move. size is in MILLIMETRES
    ## here and points in theme() -- they are not the same unit.
    ## THE EXTINCT FRACTION IS ON THE PANEL, NOT ONLY IN THE CONSOLE. The claim
    ## the panel supports is that a delay converts harm to help exactly where
    ## the untreated arm clears the pathogen, so the share of plotted points
    ## whose Tregs OFF arm reached the extinction floor is part of the result
    ## and belongs beside the fit. It is expected to be 100%, and printing it
    ## is what makes that an observation rather than an assumption: if a future
    ## dataset puts a non-extinct point on this panel, the number says so.
    .ex = 100 * mean(D$extinct)
    .lab = sprintf("n = %d\nTregs-OFF clears: %.1f%%\nSpearman %.2f\nratio %.2f (%.2f-%.2f)\nslope %.2f (%.2f-%.2f)",
                   nrow(D), .ex, .sp, median(.r), quantile(.r, 0.25), quantile(.r, 0.75),
                   .sl, .ci[1], .ci[2])
    cat(sprintf("[Supp. Fig. 3] panel B: %.1f%% of the %d plotted points have a Tregs-OFF arm\n     that reaches the extinction floor\n",
                .ex, nrow(D)))
    ## SIZE IS IN MILLIMETRES HERE AND POINTS IN theme(), and the two are not the
    ## same unit -- .pt is ggplot2's 72.27/25.4 conversion. The old literal factor
    ## put this at about 6pt against a 19pt theme, which is why it was unreadable.
    ## Expressed against BASE_SZ it now tracks the rest of the figure's text.
    ## Lower right: the points follow the diagonal, and at this type size the block
    ## reached them from the upper left.
    g = g + annotate("text", x = Inf, y = min(D$thr),   # y finite: log10(-Inf) is NaN
                     label = .lab, hjust = 1.01, vjust = 0,
                     size = BASE_SZ * 0.8 / .pt,
                     colour = "grey20", lineheight = 1.1)
  }

  ## only where both axes are times: the delay equal to the dip time
  if (identity_line) {
    if (refined) {
      g = g + geom_abline(slope = 1, intercept = 0, colour = "grey45",
                          linetype = "22", linewidth = 0.7)
    } else {
      ID = data.table(xx = lv[lv > 0]); ID[, yy := match(xx, lv)]
      g = g + geom_line(data = ID, aes(xx, yy), inherit.aes = FALSE,
                        colour = "grey45", linetype = "22", linewidth = 0.7)
    }
  }
  g
}
## n_extinct is carried here because the rate panels no longer colour extinction:
## the figure shows the neutral/help split, the table says how much of each class
## eliminated the pathogen, and no information is lost by the change.
TAB = L[, .(n = .N, n_extinct = sum(extinct), least_step = min(step),
            most_step = max(step), median_step = as.numeric(median(step)),
            spearman_vs_tdip = round(suppressWarnings(cor(step, t_dip, method = "spearman")), 3)),
        by = .(treat, crit)][order(treat, crit)]
dir.create(dirname(OUT_TAB), showWarnings = FALSE, recursive = TRUE); fwrite(TAB, OUT_TAB)
cat(sprintf("\n[Supp. Fig. 3] wrote %s\n", OUT_TAB))

## ============================================================================
## B, C AND D: HOW MUCH INTERVENTION SUFFICES
## ----------------------------------------------------------------------------
## WHAT WAS WRONG WITH THE SCATTER. Those panels put the load at the dip on x
## and the least sufficient level on y, several hundred points per treatment,
## two classes overlaid on seven jittered rows. The x variable carries almost no
## signal: the Spearman of threshold against load at the dip is 0.25, -0.02,
## 0.25, -0.07, 0.16 and -0.31 over the six treatment-by-class combinations,
## (values from an earlier dataset, not re-checked) against 0.94 for the delay panel. So the arrangement was several hundred marks
## ordered by something that does not order them.
##
## WHAT REPLACES IT. The distribution of the threshold itself, as an empirical
## cumulative fraction, one step line per class. Read it as: by the time the
## intervention has reached this strength, this fraction of the harmed
## population has stopped being harmed (purple) or has reached help (green).
##
## THE DENOMINATOR IS THE POPULATION AT RISK, NOT EACH CLASS SEPARATELY. Both
## curves are divided by the number of combinations harmed at the untreated
## level of that treatment, so neither reaches 1 and the height each plateaus at
## IS the fraction that ever reaches that class. Normalising each curve to its
## own total would hide exactly that. The vertical gap between the two curves at
## any strength is the fraction sitting in neutral and not beyond.
##
## THE X AXIS IS FOLD CHANGE IN THE DIRECTION OF MORE INTERVENTION. The three
## ladders do not run the same way -- intervention increases as the tau_t_d
## multiplier rises and as the tau_10 and tau_t_a multipliers fall -- so
## plotting the raw multiplier would have two panels reading right to left and
## one left to right. Each is converted to a fold change away from the untreated
## arm: the multiplier for tau_t_d, its reciprocal for the other two. Every
## panel then starts at 1, increases rightward, and the three are comparable.
ecdf_panel = function(tt) {
  D = DSTAR[treat == tt & is.finite(thr)]
  ## the population at risk: harmed at the untreated level of THIS treatment
  n0 = uniqueN(SM[treat == tt & level == TRT[[tt]]$lv[1] & dcomp > TOL, ..KK])
  if (!nrow(D) || !n0) return(ggplot() + labs(title = tt) + th)
  D[, fold := if (TRT[[tt]]$sgn < 0) 1 / thr else thr]
  D[, crit := factor(fifelse(bound == "neutral", "Harm to neutral", "Harm to help"),
                     levels = c("Harm to neutral", "Harm to help"))]
  setorder(D, crit, fold)
  E = D[, .(fold = fold, frac = seq_len(.N) / n0), by = crit]
  ## the curves start at the untreated arm, where nothing has converted yet, so
  ## the first step is drawn from fold = 1 rather than from the first threshold
  E = rbind(unique(E[, .(crit)])[, .(crit, fold = 1, frac = 0)], E)
  setorder(E, crit, fold)
  ggplot(E, aes(fold, frac, colour = crit)) +
    geom_step(linewidth = 1.1 * LW_MULT, direction = "hv") +
    scale_colour_manual(values = CRIT_COL, name = NULL, drop = FALSE) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = function(z) paste0(round(100 * z), "%")) +
    labs(x = TRT[[tt]]$fold_lab, y = "Combinations converted") + th
}
## the x label for each, naming the parameter and the direction
TRT$strength$fold_lab = expression("Fold reduction in " * tau[10])
TRT$deact$fold_lab    = expression("Fold increase in " * tau[td])
TRT$induce$fold_lab   = expression("Fold reduction in " * tau[ta])
USE_ECDF = USE_ECDF && !is.null(DSTAR)
rate_panel = function(tt) if (USE_ECDF) ecdf_panel(tt) else panel(tt)
if (USE_ECDF) {
  cat("\n[Supp. Fig. 3] B, C and D show the cumulative fraction converted against fold change;\n")
  cat("     ECDF=0 restores the load-at-dip scatter. Plateau heights:\n")
  print(DSTAR[is.finite(thr) & treat != "delay",
              .(n = .N), by = .(treat, bound)][order(treat, bound)])
}

## ============================================================================
## THE TWO-ROW MAIN FIGURE
## ----------------------------------------------------------------------------
## ROW 1, three panels, the population result.
##   A  what Tregs do to the pathogen: the median P_LP trajectory of the
##      combinations they help against the ones they harm, each under Tregs OFF
##      and Tregs ON. Drawn from pipeline/51_harm_explore_sample.R's cache, because its
##      contrast needs the HELPED combinations and harm_rescue_data.rds holds only the
##      harmed ones.
##   B  the delay at which harm ends against the time the untreated arm
##      eliminates the pathogen, on the identity.
##   C  the three rate treatments in ONE panel: the cumulative fraction of the
##      harmed population converted, against fold change away from the untreated
##      arm. Colour is the class reached, linetype is which rate was changed.
##
## ROW 2, five panels, one combination in detail. The untreated Tregs OFF arm in
## black, and four treated arms at the LEAST intervention that reaches help --
## the shortest delay, the largest tau_10 multiplier, the smallest tau_t_d
## multiplier, the largest tau_t_a multiplier. Five states in the causal order
## the treatments act through: Tr -> M1 -> P_LP -> E_h. S is omitted; the reason
## is at the VARS assignment inside demo_row().
##
## ONE COMBINATION FOR ALL FIVE CURVES, WHICH IS WHY THE CASE IS CHOSEN AND NOT
## ASSUMED. A row with a different baseline per treatment would make the black
## curve ambiguous and the four arms incomparable. The ranking below finds
## the combinations every treatment carries to help -- 178 of the 1,745 harmed in the current
## sweep, all of which clear the pathogen untreated -- and ranks them by how
## typical their elimination time t_dip is (|log10(t_dip / median)|).
## DEMO_ROW picks from that ranking (default 1).
##
## THE FOUR TREATMENT COLOURS. The first version of this row used Okabe-Ito,
## which passed the validator but only just: its worst all-pairs separations
## were 9.6 CVD and 15.6 normal vision, and the two failing pairs were the ones
## that matter here -- vermillion against amber, and pink against blue. At a
## line width of 0.9 with five curves overlaid, those read as one colour. The
## set below was chosen by searching a pool of saturated candidates for the
## 4-subset with the largest worst-pair separation, subject to every colour also
## being separable from the two outcome-class colours used in row 1:
##   lightness band   PASS, all four inside L 0.43-0.77
##   chroma floor     PASS, all four >= 0.1
##   CVD separation   worst all-pairs 13.2 deutan, 8.5 protan
##   normal vision    worst all-pairs 23.1, against 15.6 for Okabe-Ito
##   contrast         WARN, the orange and the teal sit below 3:1 against the
##                    page. The relief the warning asks for is present twice
##                    over: every arm carries a linetype as well as a colour,
##                    and the legend is labelled, so identity is never colour
##                    alone.
## That is a 1.5-fold improvement in the normal-vision worst pair, which is the
## separation the eye actually uses here.
##
## CROSS-FIGURE CONSTRAINT. Row 1 uses green for "harm to help" and purple for
## "harm to neutral", so a treatment curve in either of those hues would read as
## an outcome class. Each of the four was checked pairwise against both: the
## closest is cyan against the class green at 15.8 normal-vision, and every
## other pair is above 17.
##
## THE FOUR COLOURS COME FROM TWO RAMPS, a pair of treatments from each:
## delay and tau_t_d from
## mako (deep blue, teal), tau_10 and tau_t_a from plasma (orange, magenta).
DEMO_ROW = as.integer(Sys.getenv("DEMO_ROW", "1"))
EXPL_CACHE = Sys.getenv("EXPL_CACHE", "./data/harm_explore_cache.rds")
DEMO_TMAX  = as.numeric(Sys.getenv("DEMO_TMAX", if (!is.null(STAMP0$TMAX)) as.character(STAMP0$TMAX) else "25000"))
PLP_DISP   = as.numeric(Sys.getenv("PLP_DISP", "1"))   # panel A only; row 2 is unfloored

TRT_COL = c(delay = "#06D6A0", strength = "#FF467A", deact = "#9564DD", induce = "#2BBBD7")
TRT_LTY = c(delay = "solid", strength = "8242", deact = "42", induce = "11")
## THE ARM LABELS ARE ONE EXPRESSION VECTOR IN FACTOR-LEVEL ORDER. Mixing a
## plain string with a list of expressions -- c("Tregs OFF", TRT_LAB[...]) --
## does not give ggplot a usable labels argument: the first entry is dropped and
## the legend falls back to the bare level names, which is how this first drew
## "off" as a key. c() over expression() calls gives a length-5 expression
## vector, which is what the scale wants.
ARM_LEV = c("off", "delay", "strength", "deact", "induce")
ARM_LAB = c(expression("Tregs-OFF"), expression("Delay"),
            expression(tau[10] ~ "reduced"), expression(tau[td] ~ "increased"),
            expression(tau[ta] ~ "reduced"))
BASE_COL = "#1C1C1C"

## ---- ROW 1 PANEL A: the explore trajectories -------------------------------
explore_panel = function() {
  if (!file.exists(EXPL_CACHE)) {
    cat(sprintf("[Supp. Fig. 3] %s absent -- panel A is empty. Run pipeline/51_harm_explore_sample.R first.\n",
                EXPL_CACHE))
    return(ggplot() + labs(x = "Time", y = expression(P[LP])) + th)
  }
  D = as.data.table(readRDS(EXPL_CACHE)$tr)
  D[, grp  := factor(grp, levels = c("help", "harm"), labels = c("Help", "Harm"))]
  D[, cond := factor(ifelse(theta == 0L, "Tregs-OFF", "Tregs-ON"),
                     levels = c("Tregs-OFF", "Tregs-ON"))]
  D[, pd := pmax(PLP, PLP_DISP)]
  B = D[is.finite(pd), .(med = median(pd), lo = quantile(pd, .25), hi = quantile(pd, .75)),
        by = .(grp, cond, time)][time > 0]
  ggplot(B, aes(time, med, colour = grp, fill = grp, linetype = cond)) +
    geom_ribbon(data = B[cond == "Tregs-OFF"], aes(ymin = lo, ymax = hi),
                colour = NA, alpha = 0.10, show.legend = FALSE) +
    geom_ribbon(data = B[cond == "Tregs-ON"], aes(ymin = lo, ymax = hi),
                colour = NA, alpha = 0.16, show.legend = FALSE) +
    geom_line(linewidth = 0.9) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c(Help = "#1B2CC1", Harm = "#DF301C")) +
    scale_fill_manual(values   = c(Help = "#1B2CC1", Harm = "#DF301C")) +
    scale_linetype_manual(values = c(`Tregs-OFF` = "22", `Tregs-ON` = "solid")) +
    labs(x = "Time", y = bquote(P[LP] ~ " (floored at " * .(formatC(PLP_DISP, format = "g")) * ")")) +
    th + theme(legend.position = "top", legend.title = element_blank(),
               legend.box = "horizontal") +
    guides(colour = guide_legend(ncol = 1, order = 1),
           fill   = guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2))
}

## ---- ROW 1 PANEL C: the three rate treatments in one panel -----------------
## COLOUR IS THE CLASS AND LINETYPE IS THE RATE, so the neutral-to-help gap --
## which is the point of the panel -- is carried by the channel the eye reads
## first, and the three rates are separable within each class. The x axis is
## fold change away from the untreated arm in the direction of more
## intervention, so the three ladders are comparable: the multiplier for
## tau_t_d, its reciprocal for the other two.
rates_panel = function() {
  tts = intersect(c("strength", "deact", "induce"), unique(DSTAR$treat))
  E = rbindlist(lapply(tts, function(tt) {
    D = DSTAR[treat == tt & is.finite(thr)]
    n0 = uniqueN(SM[treat == tt & level == TRT[[tt]]$lv[1] & dcomp > TOL, ..KK])
    if (!nrow(D) || !n0) return(NULL)
    D = copy(D)[, fold := if (TRT[[tt]]$sgn < 0) 1 / thr else thr]
    D[, crit := factor(fifelse(bound == "neutral", "Harm to neutral", "Harm to help"),
                       levels = c("Harm to neutral", "Harm to help"))]
    setorder(D, crit, fold)
    z = D[, .(fold = fold, frac = seq_len(.N) / n0), by = crit]
    rbind(unique(z[, .(crit)])[, .(crit, fold = 1, frac = 0)], z)[, treat := tt][]
  }), fill = TRUE)
  if (!nrow(E)) return(ggplot() + th)
  E[, treat := factor(treat, levels = c("strength", "deact", "induce"))]
  setorder(E, treat, crit, fold)
  ## THE LINETYPES HAVE TO BE UNMISTAKABLE, and the first attempt's were not:
  ## "22" and "42" differ only in dash length and at this line width they read as
  ## the same pattern, so the legend could not be matched to a curve. That
  ## matters here more than in a panel with two series, because two of the three
  ## rates land close together -- tau_10 and tau_ta both run to about 20-fold and
  ## converge within a few percent -- so a reader who cannot separate the
  ## patterns sees one curve where there are two. Solid, fine dotted and
  ## long-dash-dot are distinguishable at a glance and at print size.
  LT = c(strength = "solid", deact = "8242", induce = "11")
  LAB = c(expression(tau[10] ~ "reduced"), expression(tau[td] ~ "increased"),
          expression(tau[ta] ~ "reduced"))
  ## MARKERS CARRY THE TWO CURVES A DASH PATTERN CANNOT SEPARATE. tau_t_d and
  ## tau_t_a act on the same quantity -- the Treg set point depends on the ratio
  ## tau_t_a / tau_t_d and on nothing else -- so their conversion curves lie
  ## within a percent of each other over most of the range and no linetype makes
  ## them distinguishable where they coincide. A marker every few steps does:
  ## the two symbols alternate along the shared path instead of one pattern
  ## hiding the other. tau_10 is far enough from both to need none, so it has
  ## shape NA and is drawn as a bare line.
  ##
  ## THE MARKERS ARE THINNED, not drawn at every step. The curves carry several
  ## hundred steps each; a marker on every one is a thick band, not a symbol.
  ## NMK per curve is enough to identify a line anywhere the reader picks it up.
  ##
  ## AND THEY ARE IN ANTIPHASE. Placing them at the same positions on both
  ## curves defeats the purpose: where tau_t_d and tau_t_a coincide -- which is
  ## most of the range -- the circle and the triangle land on the same pixel and
  ## one hides the other, so the panel still shows a single line. The positions
  ## are therefore taken from a common log-spaced grid of step size stp, with
  ## tau_t_d at phase 0 and tau_t_a at phase 0.5, so along a shared path the two
  ## symbols alternate and both are visible.
  ##
  ## THE POSITIONS ARE FOLD VALUES, NOT ROW INDICES. The two curves have
  ## different numbers of steps and different fold ranges, so an index-based
  ## grid gives them different spacings and the phase offset means nothing. A
  ## grid in log(fold) is the same grid for both, and frac is read off each
  ## curve at those folds by constant interpolation, which puts the marker
  ## exactly on the step it belongs to rather than near it.
  NMK = as.integer(Sys.getenv("N_MARK", "9"))
  PH  = c(deact = 0, induce = 0.5)
  mk_one = function(d, ph) {
    d = d[order(fold)]
    f = d$fold[d$fold > 1]
    if (length(f) < 2L) return(NULL)
    lo = log10(min(f)); hi = log10(max(f))
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) return(NULL)
    stp = (hi - lo) / NMK
    at  = 10^(lo + stp * (seq_len(NMK) - 0.5 + ph))
    at  = at[at >= 10^lo & at <= 10^hi]
    if (!length(at)) return(NULL)
    data.table(fold = at,
               frac = approx(d$fold, d$frac, xout = at, method = "constant",
                             f = 0, rule = 2)$y)
  }
  MK = E[as.character(treat) %in% names(PH),
         { z = mk_one(.SD, PH[[as.character(.BY$treat)]])
           if (is.null(z)) NULL else z }, by = .(treat, crit)]
  SH = c(strength = NA, deact = 16, induce = 17)
  ggplot(E, aes(fold, frac, colour = crit, linetype = treat)) +
    geom_step(linewidth = 1.0 * LW_MULT, direction = "hv") +
    ## show.legend is per aesthetic: this layer contributes the SHAPE key and
    ## nothing else, so the two class keys stay bare lines instead of picking up
    ## a marker they do not mean.
    geom_point(data = MK, aes(shape = treat), size = SIZE_U * 2.3, stroke = 0,
               show.legend = c(shape = TRUE, colour = FALSE, linetype = FALSE)) +
    scale_colour_manual(values = CRIT_COL, name = NULL, drop = FALSE) +
    ## THE LINETYPE AND SHAPE GUIDES MERGE INTO ONE. ggplot combines two guides
    ## when their title, breaks and labels are identical, so each key shows the
    ## line pattern with its marker on it rather than listing the three rates
    ## twice. limits and drop = FALSE keep tau_10 in the key even though it
    ## contributes no marker.
    scale_linetype_manual(values = LT, name = NULL, labels = LAB,
                          limits = names(LT), drop = FALSE) +
    scale_shape_manual(values = SH, name = NULL, labels = LAB,
                       limits = names(SH), drop = FALSE, na.translate = FALSE) +
    scale_x_log10(expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = function(z) paste0(round(100 * z), "%")) +
    labs(x = "Fold change away from the untreated arm",
         y = "Combinations converted") + th +
    ## a wide key, because a short one renders every dash pattern as a dash: the
    ## legend has to show the difference the curves are drawn with
    theme(legend.position = "right",
          legend.key.width = unit(BASE_SZ * 0.16, "cm")) +
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2),
           shape    = guide_legend(order = 2))
}

## ---- ROW 2: one combination, five states, five arms ------------------------
## It integrates, and this is the only part of the figure that does. Five arms
## of one combination is a few seconds; making the computation stage carry them
## would mean storing a trajectory per treatment per combination for a figure
## that draws one.
demo_row = function() {
  if (!file.exists("./data/skeleton_pathogenic_harm.rds"))
    return(list(msg = "skeleton_pathogenic_harm.rds absent"))
  H = DSTAR[bound == "help" & is.finite(thr)]
  if (!nrow(H)) return(list(msg = "no help thresholds in the data file"))
  W = dcast(H, batch + row_ind + Kg + beta_ros + r_L_P ~ treat, value.var = "thr")
  need = c("delay", "strength", "deact", "induce")
  if (!all(need %in% names(W))) return(list(msg = "not every treatment reaches help anywhere"))
  W = W[complete.cases(W[, ..need])]
  if (!nrow(W)) return(list(msg = "no combination is helped by all four treatments"))
  ## candidates are ranked by |log10(t_dip / median t_dip)| over the
  ## all-four, Tregs-OFF-clearing set; DEMO_ROW = 1 takes the most typical
  OFFq = SM[treat == "off", c(KK, "t_dip", "PL_end", "PLP_end"), with = FALSE]
  .xt = if (!is.null(STAMP0$EXT_TOL)) STAMP0$EXT_TOL else 1e-11
  OFFq[, cleared := PL_end < .xt & PLP_end < .xt]
  W = merge(W, OFFq, by = KK)[cleared == TRUE]
  if (!nrow(W)) return(list(msg = "no all-four candidate clears the pathogen untreated"))
  .m = median(W$t_dip); W[, d := abs(log10(t_dip / .m))]; setorder(W, d)
  if (DEMO_ROW > nrow(W)) return(list(msg = sprintf("DEMO_ROW %d but only %d candidates",
                                                    DEMO_ROW, nrow(W))))
  r = W[DEMO_ROW]
  cat(sprintf("\n[Supp. Fig. 3] row 2 demonstrates batch %d row_ind %d Kg %g beta_ros %g r_L_P %g\n",
              r$batch, r$row_ind, r$Kg, r$beta_ros, r$r_L_P))
  cat(sprintf("     of %d combinations helped by all four; elimination at t = %.4g\n",
              nrow(W), r$t_dip))
  cat(sprintf("     levels: delay %.4g | tau_10 x%.4g | tau_td x%.4g | tau_ta x%.4g\n",
              r$delay, r$strength, r$deact, r$induce))

  ## the model, loaded only for this row.
  ##
  ## BRANCH IS SET IN THE GLOBAL ENVIRONMENT, NOT HERE. source() evaluates in
  ## globalenv() by default, so skeleton_defs.R's
  ##     if (!exists("BRANCH")) BRANCH = Sys.getenv("BRANCH", "sterile")
  ## never sees a variable local to this function, falls through to "sterile",
  ## and sets PL0 = PLP0 = 0 -- the demo arms are then integrated with no
  ## pathogen at all and every curve in the second row is wrong while the figure
  ## still renders. The banner it prints is the only thing that shows it, which
  ## is how this was caught; the assignment below is what prevents it.
  assign("BRANCH", "pathogenic_harm", envir = globalenv())
  suppressMessages(source("./lib/skeleton_defs.R"))
  if (!identical(get("BRANCH", envir = globalenv()), "pathogenic_harm") ||
      !isTRUE(get("PL0", envir = globalenv()) > 0))
    return(list(msg = "skeleton_defs.R did not load the pathogenic branch"))
  SK = as.data.table(readRDS("./data/skeleton_pathogenic_harm.rds"))
  cc = SK[batch == r$batch & row_ind == r$row_ind & Kg == r$Kg &
          abs(beta_ros - r$beta_ros) < 1e-12 & abs(r_L_P - r$r_L_P) < 1e-12]
  if (!nrow(cc)) return(list(msg = "the demo combination is not in the skeleton"))
  cc = cc[1]
  trow = tg[.(r$batch, r$row_ind), nomatch = NULL]
  if (!nrow(trow)) return(list(msg = "the demo tau row is not in the grids"))
  p = c(STRUCTURAL, EFFERO); for (nm in KIN) p[[nm]] = trow[[nm]][1]
  p$Kg = r$Kg; p$beta_ros = r$beta_ros; p$r_L_P = r$r_L_P; p$r_LP_P = r$r_L_P
  EPSd = if (!is.null(STAMP0$EPS)) STAMP0$EPS else 0.01
  EHHI = if (!is.null(STAMP0$EH0_HI)) STAMP0$EH0_HI else 99
  e0 = cc$Eh0_off + EPSd * (EHHI - cc$Eh0_off)
  TTd = c(0, exp(seq(log(0.05), log(DEMO_TMAX), length.out = 400)))
  ig = function(y, tt, pp) tryCatch(as.data.table(
    ode(y, tt, rhs_app, pp, method = "lsodar", rootfun = rootfun_ext,
        events = EVENTS_EXT, atol = ATOL, rtol = RTOL, hmax = HMAX, maxsteps = 1e5)),
    error = function(e) NULL)
  armd = function(pp, th, d = 0) {
    if (th == 0 || d <= 0) return(ig(seed(e0), TTd, c(pp, list(theta_treg = th))))
    t1 = c(TTd[TTd < d], d)
    o1 = ig(seed(e0), t1, c(pp, list(theta_treg = 0))); if (is.null(o1)) return(NULL)
    y1 = setNames(as.numeric(unlist(o1[nrow(o1), ..STATE_ORDER])), STATE_ORDER)
    t2 = c(d, TTd[TTd > d]); if (length(t2) < 2) return(o1)
    o2 = ig(y1, t2, c(pp, list(theta_treg = 1))); if (is.null(o2)) return(NULL)
    rbind(o1, o2[-1])
  }
  mk = list(
    delay    = function() armd(p, 1, r$delay),
    strength = function() { q = p; q$tau_10  = p$tau_10  * r$strength; armd(q, 1, 0) },
    deact    = function() { q = p; q$tau_t_d = p$tau_t_d * r$deact;    armd(q, 1, 0) },
    induce   = function() { q = p; q$tau_t_a = p$tau_t_a * r$induce;   armd(q, 1, 0) })
  A = list(); A[["off"]] = armd(p, 0, 0)
  for (nm in names(mk)) A[[nm]] = mk[[nm]]()
  if (any(vapply(A, is.null, logical(1))))
    return(list(msg = paste("arms failed:", paste(names(A)[vapply(A, is.null, logical(1))],
                                                  collapse = ", "))))
  ## S IS NOT DRAWN. dS = sigma_S*Tr - tau_S*S, so S is a first-order low-pass
  ## filter of Tr with no other input: its panel is the Tr panel with the
  ## corners rounded off and carries no ordering, no crossing and no timing that
  ## the Tr panel does not already carry, and it would cost a fifth of the
  ## row's width. The chain the
  ## row shows is therefore Tr -> M1 -> P_LP -> E_h, with S implicit between the
  ## first two.
  VARS = c("Tr", "M1", "PLP", "Eh")
  TRj = rbindlist(lapply(names(A), function(nm)
    cbind(arm = nm, A[[nm]][, c("time", VARS), with = FALSE])), fill = TRUE)
  TRj = melt(TRj, id.vars = c("arm", "time"), measure.vars = VARS,
             variable.name = "state", value.name = "y")
  TRj = TRj[time > 0]
  ## P_LP IS DRAWN LINEAR AND UNFLOORED, on a log time axis. An earlier version
  ## floored it at PLP_DISP
  ## and put it on a log y axis, which shows the elimination more sharply but
  ## puts this panel alone on a different scale.
  ## PLP_DISP still governs panel A, which is the explore figure and has always
  ## been floored.
  TRj[, arm := factor(arm, levels = ARM_LEV)]
  YLAB = list(Tr = expression(T[r]), S = "S", M1 = expression(M[1]),
              PLP = expression(P[LP]), Eh = expression(E[h]))
  one_state = function(v) {
    D = TRj[state == v]
    g = ggplot(D[arm != "off"], aes(time, y, colour = arm, linetype = arm)) +
      ## the baseline is drawn FIRST and in black so the treated arms sit on top
      ## of it; it is mapped rather than fixed so it appears in the legend
      geom_line(data = D[arm == "off"], linewidth = 1.2 * LW_MULT) +
      geom_line(linewidth = 0.9 * LW_MULT) +
      ## values are indexed BY NAME so the ordering of TRT_COL cannot silently
      ## repaint an arm; labels are positional, hence ARM_LEV as the level order
      scale_colour_manual(values = c(off = BASE_COL, TRT_COL[ARM_LEV[-1]]),
                          breaks = ARM_LEV, labels = ARM_LAB, name = NULL) +
      scale_linetype_manual(values = c(off = "solid", TRT_LTY[ARM_LEV[-1]]),
                            breaks = ARM_LEV, labels = ARM_LAB, name = NULL) +
      scale_x_log10() +
      labs(x = "Time", y = YLAB[[v]]) + th +
      ## THE SAME KEY WIDTH AS PANEL C. Without it the default key is about a
      ## third as wide, and a dash pattern rendered into that space collapses to
      ## a single dash -- the row 2 keys were showing four nearly identical
      ## stubs while the curves they name are clearly different. One constant,
      ## used in both places, so the two rows cannot drift apart.
      theme(legend.position = "right",
            legend.key.width = unit(BASE_SZ * 0.16, "cm"))
    g
  }
  list(panels = lapply(VARS, one_state), key = r)
}

dir.create(dirname(OUT_FIG), showWarnings = FALSE, recursive = TRUE)
## ---- ASSEMBLY --------------------------------------------------------------
## TWO_ROW=0 falls back to the single-row figure: panel A plus the three rate
## panels drawn separately, which is what this script produced before the second
## row existed and which is still the quickest way to look at one treatment.
TWO_ROW = as.integer(Sys.getenv("TWO_ROW", "1")) == 1L
if (TWO_ROW) {
  .dm = demo_row()
  if (!is.null(.dm$msg)) {
    cat(sprintf("\n[Supp. Fig. 3] SECOND ROW SKIPPED: %s\n", .dm$msg))
    cat("     No combination qualifies for the demonstration row.\n")
  }
  ## GUIDES ARE NOT COLLECTED ACROSS THE TWO ROWS. The rows encode different
  ## things with colour -- outcome class above, treatment below -- so one
  ## collected legend would be wrong even though patchwork would happily build
  ## it. Each panel keeps its own, and the row 2 panels are given one legend
  ## between them by collecting within that row only.
  ROW1 = explore_panel() | panel("delay", identity_line = TRUE) | rates_panel()
  if (is.null(.dm$msg)) {
    ROW2 = Reduce(`|`, .dm$panels) + plot_layout(guides = "collect")
    G = ROW1 / ROW2 + plot_layout(heights = c(1, 1))
  } else {
    G = ROW1
  }
  png(OUT_FIG, width = FIG_W2, height = FIG_H2, res = FIG_RES)
  print(G + plot_annotation(tag_levels = "A"))
  invisible(dev.off())
} else {
  png(OUT_FIG, width = FIG_W, height = FIG_H, res = FIG_RES)
  print((panel("delay", identity_line = TRUE) | rate_panel("strength") |
         rate_panel("deact") | rate_panel("induce")) +
        plot_layout(guides = "collect") +
        plot_annotation(tag_levels = "A",
                        theme = theme(legend.position = "right")))
  invisible(dev.off())
}
cat(sprintf("[Supp. Fig. 3] wrote %s\n", OUT_FIG))
## AFTER the panels, not before: DELAY_FIT is filled while panel A is built, so a
## write placed above the print() would always see the NULL it was initialised to.
## The fit goes to its own file rather than into TAB, whose rows are per treatment
## and class and whose columns mean something else.
if (!is.null(DELAY_FIT)) {
  fwrite(DELAY_FIT, OUT_FIT)
  cat(sprintf("[Supp. Fig. 3] wrote %s\n", OUT_FIT))
} else {
  cat("[Supp. Fig. 3] no panel A fit written -- the data file carries no bisected\n")
  cat("     thresholds, so the delay axis is the coarse grid and a ratio taken\n")
  cat("     from it would be biased above 1 by the grid rather than by the system.\n")
}
cat("   FOR THE CAPTION: one point per combination, purple for harm to neutral\n")
cat("   and green for harm to help, the same two classes in all four panels.\n")
cat("   Combinations that never leave harm are not drawn. A puts the delay at\n")
cat("   which harm ends against the time that arm eliminates the pathogen, both\n")
cat("   times, on a continuous axis, and the dashed line is the identity. Every\n")
cat("   combination in A drives the pathogen below the extinction floor; that is\n")
cat("   a property of the data rather than a selection, and it is checked at run\n")
cat("   time. C shows the cumulative fraction of the harmed population that has\n")
cat("   reached each class by a given intervention strength, against fold change\n")
cat("   away from the untreated arm on a log axis, with linetype naming the rate\n")
cat("   that was changed. Both curves of a rate share one denominator, the\n")
cat("   number of combinations harmed untreated under it, so the height each\n")
cat("   settles at is the fraction that ever reaches that class and the vertical\n")
cat("   gap between them is the fraction left in neutral.\n")
cat("   A is the median P_LP of the combinations Tregs help against the ones\n")
cat("   they harm, each under Tregs-OFF and ON, with the interquartile range as\n")
cat("   a band. D to G are ONE combination: the untreated Tregs-OFF arm in black\n")
cat("   and four treated arms at the least intervention that reaches help.\n")
