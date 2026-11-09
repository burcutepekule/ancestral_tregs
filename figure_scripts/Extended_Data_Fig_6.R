## ============================================================================
## Extended_Data_Fig_6.R -- Extended Data Fig. 6: the unaligned marginals of
## P(rescue) along Treg phenotypic stability and along antigenic specificity.
## ----------------------------------------------------------------------------
## These are the raw, per-combination marginals: each combination contributes
## one curve on its own axis, summarised as nine nested central quantile bands
## with the mean in black and the identity as a dot-dashed red line. The aligned
## versions (Fig. 5 A, B) carry the result: the spread of midpoints that these
## panels average over is exactly what blurs the threshold in the stability
## marginal, so the unaligned panels understate the difference between the two
## axes. They are shown because the reader needs to see that the alignment is a
## re-centring of real curves and not a transform that manufactures the step.
## in : ./data/fig05_switch_prob_cells.rds, ./data/fig05_switch_prob_units.rds,
##      ./data/fig05_units_*.rds, ./data/rescue_*.rds
## out: ./figures/extended_data/Extended_Data_Fig_6.png
## RUN: Rscript figure_scripts/Extended_Data_Fig_6.R      (no integration; seconds)
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
source("./lib/fig_imperfect_treg.R")
S = imperfect_treg_shape_data("", trim = as.numeric(Sys.getenv("TRIM", "0.6")))
cat(sprintf("[ED Fig. 6] %d combinations | unaligned marginals\n", length(S$units)))
imperfect_treg_write(imperfect_treg_tworow_raw(S, letters = c("A", "B")),
                     "./figures/extended_data/Extended_Data_Fig_6.png",
                     IMP_FIG_W * 0.499, IMP_FIG_H * 0.33)
