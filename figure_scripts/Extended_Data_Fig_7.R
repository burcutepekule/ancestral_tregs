## ============================================================================
## Extended_Data_Fig_7.R -- Extended Data Fig. 7: where the probability of
## rescue is gained along each axis, what single combinations look like, and
## what the combinations are doing underneath the means of Fig. 5.
##   A, B     change in P(rescue) over one grid step, along stability and eta
##   C-F      four single combinations, both marginals against the axis value,
##            with a rule at each curve's own midpoint; ordered by the stability
##            midpoint, which is the spread the alignment of Fig. 5 removes
## The band composition at each grid value (imperfect_treg_composition() in
## lib/fig_imperfect_treg.R) is the one view that shows the mass sitting at
## exactly 0 and exactly 1, which a quantile fan cannot; it is left out because
## the text gives those fractions as numbers.
## in : ./data/fig05_switch_prob_cells.rds, ./data/fig05_switch_prob_units.rds,
##      ./data/fig05_units_*.rds, ./data/rescue_*.rds
## out: ./figures/extended_data/Extended_Data_Fig_7.png
## RUN: Rscript figure_scripts/Extended_Data_Fig_7.R      (no integration; seconds)
##      EXAMPLES=389,151,203,195 Rscript figure_scripts/Extended_Data_Fig_7.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
source("./lib/fig_imperfect_treg.R")
IDS = as.integer(strsplit(Sys.getenv("EXAMPLES", "389,151,203,195"), ",")[[1]])
P = imperfect_treg_panels("")
S = imperfect_treg_shape_data("")
stopifnot(all(IDS %in% S$units))
## the shared grey scale travels with whichever panel carries it, so it is moved
## off A and onto B to keep it at the right-hand edge of the figure
fig = (P$pC + labs(title = imp_ttl("A", IMP_T$C)) + theme(legend.position = "none") |
       P$pB + labs(title = imp_ttl("B", IMP_T$B)) + theme(legend.position = "right")) /
      imperfect_treg_examples(S, IDS, LETTERS[3:6]) +
      patchwork::plot_layout(heights = c(1.5, 0.85))
imperfect_treg_write(fig, "./figures/extended_data/Extended_Data_Fig_7.png",
                     IMP_FIG_W * 0.84, IMP_FIG_H * 0.92)
