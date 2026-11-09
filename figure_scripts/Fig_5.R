## ============================================================================
## Fig_5.R -- Fig. 5: how the probability of rescue responds to Treg phenotypic
## stability and to antigenic specificity, and the evolutionary window in which
## a cell with those two capacities could have arisen.
##   A, B  the two marginals of P(rescue), each combination shifted onto the
##         midpoint of its own rise before averaging, drawn where at least TRIM
##         of the curves are defined (the unaligned versions are Extended Data
##         Fig. 6)
##   C     time-calibrated schematic of early vertebrate divergences with the
##         genomic and immunological events that bound the window
##
## Row 1 (A, B) is computed and drawn here. Row 2 (C) is NOT computed: it is a
## schematic drawn by hand in BioRender from literature dates, kept as the image
## figures/parts/Fig_5_row2_timeline.png. The two rows are stacked, and the
## panel letter C typeset to match A and B, by figure_scripts/Fig_5_assemble.py.
## in : ./data/fig05_switch_prob_cells.rds, ./data/fig05_switch_prob_units.rds,
##      ./data/fig05_units_*.rds, ./data/fig05_units.csv, ./data/rescue_*.rds,
##      ./figures/parts/Fig_5_row2_timeline.png
## out: ./figures/parts/Fig_5_row1.png
##      ./figures/main/Fig_5.png
## RUN: Rscript figure_scripts/Fig_5.R      (no integration; seconds)
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
source("./lib/fig_imperfect_treg.R")
S = imperfect_treg_shape_data("", trim = as.numeric(Sys.getenv("TRIM", "0.6")))
cat(sprintf("[Fig. 5 row 1] %d combinations | %s\n", length(S$units),
            paste(S$lims[, sprintf("%s %+.2f to %+.2f", axis, lo, hi)], collapse = " | ")))
imperfect_treg_write(imperfect_treg_tworow_aligned(S, letters = c("A", "B")),
                     "./figures/parts/Fig_5_row1.png",
                     IMP_FIG_W * 0.499, IMP_FIG_H * 0.31)

## ---- stack row 1 over the timeline ------------------------------------------
status = system2("python3", "figure_scripts/Fig_5_assemble.py")
if (status != 0) stop("Fig_5_assemble.py exited with status ", status)
