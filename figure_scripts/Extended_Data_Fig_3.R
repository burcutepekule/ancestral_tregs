## ============================================================================
## Extended_Data_Fig_3.R -- Extended Data Fig. 3: the load balance during the
## injury, ON_mono class, the counterpart of Fig. 3C and Extended Data Fig. 2C.
##   A  sterile, the Tregs-OFF and Tregs-ON arms overlaid in one set of axes
##   B  pathogenic help, the same
## 20 combinations per branch, each from Eh0 = Eh0_off / 2.
## Drawing code: lib/fig_load_balance.R, overlay = TRUE.
##
## It is the same overlay treatment and the same short label as Fig. 3C and
## Extended Data Fig. 2C, so that the three read as one family. It is the only
## one of the three with no phase-portrait half, so it is load balance alone.
##
## Geometry is lib/fig_load_balance.R's lb_png(), the same device every overlaid
## panel is drawn on, asked for two rows instead of one.
## in : ./data/load_balance_sterile.rds, ./data/load_balance_pathogenic_help.rds,
##      ./data/selected_big.rds
## out: ./figures/extended_data/Extended_Data_Fig_3.png
## RUN: Rscript figure_scripts/Extended_Data_Fig_3.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
source("./lib/fig_load_balance.R")
g = fig_load_balance("ON_mono", letters = c("A", "B"), overlay = TRUE)
OUT = Sys.getenv("OUT", "./figures/extended_data/Extended_Data_Fig_3.png"); dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
lb_png(OUT, n_rows = 2)
grid::grid.draw(g); invisible(dev.off()); cat("wrote", OUT, "\n")
