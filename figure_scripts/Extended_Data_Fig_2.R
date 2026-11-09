SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## Extended_Data_Fig_2.R -- Extended Data Fig. 2: the pathogenic help branch,
## three panels.
##   A, B  phase portraits of two pathogenic Treg-help combinations, drawn in
##         (E_h, P_LP), the plane of the composite burden
##   C     the load balance during the injury, pathogenic bistable class, the
##         Tregs-OFF and Tregs-ON arms overlaid in one set of axes
## The sterile counterpart is Fig_3.R; this figure is built exactly as that one
## and reads as its twin. The construction, the geometry and the stacking are
## described in the header of Fig_3.R. The only differences are the branch, the
## pathogenic callout anchors and lab_cex.
## in : ./data/skeleton_pathogenic_help.rds, ./data/rescue_pathogenic_help.rds,
##      ./data/tau_grid_v6_*.rds, ./data/load_balance_sterile.rds,
##      ./data/load_balance_pathogenic_help.rds, ./data/selected_big.rds
## out: ./figures/extended_data/Extended_Data_Fig_2.png
##      ./figures/parts/Extended_Data_Fig_2_row{1,2}.png   (the two halves)
## RUN: Rscript figure_scripts/Extended_Data_Fig_2.R
## ============================================================================
source("./lib/fig_skeleton.R")
source("./lib/fig_load_balance.R")

OUT  = Sys.getenv("OUT",  "./figures/extended_data/Extended_Data_Fig_2.png")
ROW1 = Sys.getenv("ROW1", "./figures/parts/Extended_Data_Fig_2_row1.png")
ROW2 = Sys.getenv("ROW2", "./figures/parts/Extended_Data_Fig_2_row2.png")
for (p in c(OUT, ROW1, ROW2)) dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)

## ---- A, B: the pathogenic phase portraits ----------------------------------
fig_portrait_row(ROW1, "pathogenic_help")

## ---- C: the pathogenic load balance, both arms overlaid ---------------------
g = fig_load_balance("Bistable", branches = "Pathogenic", letters = "C", overlay = TRUE)
lb_png(ROW2); grid::grid.draw(g); invisible(dev.off()); cat(sprintf("[fig] wrote %s\n", ROW2))

## ---- stack -----------------------------------------------------------------
status = system2("python3", c("figure_scripts/stack_rows.py", ROW1, ROW2, OUT),
                 env = "DPI=250")
if (status != 0) stop("stack_rows.py exited with status ", status)
