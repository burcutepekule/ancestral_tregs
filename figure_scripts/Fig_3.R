SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## Fig_3.R -- Fig. 3: the sterile branch, three panels.
##   A, B  phase portraits of two sterile Treg-help combinations (one whose
##         Tregs-ON system is bistable, one where it is monostable)
##   C     the load balance of 20 sterile bistable combinations, the Tregs-OFF
##         and Tregs-ON arms overlaid in one set of axes
## The pathogenic counterpart is Extended_Data_Fig_2.R and has the same shape;
## the ON_mono class is Extended_Data_Fig_3.R, which is panel C's treatment
## alone, twice.
##
## A, B are drawn by lib/fig_skeleton.R (fig_portrait_row) and C by
## lib/fig_load_balance.R with overlay = TRUE. The two halves are rendered
## separately and stacked by figure_scripts/stack_rows.py.
##
## Panel C:
##   - one shared log x axis over the union of the two arms' ranges (about 0.018
##     to 20). On it the Tregs-OFF arm lies wholly inside the shaded overwhelmed
##     region and ends there, while the Tregs-ON arm crosses x = 1 and ends far
##     into the functional region.
##   - line type carries the arm: solid Tregs-ON, dashed Tregs-OFF, as in
##     Fig. 2 and Figs. 3A, 3B.
##   - end marker carries the arm too: triangle Tregs-OFF, circle Tregs-ON.
##   - the start marker is a cross, drawn once per pair (the two arms start from
##     the same state; checked bitwise equal in x and y for all 20 pairs).
##
## Geometry is not set here. The phase half is lib/fig_skeleton.R's
## fig_portrait_row() and the load half is lib/fig_load_balance.R's lb_png();
## both are defined once in those files and shared with the sibling figures, so
## the two halves come out at one pixel width and the panel style is identical
## across Fig. 3C, Extended Data Fig. 2C and Extended Data Figs. 3A, 3B.
## in : ./data/skeleton_sterile.rds, ./data/rescue_sterile.rds,
##      ./data/tau_grid_v6_*.rds, ./data/load_balance_sterile.rds,
##      ./data/load_balance_pathogenic_help.rds, ./data/selected_big.rds
## out: ./figures/main/Fig_3.png
##      ./figures/parts/Fig_3_row{1,2}.png   (the two halves)
## RUN: Rscript figure_scripts/Fig_3.R
## ============================================================================
source("./lib/fig_skeleton.R")
source("./lib/fig_load_balance.R")

OUT  = Sys.getenv("OUT",  "./figures/main/Fig_3.png")
ROW1 = Sys.getenv("ROW1", "./figures/parts/Fig_3_row1.png")
ROW2 = Sys.getenv("ROW2", "./figures/parts/Fig_3_row2.png")
for (p in c(OUT, ROW1, ROW2)) dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)

## ---- A, B: the sterile phase portraits -------------------------------------
## lib/fig_skeleton.R's fig_portrait_row(): its device and its hand-tuned
## anchors are defined once there and shared with Extended_Data_Fig_2.R, so the
## two figures' A, B halves cannot drift apart.
fig_portrait_row(ROW1, "sterile")

## ---- C: the sterile load balance, both arms overlaid ------------------------
g = fig_load_balance("Bistable", branches = "Sterile", letters = "C", overlay = TRUE)
lb_png(ROW2); grid::grid.draw(g); invisible(dev.off()); cat(sprintf("[fig] wrote %s\n", ROW2))

## ---- stack -----------------------------------------------------------------
status = system2("python3", c("figure_scripts/stack_rows.py", ROW1, ROW2, OUT),
                 env = "DPI=250")
if (status != 0) stop("stack_rows.py exited with status ", status)
