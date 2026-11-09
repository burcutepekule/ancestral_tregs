SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## Supplementary_Fig_2.R -- Supplementary Fig. 2: the pathogenic combinations in
## which Tregs harm, on one page. Top row, the two phase portraits, drawn by
## fig_portrait_row() as the first row of Fig. 3; bottom row, the attractor
## shifts of each Tregs-ON class and the rescue-gap distribution of the bistable
## class (as Extended Data Fig. 1).
## The two rows are rendered separately and stacked by figure_scripts/stack_rows.py,
## as Fig. 3 is. The top row is Fig. 3's 6400 x 2900 at res 250 taken to this
## page's 9600 width by res 375: the same inches, so every proportion of that row
## carries over (type, squares, brace) at full sharpness.
## in : ./data/skeleton_pathogenic_harm.rds, ./data/rescue_pathogenic_harm.rds,
##      ./data/tau_grid_v6_*.rds
## out: ./figures/supplementary_information/Supplementary_Fig_2.png
##      ./figures/parts/Supplementary_Fig_2_row{1,2}.png   (the two halves)
## RUN: Rscript figure_scripts/Supplementary_Fig_2.R
## ============================================================================
source("./lib/fig_skeleton.R")

OUT  = "./figures/supplementary_information/Supplementary_Fig_2.png"
ROW1 = "./figures/parts/Supplementary_Fig_2_row1.png"
ROW2 = "./figures/parts/Supplementary_Fig_2_row2.png"
for (p in c(OUT, ROW1, ROW2)) dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)

## ---- A, B: the phase portraits ----------------------------------------------
fig_portrait_row(ROW1, "pathogenic_harm", width = 9600, height = 4350, res = 375)

## ---- C, D, E: attractor shifts and rescue gap --------------------------------
png(ROW2, width = 9600, height = 3100, res = 250, pointsize = 38)
layout(matrix(1:3, nrow = 1))
## layout() picks cex from the grid shape; pinned so the page matches the others
par(cex = 0.66)
fig_skeleton("pathogenic_harm", NULL, panels = function(e) {
  par(mar = SUM_MAR)
  e$vio(e$shifts(e$BI, TRUE), e$vcols,
        sprintf("C   Tregs-ON bistable\n(n = %d)", nrow(e$BI)),
        "ON - OFF", gstep = 10, labs = e$labs_bi, ttl_cex = SUM_TTL, lab_cex = SUM_XLAB)
  par(mar = SUM_MAR)
  e$vio(e$shifts(e$MB, FALSE), e$vcols,
        sprintf("D   Tregs-ON monostable\n(n = %d)", nrow(e$MB)),
        "ON - OFF", gstep = 10, labs = e$labs_mb, ttl_cex = SUM_TTL, lab_cex = SUM_XLAB)
  par(mar = SUM_MAR)
  e$gapdens(e$BI, e$LOW, "E   Rescue gap\n(bistable)", xlab = SUM_GAP_XLAB,
            ttl_cex = SUM_TTL, xlab_cex = SUM_XLAB, xaxis_cex = SUM_XAX)
})
invisible(dev.off()); cat(sprintf("[fig] wrote %s\n", ROW2))

## ---- stack -------------------------------------------------------------------
status = system2("python3", c("figure_scripts/stack_rows.py", ROW1, ROW2, OUT), env = "DPI=250")
if (status != 0) stop("stack_rows.py exited with status ", status)
