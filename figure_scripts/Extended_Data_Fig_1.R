SRC = Sys.getenv("SRC", "."); setwd(SRC)
## ============================================================================
## Extended_Data_Fig_1.R -- Extended Data Fig. 1: the attractor shifts and the
## rescue-gap distribution behind the phase portraits: sterile (row 1, the
## population of Figs. 3A, 3B) and pathogenic help (row 2, the population of
## Extended Data Figs. 2A, 2B).
## in : ./data/skeleton_{sterile,pathogenic_help}.rds,
##      ./data/rescue_{sterile,pathogenic_help}.rds, ./data/tau_grid_v6_*.rds
## out: ./figures/extended_data/Extended_Data_Fig_1.png
## RUN: Rscript figure_scripts/Extended_Data_Fig_1.R
## ============================================================================
source("./lib/fig_skeleton.R")
fig_summaries(list(
  list(branch = "sterile",         effect = "help", lab = "Sterile",    tags = c("A","B","C")),
  list(branch = "pathogenic_help", effect = NULL,   lab = "Pathogenic", tags = c("D","E","F"))),
  "./figures/extended_data/Extended_Data_Fig_1.png")
