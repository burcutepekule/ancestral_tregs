## ============================================================================
## pipeline/90_export_csv.R -- write the .rds tables of ./data as the CSV files
## that are kept under git, and CHECK that each one reads back identical.
##
## This is the reverse of setup_data.R. It was run once to produce the CSV files
## of the repository; run it again after regenerating a table with a pipeline
## stage if the new table should replace the shipped CSV.
##
## Every table is written by lib/data_io.R (numbers as exact text, one row per
## column in data/schema.csv), read back, and compared with the .rds it came
## from: same columns, same classes, same bits. The script stops at the first
## table that does not come back identical.
##
## What is NOT exported, and why (see data/README.md):
##   fp_big_batch_<grid>.rds     the assembled sweep, 2.6 million combinations, 2.1 GB
##   fig05_raw_all.rds           every run of the stability-specificity sweep, 209 MB
##   harm_rescue_data.rds $tr    22 million trajectory rows; no figure of the paper reads it
##   load_balance_<branch>.rds   trajectories; pipeline/20_load_balance.R recomputes them in minutes
##   fig04b_v3_trajectories.rds  trajectories; pipeline/41_damage_matched_trajectories.R, minutes
##
## in : <FROM>/<name>.rds                 (FROM defaults to ./data)
## out: ./data/<name>.csv | .csv.gz | .partIofN.csv.gz | <name>__<member>.csv
##      ./data/schema.csv                 column classes, read by setup_data.R
##      ./data/files.csv                  one row per table: rows, columns, size
## RUN: Rscript pipeline/90_export_csv.R
##      FROM=/path/to/rds ONLY=selected_big Rscript pipeline/90_export_csv.R
## ============================================================================
SRC = Sys.getenv("SRC", "."); setwd(SRC)
source("./lib/data_io.R")
FROM = Sys.getenv("FROM", "./data")
ONLY = strsplit(Sys.getenv("ONLY", ""), ",")[[1]]

SHIP = c(
  ## the sampled kinetic parameter vectors (the only record of them: drawn without a seed)
  sprintf("tau_grid_v6_%d", 1:5),
  ## the gated set and its branch files (pipeline/03_select.R)
  "selected_big",
  "skeleton_sterile", "skeleton_pathogenic_help", "skeleton_pathogenic_harm",
  "rescue_sterile", "rescue_pathogenic_help", "rescue_pathogenic_harm",
  ## efferocytosis knockout, Fig. 2 (pipeline/06, 06b)
  "effero_knockout", "effero_knockout_units",
  "effero_fg_extra", "effero_fg_extra_units", "effero_fg_extra2", "effero_fg_extra2_units",
  ## damage-matched comparison, Fig. 4 and Extended Data Fig. 5 (pipeline/40-48)
  "fig04_gated_units", "fig04_damage_matched_gated", "fig06_first_peak", "fig06_window_damage",
  "aux_gated_fourth_arm",
  ## Treg stability and specificity, Fig. 5 and Extended Data Figs. 6, 7 (pipeline/60-63)
  "fig05_units_sterile", "fig05_units_pathogenic_help", "fig05_cells",
  "fig05_switch_prob_cells", "fig05_switch_prob_units",
  ## stronger ROS response without Tregs, Extended Data Fig. 4 (pipeline/45)
  "fig_beta_effero",
  ## basin boundary at t = 0, Supplementary Fig. 1 (pipeline/21-23)
  "corridor_tzero_sterile_round1", "corridor_tzero_pathogenic_help",
  "corridor_lo_off_sterile_round1", "corridor_lo_off_pathogenic_help",
  "corridor_lo_off_pathogenic_help_disagree",
  "corridor_q_sterile_round1", "corridor_q_pathogenic_help",
  "cache_bnd_v6_v6_10621_15_0.0015_0_En_60.46609_76.79843_107.11795_60_1000_40_20_12000_sterile_cheb1",
  ## Treg treatments where Tregs harm, Supplementary Fig. 3 (pipeline/50, 51)
  "harm_rescue_data", "harm_explore_cache")
## list members left out of the export
DROP = list(harm_rescue_data = "tr")
if (length(ONLY)) SHIP = intersect(SHIP, ONLY)

desc = function(x, o, f) {
  one = function(t, member) data.table(object = o, member = member, rows = nrow(t), columns = ncol(t))
  d = if (is.data.frame(x)) one(x, NA_character_) else
    rbindlist(lapply(names(x), function(m) if (is.data.frame(x[[m]])) one(x[[m]], m) else
      data.table(object = o, member = m, rows = 1L, columns = length(x[[m]]))))
  d[, `:=`(csv_files = length(f), csv_MB = round(sum(file.size(f)) / 1e6, 2))][]
}
FILES_CSV = "./data/files.csv"
OUT = if (file.exists(FILES_CSV)) fread(FILES_CSV, na.strings = "")[!object %in% SHIP] else NULL
for (o in SHIP) {
  rds = file.path(FROM, paste0(o, ".rds"))
  if (!file.exists(rds)) stop(rds, " not found")
  x = readRDS(rds)
  for (m in DROP[[o]]) x[[m]] = NULL
  t0 = Sys.time()
  f = write_data_csv(x, file.path("./data", o))
  y = read_data_csv(file.path("./data", o))
  ok = data_identical(x, y)
  cat(sprintf("%-46s %7.1f MB rds -> %6.1f MB csv in %d file(s) | %5.1f s | reads back identical: %s\n",
              substr(o, 1, 46), file.size(rds) / 1e6, sum(file.size(f)) / 1e6, length(f),
              as.numeric(difftime(Sys.time(), t0, units = "secs")), ok))
  if (!ok) stop(o, " does not read back identical -- not exported correctly")
  OUT = rbind(OUT, desc(x, o, f), fill = TRUE)
  fwrite(OUT, FILES_CSV)
  rm(x, y); invisible(gc())
}
cat(sprintf("\n%d tables exported, %.0f MB of CSV in all; every one reads back identical.\n",
            length(SHIP), sum(unique(OUT, by = "object")$csv_MB)))
