#!/bin/bash
## ============================================================================
## make_figures.sh -- STEP 2: render every computed figure of the paper.
##
##   Rscript setup_data.R       # step 1, once: CSV -> .rds  (about 2 min)
##   bash make_figures.sh       # step 2: all figures        (about 9 min the first time, 3 min after)
##   FIGS="Fig_4 Extended_Data_Fig_5" bash make_figures.sh   # some of them
##
## Two inputs are time courses that are cheap to compute and large to store, so
## they are not shipped; they are computed here the first time they are needed:
##   data/load_balance_<branch>.rds     pipeline/20_load_balance.R              (Fig. 3C, Extended Data Figs. 2C, 3)
##   data/fig04b_v3_trajectories.rds    pipeline/41_damage_matched_trajectories.R   (Fig. 4B, 4C)
## Everything else is read from the shipped data; the figure scripts integrate
## nothing except a handful of example trajectories.
##
## Fig. 1 (the model diagram) and Fig. 5C (the timeline) are drawn by hand and
## are kept as images: figures/main/Fig_1.png, figures/parts/Fig_5_row2_timeline.png.
##
## out: figures/main/Fig_{2..5}.png
##      figures/extended_data/Extended_Data_Fig_{1..7}.png
##      figures/supplementary_information/Supplementary_Fig_{1..3}.png
##      source_data/*.csv      the numbers behind some panels, written by the figure scripts
##      logs/<figure>.log      what each script printed
## Needs R (packages in README.md) and python3 with Pillow (three figures are
## stacked from two separately rendered rows).
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")"
export SRC="$PWD"
export N_CORES=${N_CORES:-6}
FIGS=${FIGS:-"Fig_2 Fig_3 Fig_4 Fig_5 Extended_Data_Fig_1 Extended_Data_Fig_2 Extended_Data_Fig_3 Extended_Data_Fig_4 Extended_Data_Fig_5 Extended_Data_Fig_6 Extended_Data_Fig_7 Supplementary_Fig_1 Supplementary_Fig_2 Supplementary_Fig_3"}
mkdir -p logs

[ -f data/selected_big.rds ] || { echo "data/selected_big.rds not found -- run 'Rscript setup_data.R' first" >&2; exit 1; }

stage() {  # stage <log name> <command...>
  local log="logs/$1.log"; shift
  local t0; t0=$(date +%s)
  "$@" > "$log" 2>&1 || { echo "FAILED (see $log)"; exit 1; }
  echo "ok ($(( $(date +%s) - t0 )) s)"
}

## ---- the two recomputed inputs ------------------------------------------------
case " $FIGS " in *" Fig_3 "*|*" Extended_Data_Fig_2 "*|*" Extended_Data_Fig_3 "*)
  for B in sterile pathogenic_help; do
    if [ ! -f "data/load_balance_$B.rds" ]; then
      printf "pipeline/20_load_balance.R (%s) ... " "$B"
      stage "20_load_balance_$B" env BRANCH="$B" Rscript pipeline/20_load_balance.R
    fi
  done;;
esac
case " $FIGS " in *" Fig_4 "*)
  if [ ! -f data/fig04b_v3_trajectories.rds ]; then
    printf "pipeline/41_damage_matched_trajectories.R ... "
    stage 41_damage_matched_trajectories Rscript pipeline/41_damage_matched_trajectories.R
  fi;;
esac

## ---- the figures ----------------------------------------------------------------
for f in $FIGS; do
  printf "%s ... " "$f"
  stage "$f" Rscript "figure_scripts/$f.R"
done
echo "figures are in figures/main, figures/extended_data and figures/supplementary_information"
