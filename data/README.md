# Data

Everything in this folder is a plain CSV table (`.csv`, or gzip-compressed `.csv.gz`
when it is larger than 8 MB). You can open the files with any tool. The R scripts of
the repository do not read them directly: `Rscript setup_data.R` (run once, from the
repository root) rebuilds one `<name>.rds` per table next to the CSV files, and the
scripts read those.

## How to read the files

| if you see | it means |
|---|---|
| `name.csv`, `name.csv.gz` | one table |
| `name.part1of4.csv.gz`, … | one table cut by rows into parts (GitHub refuses files above 100 MB); stack the parts in order |
| `name__member.csv` | the R object `name` is a list; this file is its member `member` |
| `name__member.list.csv` | a member that is a list of settings: one row, one column per setting |
| a cell like `0.31\|0.42\|97.5` | a list column: several values in one cell, joined by `\|` (only in `selected_big` and the `.list` files) |
| an empty cell | missing (`NA`) |

Three helper files describe the rest:

- `schema.csv`: one row per column of every file, with its R class (`integer`, `numeric`,
  `logical`, `character`, `list:<class>`). `setup_data.R` uses it to restore the tables exactly.
- `files.csv`: one row per table with its number of rows and columns and its size.
- this file.

**Numbers are exact.** A number is written with 15 significant digits when that is enough
to read it back as the same double, and with 17 digits otherwise. This matters because the
scripts join tables on the non-integer keys `beta_ros` and `r_L_P`; ordinary CSV writing
(15 digits) changes those keys in the last bit and silently drops rows from joins. When the
files were made, every table was read back and compared with its original, bit for bit
(`pipeline/90_export_csv.R`); `lib/data_io.R` holds the reader and the writer.

## The key that links the tables

A **kinetic vector** is one draw of the 22 sampled rate constants. It is identified by
`tau_id = 10000 × batch + row_ind`, where `batch` is the parameter grid (1 to 5) and
`row_ind` the row in it.

A **parameter combination** is a kinetic vector together with three swept parameters:
`Kg` (efferocytosis half-saturation load), `beta_ros` (ROS production rate, β_r in the
paper) and `r_L_P` (pathogen growth rate, r_P in the paper; 0 is the sterile branch). Its key is
`(tau_id, Kg, beta_ros, r_L_P)`.

Suffixes used throughout: `_off` and `_on` are the Tregs-OFF and Tregs-ON conditions of
the same combination; `lo`, `hi` and `sad` are the low-burden attractor, the high-burden
attractor and the saddle between them; `Eh0` is the initial epithelial health E_h(0).
State names are those of the model: `PL`, `PLP`, `CLP`, `Eh`, `Ea`, `En`, `M0`, `M1`, `R`,
`T0`, `Tr`, `D`, `S`.

## The tables, one by one

Sizes are of the CSV files. "Stage" is the script in `pipeline/` that wrote the table.

### Parameters

| table | rows | what it is | stage |
|---|---|---|---|
| `tau_grid_v6_1` … `_5` | 2,000 each | the sampled kinetic vectors (Latin hypercube). Columns: the 22 rate constants, the key (`batch`, `row_ind`, `tau_id`), and a few derived diagnostics. **These files are the only record of the parameter sets**: they were drawn without a random seed. Grids 1 to 3 were swept; Fig. 2 draws from all five. | `00` |

### The gated set (the result of the sweep and the selection)

| table | rows | what it is | stage |
|---|---|---|---|
| `selected_big` (4 parts, 156 MB) | 146,825 | **the gated set**: every parameter combination that passes gate criteria 1 to 3 (Tregs-OFF system bistable with attractors more than 20 burden units apart, one threshold, increasing), with all 193 columns of the sweep: thresholds `Eh0_off`, `Eh0_on`; the class of the Tregs-ON system `rescue_mode` (`bistable`, `mono_B` = ON_mono); the effect of Tregs `treg_effect` (`help`, `neutral`, `harm`); the rescue gap `rescue_gap` and its kind `rescue_gap_kind`; all 13 state coordinates of the low attractor, high attractor and saddle in both conditions; leading eigenvalues; solver settings; and every fixed point found (`fp_off_*`, `fp_on_*`, list columns). Criterion 4 of the gate is `rescue_gap_kind` being `measured` or `on_lo_basin_throughout`, which keeps 141,430 rows. | `01`–`03` |
| `skeleton_sterile` | 89,328 | the sterile branch (`r_L_P = 0`) of the gated set, 68 columns: thresholds, classes, total burden of each attractor and saddle, their shifts ON − OFF (`d_*`), eigenvalues (`e_*`), and E_h, M1, P_LP coordinates | `03` |
| `skeleton_pathogenic_help` | 17,495 | the same for pathogenic combinations in which Tregs lower the burden of the low-burden attractor | `03` |
| `skeleton_pathogenic_harm` | 5,246 | the same for pathogenic combinations in which Tregs raise it | `03` |
| `rescue_<branch>` (three files) | as above | an 18-column extract of the same rows: thresholds, directions, class, effect and rescue gap | `03` |

### Fig. 2: efferocytosis knockout

| table | rows | what it is | stage |
|---|---|---|---|
| `effero_knockout` (2 parts) | 1,763,998 | one row per integration: 10,500 units (kinetic vector, `beta_ros`, `r_L_P`) × 8 arms (seven levels of `Kg` with efferocytosis intact, `eff_ko = 1`, and one knockout arm, `eff_ko = 0`) × 21 initial conditions `Eh0`. Columns: the final state (`Eh`, `PLP`, `En`, `M1`, `R`, `Ea`) and the total burden `comp`. | `06` |
| `effero_knockout_units` | 10,500 | the units that were drawn | `06` |
| `effero_fg_extra`, `effero_fg_extra2` | 671,996; 386,400 | the same for the additional pathogenic units used in panels F and G only | `06`, `06b` |
| `effero_fg_extra_units`, `effero_fg_extra2_units` | 4,000; 2,600 | the additional units | `06b` |

### Fig. 4 and Extended Data Fig. 5: damage-matched comparison

| table | rows | what it is | stage |
|---|---|---|---|
| `fig04_gated_units` | 75,000 | the units (combination × initial condition) drawn for the experiment | `43` |
| `fig04_damage_matched_gated` | 52,010 | one row per unit that was run: the search for β_r\* (the ROS rate at which Tregs-ON matches the cumulative damage of Tregs-OFF), its status, and summaries of each arm (areas under pathogen and damage, peaks, clearance times, final states). `i` is the row of the unit in `fig04_gated_units`. | `40` |
| `fig06_first_peak` | 11,940 | the first peak of the pathogen in each arm, per pair | `42` |
| `fig06_window_damage` | 48,842 | cumulative damage of each arm over the matching window, per pair | `44` |
| `aux_gated_fourth_arm` | 1,148 | the fourth arm (Tregs-OFF at β_r\*), used for numbers quoted in the Methods | `48` |

### Fig. 5 and Extended Data Figs. 6, 7: Treg stability and specificity

| table | rows | what it is | stage |
|---|---|---|---|
| `fig05_units_sterile`, `fig05_units_pathogenic_help` | 405; 572 | the 977 units drawn (stratified) for the stability–specificity sweep | `60` |
| `fig05_units.csv`, `fig05_sample_ledger.csv` | 977; — | the same units in one file, and the ledger of the stratified draw. These two were CSV from the start and are read as they are (they are not in `schema.csv`). | `60` |
| `fig05_cells` | 116,886 | one row per (unit, φ, η) cell of the 11 × 11 grid for the 966 retained units: mean and variance of every state, the rescue outcome, and the Tregs-OFF reference | `62` |
| `fig05_switch_prob_cells` | 116,886 | the probability of rescue in each cell (`p_open`, from `n_open` of `n_rep` replicates) | `63` |
| `fig05_switch_prob_units` | 966 | per-unit summaries of the probability surface (range and area along each axis) | `63` |

### Extended Data Fig. 4

| table | rows | what it is | stage |
|---|---|---|---|
| `fig_beta_effero` | 109,899 | 12,211 parameter sets × 9 levels of β_r without Tregs: load, pathogen and epithelium at the two stable steady states, and the weights used for the expected values | `45` |

### Supplementary Fig. 1: basin boundary at t = 0

| table | rows | what it is | stage |
|---|---|---|---|
| `corridor_tzero_sterile_round1`, `corridor_tzero_pathogenic_help` | 118,928; 138,984 | the boundary scans: one row per (combination, state coordinate), with the position of the basin boundary in both conditions | `21` |
| `corridor_lo_off_*` (three files) | — | coordinates of the Tregs-OFF attractors used to place the boundary, and the rows where two sources disagree | `22` |
| `corridor_q_sterile_round1`, `corridor_q_pathogenic_help` | 118,928; 138,984 | the scans with the boundary displacement expressed as a position between the attractors (`q`), which panels B and C plot | `23` |
| `cache_bnd_…_cheb1` | 360 | the boundary of the one example combination of panel A, re-located at 60 times | written by `lib/corridor_defs.R` the first time the figure script ran |

`_round1`: the sterile cluster array ran out of wall time and only its first round was
recovered (see the figure methods).

### Supplementary Fig. 3: Treg treatments where Tregs harm

| table | rows | what it is | stage |
|---|---|---|---|
| `harm_rescue_data__sm` | 74,382 | one row per (combination, treatment, level): 1,771 harmed combinations × 42 arms (Tregs OFF, and four treatments at 8 to 12 levels each), with peaks, dips, final states and the change in total burden | `50` |
| `harm_rescue_data__thresh` | 7,591 | bisected thresholds: the least treatment level that turns harm into neutral or help | `50` |
| `harm_rescue_data__settings.list` | 1 | the solver settings and treatment levels of that run | `50` |
| `harm_explore_cache__tr` | 420,000 | pathogen time courses of 250 help and 250 harm combinations, Tregs OFF and ON (panel A) | `51` |
| `harm_explore_cache__key.list` | 1 | the settings of that sample | `51` |

## What is not in this folder, and why

| file | size | why it is not here | how to get it |
|---|---|---|---|
| `fp_big_batch_1…3.rds` | 2.1 GB | the assembled sweep: all 2,642,031 parameter combinations before the gate. Too large for git. | deposited at [DOI], or re-run stages `01`–`02` (a cluster job) |
| `fig05_raw_all.rds` | 209 MB | every single run of the stability–specificity sweep. `fig05_cells` is its per-cell summary. | deposited at [DOI], or re-run stages `61`–`62` (27,000 core-hours) |
| the `tr` member of `harm_rescue_data` | 1.5 GB in memory | 22 million trajectory rows. No figure of the paper reads it. | re-run stage `50` (about 3.6 hours on a 32-core node) |
| `load_balance_sterile.rds`, `load_balance_pathogenic_help.rds` | 95 MB | time courses, cheap to recompute | `make_figures.sh` runs stage `20` (about 5 min) |
| `fig04b_v3_trajectories.rds` | 135 MB | time courses, cheap to recompute | `make_figures.sh` runs stage `41` (a few minutes) |

## A note on file names

Data file names were fixed early and kept, so that the code that wrote the data and the
code that reads it agree. A few carry a figure number from a working draft that is not the
number in the paper: `fig06_first_peak` and `fig06_window_damage` belong to Fig. 4, and
"corridor" in a file name is the basin-boundary analysis of Supplementary Fig. 1. The
tables above, not the file names, say which figure a table belongs to.
