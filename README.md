# Cell-mediated immune regulation may have evolved before adaptive immunity in jawed vertebrates: a modeling study

Code and data for the paper by Burcu Tepekule, Lucas Mathieu, Viktor Müller,
Roger D. Kouyos and C. Jessica E. Metcalf.

The paper studies a 13-equation model of an injured gut epithelium: pathogens, epithelial
cells, macrophages, reactive oxygen species and regulatory T cells (Tregs). It samples the
model's parameters widely, finds the parameter combinations in which the tissue has two
possible fates (it recovers, or it stays inflamed), and asks what Tregs change. This
repository holds the model, the parameter sweep, the selection of parameter combinations,
the data behind every figure, and one script per figure.

**Contents**

1. [Three commands](#1-three-commands)
2. [What the paper says is here, and where to find it](#2-what-the-paper-says-is-here-and-where-to-find-it)
3. [What you need installed](#3-what-you-need-installed)
4. [Map of the repository](#4-map-of-the-repository)
5. [The figures, one by one](#5-the-figures-one-by-one)
6. [The data](#6-the-data)
7. [The pipeline, stage by stage](#7-the-pipeline-stage-by-stage)
8. [The model in brief](#8-the-model-in-brief)
9. [What was checked](#9-what-was-checked)
10. [Names you will meet in the code](#10-names-you-will-meet-in-the-code)
11. [Citation, license, contact](#11-citation-license-contact)

---

## 1. Three commands

Run them from the repository root.

```sh
Rscript setup_data.R       # 1. unpack the data: CSV -> the .rds files R reads   (about 2 min)
bash make_figures.sh       # 2. render every figure of the paper                 (about 9 min on 6 cores)
open figures/              # 3. look at them (macOS; any file browser will do)
```

Step 1 computes nothing; it only converts the shipped CSV files. Step 2 reads those data
and draws Figs. 2 to 5, Extended Data Figs. 1 to 7 and Supplementary Figs. 1 to 3. Fig. 1
(the model diagram) is drawn by hand and is kept as an image.

To render one figure: `Rscript figure_scripts/Fig_2.R`.

## 2. What the paper says is here, and where to find it

The *Data and code availability* statement of the paper lists what this repository
contains. Each item, and its place:

| the paper says | where it is |
|---|---|
| the ODE system | `lib/ODE_SYSTEM_v6.R` (`make_ode_system()`, the right-hand side). The same right-hand side with its constants and initial conditions, as used by the sweep and the figures: `rhs_app` in `pipeline/01_sweep.R` and `lib/skeleton_defs.R`. Solver settings and the pathogen extinction floor: `lib/extinction_floor.R`. |
| the sampling | `pipeline/00_parameter_grids.R` (Latin hypercube over 22 rate constants). The sampled vectors themselves: `data/tau_grid_v6_1.csv` … `_5.csv`. |
| the sweep pipeline | `pipeline/01_sweep.R` (one parameter combination at a time: integrate from a line of initial conditions, find the fixed points, locate the thresholds), `pipeline/02_assemble.R` (merge and check the output), launchers in `pipeline/cluster/`. |
| the fixed-point finder | `lib/FP_finder_v6.R`, with `lib/attractor_cluster.R` (when two fixed points count as one attractor). |
| the selection stages | `pipeline/03_select.R` (the bistability gate, the classes, the branches) and `pipeline/census_core.R` (every count quoted in the Methods). The result, the gated set: `data/selected_big.*.csv.gz`. |
| one script per figure | `figure_scripts/`, named as the figures are named in the paper. |
| … to reproduce every figure | `bash make_figures.sh`; see [section 5](#5-the-figures-one-by-one). |

## 3. What you need installed

- **R** (4.6.1 was used) with these packages:

  ```r
  install.packages(c("data.table", "deSolve", "ggplot2", "patchwork", "scales", "gtable",
                     "ggh4x", "ggnewscale", "ggrepel", "viridisLite", "MASS"))
  ```

  Versions used: `data.table` 1.18.4, `deSolve` 1.42, `ggplot2` 4.0.3, `patchwork` 1.3.2,
  `scales` 1.4.0, `gtable` 0.3.6, `ggh4x` 0.3.1, `ggnewscale` 0.5.2, `ggrepel` 0.9.8,
  `viridisLite` 0.4.3, `MASS` 7.3.65.
- **Python 3 with Pillow** (`pip install pillow`), used only to stack two separately
  rendered rows into one image for Fig. 3, Fig. 5, Extended Data Fig. 2 and
  Supplementary Fig. 2.
- About 1.5 GB of free disk space (the CSV files, plus the `.rds` copies step 1 makes).
- A cluster is needed only to redo the expensive stages of [section 7](#7-the-pipeline-stage-by-stage).
  Everything in section 1 runs on a laptop.

Every script works from the repository root. Either start it there, or set the
environment variable `SRC` to the path of the repository.

## 4. Map of the repository

```
README.md             this file
setup_data.R          step 1: data/*.csv  ->  data/*.rds
make_figures.sh       step 2: all figures

figure_scripts/       one script per figure, named as in the paper
figures/              the rendered figures
  main/                 Fig_1.png … Fig_5.png
  extended_data/        Extended_Data_Fig_1.png … _7.png
  supplementary_information/   Supplementary_Fig_1.png … _3.png
  parts/                Fig_5_row2_timeline.png, the hand-drawn row of Fig. 5
source_data/          small tables with the numbers behind some panels, written by the figure scripts

lib/                  the model and the shared code
  ODE_SYSTEM_v6.R       the ODE system
  FP_finder_v6.R        the fixed-point finder
  attractor_cluster.R   when two fixed points are one attractor
  skeleton_defs.R       constants, initial conditions, the right-hand side used for integration
  extinction_floor.R    solver settings and the pathogen extinction floor
  plp_clean.R           the one rule for a negative pathogen load in sweep output
  corridor_defs.R       the basin-boundary machinery (Supplementary Fig. 1)
  fig_*.R               drawing code shared by figures that come as a family
  data_io.R             reads and writes the CSV files without losing a bit

pipeline/             everything that produces data, numbered in run order (section 7)
  cluster/              Slurm and local launchers for the expensive stages
  checks/               scripts that recompute the numbers quoted in the Methods

data/                 the data as CSV, with its own README.md
docs/                 methods.pdf (the model, the numerics, the sweep, the gate) and
                      figure_methods.pdf (how each figure was made, one section per figure)
```

## 5. The figures, one by one

All scripts are in `figure_scripts/`; all data files are in `data/` (the table names the
`.rds` file the script reads; the CSV of the same name is what is shipped). "Stage" is the
script in `pipeline/` that produced the data.

| figure | what it shows | script | data it reads | stage |
|---|---|---|---|---|
| **Fig. 1** | diagram of the model | none: drawn by hand, `figures/main/Fig_1.png` | — | — |
| **Fig. 2** | efferocytosis knockout: what efferocytosis-based regulation of macrophages buys | `Fig_2.R` | `effero_knockout`, `effero_fg_extra`, `effero_fg_extra2` | `06`, `06b` |
| **Fig. 3** | Tregs shift both attractors and widen the range of injury the tissue recovers from (sterile). A, B: phase portraits. C: load balance. | `Fig_3.R` | `skeleton_sterile`, `rescue_sterile`, `load_balance_sterile`†, `selected_big` | `03`, `20` |
| **Fig. 4** | damage-matched comparison: with Tregs the host can afford a stronger ROS response at equal damage | `Fig_4.R` | `fig04_gated_units`, `fig04_damage_matched_gated`, `fig06_window_damage`, `fig04b_v3_trajectories`†, `selected_big` | `43`, `40`, `44`, `41` |
| **Fig. 5** | probability of rescue against Treg stability and specificity (A, B); evolutionary timeline (C) | `Fig_5.R`, `Fig_5_assemble.py` | `fig05_switch_prob_cells`, `fig05_switch_prob_units`, `fig05_units_*`, `rescue_*`. Panel C is the hand-drawn `figures/parts/Fig_5_row2_timeline.png`. | `60`–`63` |
| **Extended Data Fig. 1** | attractor shifts and rescue gaps, sterile and pathogenic help | `Extended_Data_Fig_1.R` | `skeleton_*`, `rescue_*` (sterile, pathogenic help) | `03` |
| **Extended Data Fig. 2** | as Fig. 3, for the pathogenic help branch | `Extended_Data_Fig_2.R` | `skeleton_pathogenic_help`, `rescue_pathogenic_help`, `load_balance_pathogenic_help`†, `selected_big` | `03`, `20` |
| **Extended Data Fig. 3** | load balance when the Tregs-ON system has a single attractor | `Extended_Data_Fig_3.R` | `load_balance_*`†, `selected_big` | `20` |
| **Extended Data Fig. 4** | what a stronger ROS response buys and costs without Tregs | `Extended_Data_Fig_4.R` | `fig_beta_effero` | `45` |
| **Extended Data Fig. 5** | three classifications of the damage-matched pairs | `Extended_Data_Fig_5.R` | as Fig. 4, and `fig06_first_peak` | `42` |
| **Extended Data Fig. 6** | the unaligned versions of Fig. 5 A, B | `Extended_Data_Fig_6.R` | as Fig. 5 | `60`–`63` |
| **Extended Data Fig. 7** | where the probability of rescue is gained, and single combinations | `Extended_Data_Fig_7.R` | as Fig. 5 | `60`–`63` |
| **Supplementary Fig. 1** | displacement of the basin boundary at t = 0 | `Supplementary_Fig_1.R` | `corridor_q_sterile_round1`, `corridor_q_pathogenic_help`, `cache_bnd_…`, `selected_big` | `21`–`23` |
| **Supplementary Fig. 2** | the pathogenic combinations in which Tregs harm | `Supplementary_Fig_2.R` | `skeleton_pathogenic_harm`, `rescue_pathogenic_harm` | `03` |
| **Supplementary Fig. 3** | Treg treatments on the combinations Tregs harm | `Supplementary_Fig_3.R` | `harm_rescue_data`, `harm_explore_cache`, `skeleton_pathogenic_harm` | `50`, `51` |

† Not shipped: time courses that take minutes to compute and about 230 MB to store.
`make_figures.sh` computes them the first time (stages `20` and `41`), and they come out
bit-for-bit as the originals ([section 9](#9-what-was-checked)).

All figure scripts also read the parameter grids `tau_grid_v6_*`. No figure script
integrates the model beyond a few example trajectories (the phase portraits of Fig. 3,
Extended Data Fig. 2 and Supplementary Fig. 2, panel D of Fig. 4, panel A of
Supplementary Fig. 1, panels D–G of Supplementary Fig. 3).

How each figure was made, down to the population, the selection, the solver settings and
every number in its caption, is in `docs/figure_methods.pdf`, one section per figure.

## 6. The data

**What is in git: the data behind every figure, as CSV.** 39 tables, 432 MB, no file
above 40 MB. This includes the **gated set** (`selected_big`: the 146,825 parameter
combinations that pass the bistability gate, with all their fixed points and thresholds)
and the output of every expensive experiment. `data/README.md` describes each table.

**What is not in git: the raw sweep.** The sweep evaluated 2,642,031 parameter
combinations; its assembled output is 2.1 GB and does not belong in a git repository.
It is deposited at [DOI]. Three stages need it: `03_select.R`, `census_core.R` and
`45_beta_effero_tradeoff.R`. Their outputs are shipped, so nothing in section 1 needs it.
`data/README.md` lists the few other files that are left out, and why.

**Why CSV, and why step 1.** CSV is an open format that any tool reads. But an ordinary
CSV keeps 15 significant digits of each number, and the scripts join tables on
non-integer keys, so a plain export silently changes results. The files here are written
with exact numbers (up to 17 digits where needed), and `setup_data.R` turns them back into
R tables that are identical, bit for bit, to the ones the paper was made from. That was
checked for every table when the files were written (`pipeline/90_export_csv.R`).

If you only want to look at the data, you do not need R: open the CSV files.

## 7. The pipeline, stage by stage

You do not need this section to reproduce the figures. It is the route from nothing to
the data, for anyone who wants to re-run or change a stage. Scripts are in `pipeline/`;
run them from the repository root. Each script's header lists its inputs, outputs and
switches.

**From the shipped data?** *yes*: the stage runs after `setup_data.R` with nothing else.
*raw sweep*: it needs the 2.1 GB sweep output from [DOI].

### The core: sampling, sweep, selection

| stage | what it does | cost | from the shipped data? |
|---|---|---|---|
| `00_parameter_grids.R` | draws the kinetic parameter vectors (5 grids × 2,000) | seconds | not to be re-run: the grids were drawn without a seed, the stored ones are the record, and the script refuses to overwrite them |
| `01_sweep.R` | for each kinetic vector, 441 combinations of `Kg` × `beta_ros` × `r_L_P`: integrates Tregs OFF and Tregs ON from a line of initial conditions, finds fixed points, bisects the thresholds | about 1.7 core-hours per kinetic vector; 6,000 were swept on a Slurm cluster (`cluster/sweep_array.sh`; on one machine `cluster/sweep_local.sh`) | yes |
| `02_assemble.R` | merges the sweep's chunks, applies the convergence checks | minutes | needs the chunks of `01` |
| `03_select.R` | the gate, the classes, the branches → `selected_big`, `skeleton_<branch>`, `rescue_<branch>` | seconds | raw sweep |
| `census_core.R` | prints every count of the Methods | minutes | raw sweep |

To see the sweep run on your machine without the cluster:
`Rscript pipeline/checks/sweep_smoke_test.R` re-runs stages `01` and `02` on two
combinations of the gated set and compares the result with the stored rows (under a minute).

### The experiments behind single figures

| stage | figure | what it does | cost | from the shipped data? |
|---|---|---|---|---|
| `06_effero_knockout.R`, `06b_…` | Fig. 2 | 10,500 units × 8 arms × 21 initial conditions | about 195 core-hours | yes |
| `20_load_balance.R` | Fig. 3C, ED Figs. 2C, 3 | time courses of the load balance | 2 min per branch on 6 cores | yes (`make_figures.sh` runs it) |
| `21_corridor_tzero.R` | Supp. Fig. 1 | scans of the basin boundary along each state | about 15 core-minutes per combination (cluster) | yes |
| `22_corridor_attractors.R`, `23_corridor_position.R` | Supp. Fig. 1 | attractor coordinates, boundary position | seconds | yes |
| `24_numbers_load_balance_boundary.R` | Fig. 3C, ED Figs. 2C, 3, Supp. Fig. 1 | the numbers quoted in the Methods | seconds | yes, after `20` |
| `43_damage_matched_gated_units.R` | Fig. 4 | draws the 75,000 units | seconds | yes |
| `40_damage_matched.R` | Fig. 4, ED Fig. 5 | finds the ROS rate at which Tregs-ON matches the damage of Tregs-OFF | about 3.7 core-minutes per unit; 400 cluster tasks × 8 h | yes |
| `41_damage_matched_trajectories.R` | Fig. 4B, 4C | time courses of 1,000 pairs | 2.5 min on 6 cores | yes (`make_figures.sh` runs it) |
| `42_first_peak.R`, `44_window_damage.R`, `48_gated_fourth_arm.R` | Fig. 4, ED Fig. 5 | first peak, window damage, the fourth arm | minutes | yes |
| `45_beta_effero_tradeoff.R` | ED Fig. 4 | steady states at nine ROS levels without Tregs | minutes | raw sweep |
| `50_harm_rescue_compute.R` | Supp. Fig. 3 | four Treg treatments on 1,771 harmed combinations | 3.6 h on a 32-core node | yes |
| `51_harm_explore_sample.R` | Supp. Fig. 3A | time courses of 250 help and 250 harm combinations | minutes | yes |
| `60_imperfect_sample.R` | Fig. 5, ED Figs. 6, 7 | the stratified draw of units | seconds | yes |
| `61_imperfect_sweep.R` | Fig. 5, ED Figs. 6, 7 | 11 × 11 grid of stability × specificity, ten replicates | about 27,000 core-hours (cluster) | yes |
| `62_imperfect_read.R`, `63_imperfect_rescue_prob.R` | Fig. 5, ED Figs. 6, 7 | merge the runs, compute the probability of rescue | minutes | need the per-run output of `61` |
| `90_export_csv.R` | — | writes `.rds` tables as the CSV files under git, and checks them | minutes | yes |

A few further scripts (`25_corridor_persistence.R`, `45_short_window.R`, `checks/*`)
compute numbers and diagnostics quoted in the Methods.

The cluster launchers in `pipeline/cluster/` are written for Slurm. They take the
repository path from `SRC` and the scratch directory from `SCRATCH_BASE`; set both for
your system, and adapt the `module load` and `conda activate` lines.

### What to re-run after a change

| you changed | re-run |
|---|---|
| the model equations, solver tolerances, extinction floor or attractor clustering | stage `01` and everything after it |
| the gate (`GAP_CUT` in `03_select.R`) | stage `03` and everything after it |
| a setting of one experiment | that experiment's stages and its figure script |

## 8. The model in brief

Thirteen ordinary differential equations for two compartments, the gut lumen and the
lamina propria, separated by the epithelium:

| state | meaning |
|---|---|
| `PL`, `PLP` | pathogens in the lumen and in the lamina propria |
| `CLP` | commensals in the lamina propria |
| `Eh`, `Ea`, `En` | healthy, apoptotic and necrotic epithelial cells |
| `M0`, `M1` | resting and activated macrophages |
| `R` | reactive oxygen species |
| `T0`, `Tr` | resting and active regulatory T cells |
| `D`, `S` | a danger signal and a suppressive signal |

Tregs OFF and Tregs ON are the same equations with Treg induction switched off or on
(`theta_treg` = 0 or 1). Twenty-two rate constants are sampled; three more parameters are
swept on a grid: `Kg` (the dead-cell load at which efferocytosis is half saturated),
`beta_ros` (ROS production rate) and `r_L_P` (pathogen growth rate; 0 gives the sterile
branch). The equations, the constants, the numerical methods, the gate and the classes
are written out in `docs/methods.pdf`.

## 9. What was checked

These checks were run on the repository as it is here, starting from the CSV files alone.

- **The CSV files lose nothing.** Each of the 39 tables, read back from CSV, is identical
  to the table it was written from: same columns, same classes, same bits.
- **The figures come out.** `setup_data.R` then `make_figures.sh` renders all 14 computed
  figures. Five are identical, pixel for pixel, to the images in the manuscript (Figs. 2
  and 5, Extended Data Figs. 5 and 6, Supplementary Fig. 3). The other nine show the same
  content and differ from the manuscript images only by shifts of one or two pixels in
  the placement of lines and text. All 14 are identical, pixel for pixel, to what the
  authors' working scripts produce from the original data files.
- **The numbers come out.** The tables that the scripts write to `source_data/` (the
  counts and statistics behind Fig. 2, the load-balance figures, Supplementary Figs. 1
  and 3) are byte-identical to the ones the Methods were written from.
- **Recomputed inputs equal the originals.** The load-balance series (stage `20`) and the
  6 million trajectory rows of stage `41`, recomputed here, are bit-for-bit identical to
  the files the manuscript figures were drawn from.
- **Stages re-run from the shipped data reproduce their stored output**: the unit draw of
  stage `43`, the sample of stage `51` and the boundary post-processing of stages `22`
  and `23`, all identical.
- **The sweep runs and agrees.** Two combinations of the gated set, re-run through stages
  `01` and `02`, agree with the stored rows in every class label and in all 134 numeric
  columns to better than 10⁻¹⁰ (`pipeline/checks/sweep_smoke_test.R`).

All of this was done on one machine (R 4.6.1, macOS, Apple silicon). On other systems,
expect the numbers to agree to solver precision and the images to differ in font
rendering.

Not re-run for this repository, because of their cost: the full sweep and the cluster
experiments (stages `01`, `06`, `21`, `40`, `50`, `61`). Their stored outputs are the data
of record; `docs/methods.pdf`, section "Software and reproducibility", reports the
spot checks made on each.

## 10. Names you will meet in the code

| in the code | in the paper |
|---|---|
| `beta_ros` | β_r, the ROS production rate (innate aggression) |
| `r_L_P` | r_P, the pathogen growth rate |
| `Kg` | K_g, the efferocytosis half-saturation load |
| `tau_id`, `batch`, `row_ind` | identifiers of a kinetic parameter vector (`tau_id = 10000 × batch + row_ind`) |
| `theta_treg` | the Treg switch θ (0 = Tregs OFF, 1 = Tregs ON) |
| `phi_instability`, `treg_eff` | φ and η, Treg instability and antigenic specificity (Fig. 5) |
| `comp`, "composite", "total burden" | the total burden of a state |
| `mono_B` | the class ON_mono (the Tregs-ON system has one attractor, the low-burden one) |
| `rescue_gap` | the rescue gap: how far Tregs move the threshold of recoverable injury |
| "corridor" | the basin-boundary analysis of Supplementary Fig. 1 |
| "skeleton" | the attractors, saddle and thresholds of a parameter combination |
| "imperfect" (Treg) | the stability and specificity analysis of Fig. 5 |
| "harm" / "help" | parameter combinations in which Tregs raise / lower the burden of the low-burden attractor |

Data file names were fixed early in the project. A few carry a figure number from a
working draft (`fig06_first_peak` belongs to Fig. 4); the tables in section 5 and in
`data/README.md` are the reference, not the file names.

## 11. Citation, license, contact

If you use this code or these data, please cite the paper:

> Tepekule B., Mathieu L., Müller V., Kouyos R. D., Metcalf C. J. E. *Cell-mediated immune
> regulation may have evolved before adaptive immunity in jawed vertebrates: a modeling
> study.* [journal, year, DOI]

Archived version of this repository: [DOI].

License: [to be stated by the authors before release].

Questions: Burcu Tepekule, burcu.tepekule@unibas.ch
