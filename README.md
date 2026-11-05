# credMortality: reproducibility README

This folder contains the code and intermediate results used to reproduce the code-generated outputs for the manuscript

**“Forecasting Sub-Population Mortality Using Credibility Theory”**

The supplement is organized as an R package together with a `targets` pipeline. The code generates synthetic mortality data from publicly available Human Mortality Database (HMD) life-table information for Italy, then fits the mortality models and produces the figures and saved result files used in the manuscript.

The repeated simulation study uses a reproducible vector of 30 simulation seeds. Each existing simulation seed is also set immediately before every `StMoMo` fit, including the global and separate-population fits. This controls the random initialization used by the underlying `gnm` fitting algorithm. Results may still show small platform- or package-version-dependent numerical differences, so the software versions used for the supplied run are recorded below.
The authors acknowledge these small differences but they remark that they do not impact on the conclusions or the ranking of the models in the manuscript.

## 1. What this supplement reproduces

This supplement reproduces the code-generated manuscript outputs currently referenced in the current submitted manuscript PDF.

In manuscript order, the main code-generated outputs are:

- **Section 4.3.1 / Figure 2**: four-panel credibility illustration on one simulated data set
- **Section 4.3.2 / Figure 3**: nine sample-size plots
- **Section 4.3.3 / Figure 4**: three expected quadratic forecast error fan plots
- **Section 4.4.2 / Figure 6**: six rolling-window out-of-sample performance panels for sub-populations 1–3
- **Appendix B.3 / Figure 7a**: the mortality-trend panel used in the appendix Lexis figure
- **Appendix B.4 / Figures 8 and 9**: appendix figures on comparable scales
- the main saved simulation outputs:
  - `output/results_age_breakets_smooth.csv`
  - `output/lc_msep_results.rds`

## 2. Quick way to inspect the materials without rerunning the pipeline

For many reproducibility checks, the quickest route is to read the
manuscript-aligned `report.Rmd`, or render it to an HTML report that displays
all supplied manuscript figures in order:

```r
rmarkdown::render("report.Rmd")
```

The report explains which pipeline target and saved PDF correspond to each
manuscript figure. For many checks it is sufficient to inspect:

- `report.Rmd` (or the rendered `report.html`)
- the output mapping in Sections 8--10 of this README
- the shipped `figures/` folder
- the shipped `output/` folder
- the shipped `_targets/` cache

A full rerun is only needed when the goal is to recompute everything from source.

## 3. Data source and provenance

The synthetic data used throughout the supplement are generated from female mortality and exposure data for Italy obtained through the Human Mortality Database (HMD). The exact input is versioned in the package as the documented `StMoMoData` object `ITAFemaleData`, created from the Italy snapshot downloaded on **10 February 2026**. It covers calendar years **1872--2022**, which the simulator re-indexes as periods 0--150. After installing the package, inspect it with `?ITAFemaleData` and load it with `data("ITAFemaleData", package = "credMortality")`.

Accordingly:

- the supplement does **not** use confidential data
- the supplement uses **synthetic data generated from HMD-derived mortality schedules**
- reproduction uses only the bundled snapshot and requires neither network access nor HMD login credentials

All 92 bundled files were compared byte-for-byte with the `ITA/` subtree of `hmd_countries_20260210.zip`; there were no missing, extra, or mismatched files. Checksums and the validation result are recorded in `DATA_PROVENANCE.md`. Please see <https://www.mortality.org/> for the source database and its data-use conditions.

## 4. Included files and folders

```text
credMortality/
|-- DESCRIPTION                  R package metadata and dependencies
|-- NAMESPACE                    Exported functions and package imports
|-- LICENSE, LICENSE.md          GPL-3 licensing files
|-- README.md                    Reproduction guide and output mapping
|-- report.Rmd                   Renderable manuscript-aligned figure report
|-- DATA_PROVENANCE.md           Frozen HMD snapshot provenance and checksums
|-- credMortality.Rproj          Optional RStudio project file
|-- _targets.R                   Master reproduction pipeline
|-- R/                           Package source code
|   |-- data.R                   Documentation for the bundled dataset
|   |-- data_generators.R        Synthetic mortality-data generation
|   |-- data_preprocessing.R     Data preparation helpers
|   |-- utils_credibility.R      Model fitting, forecasting, and evaluation
|   `-- zzz.R                    Package startup configuration
|-- data/ITAFemaleData.rda       Frozen Italy female HMD-derived input
|-- data-raw/ITAFemaleData.R     Checked recipe for rebuilding the input
|-- man/                         Generated package help files
|-- scripts/                     Bounded verification and restart scripts
|-- figures/                     34 generated manuscript and appendix PDFs
|-- output/                      Regenerated simulation CSV and MSEP RDS
`-- _targets/                    Included targets cache and run metadata
```

Local session files, credentials, rendered reports, reviewer correspondence,
and internal audit workspaces are intentionally excluded. The generated
manuscript-facing figures and outputs are listed individually in Sections 8
and 11; Section 8.4 also identifies the two retained auxiliary
super-population panels.

## 5. Software requirements

The code was written for R and uses the package-based workflow plus the `targets` pipeline.

At minimum, the following R packages are required by the current code base:

- `targets`
- `data.table`
- `dplyr`
- `tidyr`
- `ggplot2`
- `forecast`
- `fs`
- `StMoMo`
- `rpart`

The package can be installed locally from this folder.

## 6. Data access and helper scripts for quick checks

No manual data download, credential setup, or edit to `_targets.R` is required. Installing the package installs `ITAFemaleData` in its standard `data/` database, and every simulation entry point reads that same versioned object.

`DATA_PROVENANCE.md` records the snapshot checksums. The scripts in `scripts/`
provide bounded verification and a restartable entry point for the repeated
simulation study.

### Included helper scripts

Run each script from the package root:

| Script | Purpose | Command |
|---|---|---|
| `scripts/verify_frozen_hmd.R` | Regenerates the seed-1 synthetic data from the bundled HMD snapshot and compares it with the supplied `dt_0` target. | `Rscript scripts/verify_frozen_hmd.R` |
| `scripts/verify_one_simulation_cell.R` | Reruns seed 61413 (the first simulation-design seed) at rolling-window setting 1 and compares all 180 rows with the corresponding supplied CSV rows. | `Rscript scripts/verify_one_simulation_cell.R` |
| `scripts/run_simulation_study.R` | Installs the current local package and restartably regenerates the repeated simulation branches, summary CSV, BIC table, and rolling-window figures. | `Rscript scripts/run_simulation_study.R` |

`run_simulation_study.R` uses two workers by default. Set the
`SIMULATION_WORKERS` environment variable to a different positive integer when
appropriate for the available memory. The two verification scripts are bounded
spot checks and do not rebuild the full pipeline.

### Regenerate only the repeated simulation study

The restartable script below installs the current local package and rebuilds
only the repeated simulation-study branches, their CSV, and their figures. It
uses the existing `_targets` store, so rerunning it after an interruption skips
completed branches.

```r
# Run from the package root with two parallel workers.
Sys.setenv(SIMULATION_WORKERS = 2)
source("scripts/run_simulation_study.R")
```

From a shell, the equivalent command is
`Rscript scripts/run_simulation_study.R`. Increase `SIMULATION_WORKERS`
conservatively according to available memory. After completion, run
`scripts/verify_one_simulation_cell.R` to compare the first simulation seed
(61413) at rolling-window setting 1 with the supplied CSV.

### Spot-check one simulation replication

The following example reruns the first seed in the simulation design for the
first rolling-window setting, extracts the corresponding rows from the supplied
results, and compares all 180 result rows. Run it from the package root after
installing and loading the local package as described in Section 7.1.

```r
library(credMortality)
library(data.table)

seed_ix <- 61413L # first element of sim_seed_values
years_ix <- 1L

supplied <- fread("output/results_age_breakets_smooth.csv")
supplied <- supplied[
  supplied[["seed_ix"]] == seed_ix & supplied[["years_ix"]] == years_ix,
]

fresh <- run_simulation_study_cell(
  seed_ix = seed_ix,
  stmomoseed = seed_ix,
  years_ix = years_ix,
  ages_fit = 15:85,
  years_fit_basic = 118:140,
  N_groups = 3,
  forecasting_horizon = 1,
  ageb_width = 6,
  exposure_sup = 94500,
  exposure_sub = c(5000, 500),
  global_mortality_model = c("lc", "apc", "rh")
)

setDT(fresh)
setcolorder(fresh, names(supplied))
comparison_key <- c("predictor", "model", "age_breakets")
setorderv(supplied, comparison_key)
setorderv(fresh, comparison_key)

comparison <- all.equal(
  as.data.frame(supplied), as.data.frame(fresh),
  tolerance = 1e-8, check.attributes = FALSE
)
stopifnot(nrow(supplied) == 180L, isTRUE(comparison))
comparison
```

This returns `TRUE`. A small tolerance is used because fitted values can show
minor platform- or package-version-dependent numerical differences.

## 7. Installation and execution

Open R in the project root, i.e. in the folder that contains `DESCRIPTION` and `_targets.R`. This supplement includes the `_targets/` folder from the supplied successful run. It is included intentionally so that users can inspect intermediate objects, perform spot checks, and regenerate selected downstream outputs without rerunning the entire pipeline from scratch.

### 7.1 Install the package and required tooling

```r
install.packages(c("devtools", "targets"))
devtools::install_local(".", upgrade = "never")
library(credMortality)
```

If additional package installation is needed on a fresh machine, install the missing dependencies reported by R.

### 7.2 Inspect the shipped materials without a full rerun

Render the manuscript-aligned report to display the supplied figures and their
target/file mappings in manuscript order:

```r
rmarkdown::render("report.Rmd")
```

Alternatively, read the `report.Rmd` source directly or use the mapping in
Section 8 and inspect the corresponding files in `figures/` or `output/`. None
of these inspection routes rebuilds the analysis. The shipped `_targets/`
cache also allows individual intermediate objects to be inspected without
rebuilding the full pipeline.

### 7.3 Retrieve outputs after running

For file targets, the saved path can be queried by

```r
targets::tar_read(<target_name>)
```

For example:

```r
targets::tar_read(fig2a_pdf)
targets::tar_read(lc_mse_group_2_pdf)
targets::tar_read(lc_msep_results_file)
```

The corresponding files will also appear directly in `figures/` and `output/`.

### 7.4 Clean rerun from scratch

Because the submission folder currently includes a `_targets/` cache, a user who wants a clean rerun should destroy the existing cache first:

```r
library(targets)
tar_destroy(destroy = "all")
```

Then run the full pipeline:

```r
library(targets)
tar_make()
```

This is the main reproduction command for the supplement.

## 8. Mapping from manuscript outputs to pipeline targets

### Important note on sub-population numbering

For some simulation-based figures, the manuscript labels the sub-populations in a different order than the internal file numbering used by the code. In particular:

- file names with `group_2` correspond to **Sub-population 1** in the manuscript
- file names with `group_3` correspond to **Sub-population 2**
- file names with `group_1` correspond to **Sub-population 3**

The same convention applies to the expected quadratic forecast error files:

- `g2` = manuscript Sub-population 1
- `g3` = manuscript Sub-population 2
- `g1` = manuscript Sub-population 3

This README follows the **manuscript numbering and section structure first** and gives the corresponding target/file explicitly.

---

### 8.1 Section 4.3.1 / Figure 2: credibility illustration on one simulated data set

| Manuscript panel | Target | Saved file |
|---|---|---|
| Figure 2a: Estimated sub-population effects | `fig2a_pdf` | `figures/lc_allthetas_zs.pdf` |
| Figure 2b: Estimated variances | `fig2b_pdf` | `figures/lc_variances.pdf` |
| Figure 2c: Credibility weights | `fig2c_pdf` | `figures/lc_simulation_zs.pdf` |
| Figure 2d: Log mortality by age | `fig2d_pdf` | `figures/lc_log_mortality_rates.pdf` |

---

### 8.2 Section 4.3.2 / Figure 3: including more calendar years in the data

| Manuscript panel | Target | Saved file |
|---|---|---|
| Figure 3a: Sub-population 1, age 15 | `additional_analysis_2_15_pdf` | `figures/an215.pdf` |
| Figure 3b: Sub-population 1, age 55 | `additional_analysis_2_55_pdf` | `figures/an255.pdf` |
| Figure 3c: Sub-population 1, age 80 | `additional_analysis_2_80_pdf` | `figures/an280.pdf` |
| Figure 3d: Sub-population 2, age 15 | `additional_analysis_3_15_pdf` | `figures/an315.pdf` |
| Figure 3e: Sub-population 2, age 55 | `additional_analysis_3_55_pdf` | `figures/an355.pdf` |
| Figure 3f: Sub-population 2, age 80 | `additional_analysis_3_80_pdf` | `figures/an380.pdf` |
| Figure 3g: Sub-population 3, age 15 | `additional_analysis_1_15_pdf` | `figures/an115.pdf` |
| Figure 3h: Sub-population 3, age 55 | `additional_analysis_1_55_pdf` | `figures/an155.pdf` |
| Figure 3i: Sub-population 3, age 80 | `additional_analysis_1_80_pdf` | `figures/an180.pdf` |

---

### 8.3 Section 4.3.3 / Figure 4: mean squared error of prediction

| Manuscript panel | Target | Saved file |
|---|---|---|
| Figure 4a: Sub-population 1 | `lc_msep_g2_pdf` | `figures/lc_msep_55_g2.pdf` |
| Figure 4b: Sub-population 2 | `lc_msep_g3_pdf` | `figures/lc_msep_55_g3.pdf` |
| Figure 4c: Sub-population 3 | `lc_msep_g1_pdf` | `figures/lc_msep_55_g1.pdf` |

Underlying saved data object:

| Output | Target | Saved file |
|---|---|---|
| expected quadratic forecast error study object | `lc_msep_results_file` | `output/lc_msep_results.rds` |

---

### 8.4 Section 4.4.2 / Figure 6: rolling-window out-of-sample performance

The current manuscript version includes only the sub-population panels.

| Manuscript panel | Target | Saved file |
|---|---|---|
| Figure 6a: Sub-population 1, MSE | `lc_mse_group_2_pdf` | `figures/lc_mse_age_classes_group_2.pdf` |
| Figure 6b: Sub-population 1, Deviance | `lc_deviance_group_2_pdf` | `figures/lc_deviance_age_classes_group_2.pdf` |
| Figure 6c: Sub-population 2, MSE | `lc_mse_group_3_pdf` | `figures/lc_mse_age_classes_group_3.pdf` |
| Figure 6d: Sub-population 2, Deviance | `lc_deviance_group_3_pdf` | `figures/lc_deviance_age_classes_group_3.pdf` |
| Figure 6e: Sub-population 3, MSE | `lc_mse_group_1_pdf` | `figures/lc_mse_age_classes_group_1.pdf` |
| Figure 6f: Sub-population 3, Deviance | `lc_deviance_group_1_pdf` | `figures/lc_deviance_age_classes_group_1.pdf` |

Underlying saved results file:

| Output | Target | Saved file |
|---|---|---|
| Rolling-window aggregated results | `sim_age_bracket_results_file` | `output/results_age_breakets_smooth.csv` |

The pipeline also retains two auxiliary super-population panels that are not
part of the current manuscript:

| Output | Target | Saved file |
|---|---|---|
| Super-population deviance | `lc_deviance_superpop_pdf` | `figures/lc_deviance_age_classes.pdf` |
| Common-scale super-population deviance | `lc_deviance_superpop_app` | `figures/lc_deviance_age_classes_app.pdf` |

---

### 8.5 Appendix B.3 / Figure 7

| Manuscript panel | Target | Saved file |
|---|---|---|
| Figure 7a: Mortality-trend panel | `fig6_pdf` | `figures/mortality_rates_simulation.pdf` |

The right-hand panel of Figure 7 is a manuscript-side Lexis diagram and is not generated by this R supplement.

---

### 8.6 Appendix B.4 / Figure 8: Figure 4 on comparable scales

The pipeline generates the common-scale versions of the expected quadratic
forecast error fan plots used in the manuscript appendix.

| Manuscript panel | Target | Saved file |
|---|---|---|
| Figure 8a: Sub-population 1 | `lc_msep_g2_app` | `figures/lc_msep_55_g2_appendix.pdf` |
| Figure 8b: Sub-population 2 | `lc_msep_g3_app` | `figures/lc_msep_55_g3_appendix.pdf` |
| Figure 8c: Sub-population 3 | `lc_msep_g1_app` | `figures/lc_msep_55_g1_appendix.pdf` |

---

### 8.7 Appendix B.4 / Figure 9: Figure 6 on comparable scales

The pipeline generates the common-scale versions of the rolling-window plots
used in the manuscript appendix.

| Manuscript panel | Target | Saved file |
|---|---|---|
| Figure 9a: Sub-population 1, MSE | `lc_mse_group_2_app` | `figures/lc_mse_age_classes_group_2_app.pdf` |
| Figure 9b: Sub-population 1, Deviance | `lc_deviance_group_2_app` | `figures/lc_deviance_age_classes_group_2_app.pdf` |
| Figure 9c: Sub-population 2, MSE | `lc_mse_group_3_app` | `figures/lc_mse_age_classes_group_3_app.pdf` |
| Figure 9d: Sub-population 2, Deviance | `lc_deviance_group_3_app` | `figures/lc_deviance_age_classes_group_3_app.pdf` |
| Figure 9e: Sub-population 3, MSE | `lc_mse_group_1_app` | `figures/lc_mse_age_classes_group_1_app.pdf` |
| Figure 9f: Sub-population 3, Deviance | `lc_deviance_group_1_app` | `figures/lc_deviance_age_classes_group_1_app.pdf` |

## 9. Manuscript elements that are not currently generated by this R pipeline

The current submitted manuscript also contains manuscript-side figures or tables that are not direct outputs of the present R pipeline. These include, for example:

- Figure 1
- Figure 5
- Figure 7b
- Table 1

These items should therefore be understood as manuscript-side assets or manually typeset content unless and until dedicated generation code is added to the supplement.

## 10. Table in the manuscript

The manuscript contains a model-comparison table in Section 4.2.

The identification constraints are textual model definitions and are therefore not generated by the R code. However, the average BIC values used in the table can be reconstructed from the rolling-window simulation output file `output/results_age_breakets_smooth.csv`.

### Reconstructing average BIC values

```r
library(data.table)

dt <- fread("output/results_age_breakets_smooth.csv")

bic_dt <- unique(dt[, .(predictor, seed_ix, years_ix, bic)])

bic_summary <- bic_dt[, .(
  average_bic = mean(bic),
  average_bic_div_1000 = mean(bic) / 1000
), by = predictor]

bic_summary[order(average_bic)]
```

This reproduces the empirical BIC summary from the saved simulation output.

## 11. Useful saved objects for spot checks

Two saved files are particularly useful for checking the reproducibility of the simulation study without rerunning every downstream plotting step:

### 11.1 Rolling-window results

```text
output/results_age_breakets_smooth.csv
```

This file contains:

- predictor
- seed index
- rolling-window index
- MSE columns
- deviance columns
- BIC
- log-likelihood
- number of observations
- number of parameters
- number of failed separate-model fits

### 11.2 expected quadratic forecast error study object

```text
output/lc_msep_results.rds
```

This file contains the saved expected quadratic forecast error study object used for the expected quadratic forecast error fan plots.

To inspect it:

```r
x <- readRDS("output/lc_msep_results.rds")
names(x)
names(x$groups)
```

## 12. Session information


Run the following in a fresh R session after installing/loading everything required by the supplement:

```r
library(credMortality)
library(targets)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forecast)
library(fs)
library(StMoMo)
library(demography)
library(rpart)

sessionInfo()

R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8
[4] LC_NUMERIC=C                    LC_TIME=Danish_Denmark.utf8    

time zone: Europe/Copenhagen
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] credMortality_0.0.0.9000

loaded via a namespace (and not attached):
 [1] generics_0.1.4      lattice_0.22-7      digest_0.6.39       magrittr_2.0.4      evaluate_1.0.5      grid_4.5.2         
 [7] RColorBrewer_1.1-3  fastmap_1.2.0       maps_3.4.3          Matrix_1.7-4        gnm_1.1-5           nnet_7.3-20        
[13] forecast_9.0.0      spam_2.11-3         StMoMo_0.4.1        viridisLite_0.4.3   scales_1.4.0        fanplot_4.0.1      
[19] cli_3.6.5           rlang_1.1.7         otel_0.2.0          rootSolve_1.8.2.4   tools_4.5.2         parallel_4.5.2     
[25] dplyr_1.2.0         colorspace_2.1-2    ggplot2_4.0.2       curl_7.0.0          rpart_4.1.24        vctrs_0.7.1        
[31] R6_2.6.1            zoo_1.8-15          lifecycle_1.0.5     tseries_0.10-59     relimp_1.0-5        MASS_7.3-65        
[37] pkgconfig_2.0.3     urca_1.3-4          pillar_1.11.1       gtable_0.3.6        glue_1.8.0          data.table_1.18.2.1
[43] quantmod_0.4.28     Rcpp_1.1.1          fields_17.1         xfun_0.56           tibble_3.3.1        lmtest_0.9-40      
[49] tidyselect_1.2.1    knitr_1.51          rstudioapi_0.18.0   farver_2.1.2        qvcalc_1.0.4        htmltools_0.5.9    
[55] nlme_3.1-168        rmarkdown_2.30      xts_0.14.1          timeDate_4052.112   dotCall64_1.2       fracdiff_1.5-3     
[61] compiler_4.5.2      S7_0.2.1            quadprog_1.5-8      TTR_0.24.4 

```

## 13. Notes on reproducibility

The submitted supplement intentionally includes the `_targets/` folder. This
folder is the `targets` cache from the supplied successful run of the pipeline
and is included to facilitate reproducibility checks for computationally
heavier parts of the analysis.

Including the cache allows users to:

- inspect intermediate objects produced by the pipeline,
- perform spot checks of specific steps or outputs without rerunning the entire study,
- regenerate selected downstream figures and summary files more quickly.

Users who wish to rely on the shipped intermediate results may run the pipeline directly and allow `targets` to reuse the existing cache where appropriate.

Users who prefer a completely fresh rerun from source can first delete the cache and then rerun the full pipeline:

```r
library(targets)
tar_destroy(destroy = "all")
tar_make()
```

In many cases, however, the output mapping in this README and the shipped files
are enough for checking the manuscript-to-file correspondence and avoid
unnecessary recomputation.

## 14. Minimal reproduction recipe

The shortest path to regenerate all code-generated outputs:

```r
install.packages(c("devtools", "targets"))
devtools::install_local(".", upgrade = "never")

library(targets)

# optional but recommended for a clean rerun
tar_destroy(destroy = "all")

tar_make()
```

After this completes, the reproduced figure files are in `figures/` and the saved result objects are in `output/`.
