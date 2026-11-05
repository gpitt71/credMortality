# Regenerate only the repeated simulation study and its derived outputs.
#
# From the package root, run:
#   Rscript scripts/run_simulation_study.R
#
# Set SIMULATION_WORKERS to control parallelism. In PowerShell:
#   $env:SIMULATION_WORKERS = "4"
#   Rscript scripts/run_simulation_study.R
#
# The targets cache is retained. If the process is interrupted, launch this
# script again and completed branches will be skipped.

if (!file.exists("DESCRIPTION") || !file.exists("_targets.R")) {
  stop("Run this script from the credMortality package root.", call. = FALSE)
}

workers <- suppressWarnings(as.integer(
  Sys.getenv("SIMULATION_WORKERS", unset = "2")
))
if (length(workers) != 1L || is.na(workers) || workers < 1L) {
  stop("SIMULATION_WORKERS must be a positive integer.", call. = FALSE)
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install the devtools package before running this script.", call. = FALSE)
}
if (!requireNamespace("targets", quietly = TRUE)) {
  stop("Install the targets package before running this script.", call. = FALSE)
}
if (!requireNamespace("future", quietly = TRUE)) {
  stop("Install the future package before running this script.", call. = FALSE)
}
if (!requireNamespace("tidyselect", quietly = TRUE)) {
  stop("Install the tidyselect package before running this script.", call. = FALSE)
}

# The future workers evaluate _targets.R in new R sessions, so install the
# current local package first to ensure every worker uses ITAFemaleData and the
# current analysis functions.
devtools::install_local(
  ".",
  upgrade = "never",
  force = TRUE,
  quiet = TRUE
)

message(
  "Running the repeated simulation study with ",
  workers,
  " worker(s). Existing completed target branches will be reused."
)

targets::tar_make_future(
  names = tidyselect::all_of(c(
    "sim_age_bracket_results_file",
    "sim_bic_table",
    "lc_mse_group_1_pdf",
    "lc_mse_group_1_app",
    "lc_mse_group_2_pdf",
    "lc_mse_group_2_app",
    "lc_mse_group_3_pdf",
    "lc_mse_group_3_app",
    "lc_deviance_superpop_pdf",
    "lc_deviance_superpop_app",
    "lc_deviance_group_1_pdf",
    "lc_deviance_group_1_app",
    "lc_deviance_group_2_pdf",
    "lc_deviance_group_2_app",
    "lc_deviance_group_3_pdf",
    "lc_deviance_group_3_app"
  )),
  workers = workers,
  # Avoid an unnecessary outer callr subprocess. Parallel branches still run
  # in isolated future workers, and this is more reliable when sourced from
  # RStudio on Windows.
  callr_function = NULL
)

# With callr_function = NULL, targets intentionally leaves a process marker
# while the hosting RStudio session remains alive. The make call has returned,
# so the pipeline is finished and the marker can be safely released.
targets::tar_unblock_process()

message("Simulation-study targets completed.")
message("Results: output/results_age_breakets_smooth.csv")
message(
  "Optional spot check: Rscript -e \"devtools::load_all('.'); ",
  "source('scripts/verify_one_simulation_cell.R')\""
)
