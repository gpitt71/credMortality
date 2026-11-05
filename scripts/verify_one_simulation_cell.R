library(credMortality)

# The first seed in sim_seed_values (set.seed(42); sample.int(100000)[1:30]).
seed_ix <- 61413L
years_ix <- 1L

reference <- data.table::fread("output/results_age_breakets_smooth.csv")
keep <- reference[["seed_ix"]] == seed_ix &
  reference[["years_ix"]] == years_ix
reference <- reference[keep, ]

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

data.table::setDT(reference)
data.table::setDT(fresh)
data.table::setcolorder(fresh, names(reference))
ordering <- c("predictor", "model", "age_breakets")
data.table::setorderv(reference, ordering)
data.table::setorderv(fresh, ordering)

comparison <- all.equal(
  as.data.frame(reference),
  as.data.frame(fresh),
  tolerance = 1e-8,
  check.attributes = FALSE
)
equal_values <- isTRUE(comparison)
numeric_columns <- names(reference)[vapply(reference, is.numeric, logical(1))]
column_differences <- vapply(
  numeric_columns,
  function(column) {
    max(abs(reference[[column]] - fresh[[column]]), na.rm = TRUE)
  },
  numeric(1)
)
column_differences[!is.finite(column_differences)] <- 0
max_abs_difference <- max(column_differences)

cat("seed_ix=", seed_ix, " years_ix=", years_ix, "\n", sep = "")
cat("rows_reference=", nrow(reference), "\n", sep = "")
cat("rows_fresh=", nrow(fresh), "\n", sep = "")
cat("equal_within_tolerance=", equal_values, "\n", sep = "")
cat("max_absolute_numeric_difference=", max_abs_difference, "\n", sep = "")

if (!equal_values) {
  print(comparison)
  for (column in names(reference)) {
    x <- reference[[column]]
    y <- fresh[[column]]
    column_equal <- isTRUE(all.equal(x, y))
    if (!column_equal) {
      if (is.numeric(x) && is.numeric(y)) {
        delta <- abs(x - y)
        cat(
          column,
          ": differing=", sum(delta > 0, na.rm = TRUE),
          ", max_abs_diff=", max(delta, na.rm = TRUE),
          "\n",
          sep = ""
        )
      } else {
        cat(column, ": values or attributes differ\n", sep = "")
      }
    }
  }
  quit(status = 1L)
}
