# Build the frozen Italy female mortality dataset shipped with credMortality.
#
# Set HMD_ITA_SNAPSHOT to the local ITA-2 directory before running this file.
# The script reads local files only and never downloads data.

snapshot_root <- Sys.getenv("HMD_ITA_SNAPSHOT", unset = "")
if (!nzchar(snapshot_root)) {
  stop("Set HMD_ITA_SNAPSHOT to the local ITA-2 directory.", call. = FALSE)
}

stats_dir <- file.path(snapshot_root, "STATS")
mx_path <- file.path(stats_dir, "Mx_1x1.txt")
exposure_path <- file.path(stats_dir, "Exposures_1x1.txt")

expected_hashes <- c(
  Mx_1x1 = "07ddeb64683310cc61d9686ce25b8a389a0d21621e25298e307380e6f1b5fbaa",
  Exposures_1x1 = "6c8469db3fcf3e0be9a5cb73c0270857909c1a1e68a5db5816fe3cd8f0e819f9"
)
actual_hashes <- c(
  Mx_1x1 = digest::digest(mx_path, algo = "sha256", file = TRUE),
  Exposures_1x1 = digest::digest(exposure_path, algo = "sha256", file = TRUE)
)
if (!identical(actual_hashes, expected_hashes)) {
  stop("The local ITA-2 input files do not match the frozen snapshot.", call. = FALSE)
}

mx <- HMDHFDplus::readHMD(mx_path, fixup = TRUE)
exposure <- HMDHFDplus::readHMD(exposure_path, fixup = TRUE)
years <- sort(unique(mx[["Year"]]))
ages <- unique(mx[["Age"]])
n_ages <- length(ages)
n_years <- length(years)

female_mx <- matrix(mx[["Female"]], nrow = n_ages, ncol = n_years)
female_exposure <- matrix(
  exposure[["Female"]],
  nrow = n_ages,
  ncol = n_years
)
dimnames(female_mx) <- dimnames(female_exposure) <- list(ages, years)
female_mx[female_mx < 0] <- NA
female_exposure[female_exposure < 0] <- NA

ITAFemaleData <- structure(
  list(
    Dxt = female_exposure * female_mx,
    Ext = female_exposure,
    ages = ages,
    years = years,
    type = "central",
    series = "female",
    label = "Italy"
  ),
  class = "StMoMoData"
)

usethis::use_data(ITAFemaleData, overwrite = TRUE, compress = "xz")
