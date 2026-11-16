library(targets)
library(credMortality)

fresh <- data_generator_hmd_lt(
  seed_input = 1,
  exposure_sup = 94500,
  exposure_sub = c(5000, 500)
)
reference <- tar_read(dt_0)

fresh_digest <- digest::digest(fresh, algo = "sha256")
reference_digest <- digest::digest(reference, algo = "sha256")
same_dimensions <- identical(dim(fresh), dim(reference))
same_columns <- identical(names(fresh), names(reference))
value_comparison <- all.equal(
  as.data.frame(fresh),
  as.data.frame(reference),
  check.attributes = FALSE
)
equal_values <- same_dimensions && same_columns && isTRUE(value_comparison)
equal_hashes <- identical(fresh_digest, reference_digest)

cat("all_equal=", equal_values, "\n", sep = "")
cat("hashes_equal=", equal_hashes, "\n", sep = "")
cat("fresh_sha256=", fresh_digest, "\n", sep = "")
cat("reference_sha256=", reference_digest, "\n", sep = "")

if (!equal_values) {
  print(value_comparison)
  quit(status = 1L)
}
