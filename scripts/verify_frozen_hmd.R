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
equal_values <- isTRUE(all.equal(fresh, reference))

cat("all_equal=", equal_values, "\n", sep = "")
cat("fresh_sha256=", fresh_digest, "\n", sep = "")
cat("reference_sha256=", reference_digest, "\n", sep = "")

if (!equal_values || !identical(fresh_digest, reference_digest)) {
  quit(status = 1L)
}
