# DGP-only diagnostic. This seed is separate from the authors' Monte Carlo seed
# and is not intended to become the final thesis Monte Carlo seed.
DIAGNOSTIC_SEED <- 20260812L
n_check <- 200000L
p <- 100L

if (getRversion() != "3.4.3") {
  stop("This validation must run under R 3.4.3; found ", R.version.string, ".")
}

required_versions <- c(
  quantreg = "5.34",
  hdm = "0.2.0",
  hqreg = "1.4",
  mvtnorm = "1.0-6",
  doSNOW = "1.0.16"
)
for (package in names(required_versions)) {
  actual <- packageDescription(package, fields = "Version")
  if (is.na(actual) || actual != required_versions[[package]]) {
    stop(package, ": expected version ", required_versions[[package]],
         "; found ", ifelse(is.na(actual), "<not installed>", actual), ".")
  }
}

cat(R.version.string, "\n")

extension_dir <- file.path(getwd(), "thesis_extension")
pilot_dir <- file.path(extension_dir, "pilot")
if (!file.exists(file.path(extension_dir, "src", "dgp_kappa.R"))) {
  stop("Run this validation from the thesis repository root.")
}

source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "config", "kappa_candidates.R"))

set.seed(DIAGNOSTIC_SEED)

# Preserve the authors' primitive random-number generation order exactly:
# (u, epsilon), x, z1, z2, v1, v2. Generate the added w only afterward.
sigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)
error_draws <- mvtnorm::rmvnorm(
  n = n_check,
  mean = c(0, 0),
  sigma = sigma
)
u <- error_draws[, 1]
epsilon <- error_draws[, 2]
x <- matrix(rnorm(n_check * p), ncol = p)
z1 <- rnorm(n_check, 0, 1)
z2 <- rnorm(n_check, 0, 1)
v1 <- rnorm(n_check, 0, 1)
v2 <- rnorm(n_check, 0, 1)
w <- rnorm(n_check, 0, 1)

# x, v1, and v2 are generated to preserve the complete baseline draw order.
# They are not needed for this treatment-index-only diagnostic.
rm(x, v1, v2, error_draws)

d_original <- z1 + z2 + epsilon
D_original <- pnorm(d_original)
baseline <- make_treatment_kappa(z1, z2, epsilon, w, 1)
max_abs_d_identity <- max(abs(baseline$d_latent - d_original))
max_abs_D_identity <- max(abs(baseline$D - D_original))
identity_pass <- identical(max_abs_d_identity, 0) &&
  identical(max_abs_D_identity, 0)

if (!identity_pass) {
  stop(
    "FAIL: kappa=1 baseline identity test; max |d_new-d_original| = ",
    format(max_abs_d_identity, scientific = TRUE),
    "; max |D_new-D_original| = ",
    format(max_abs_D_identity, scientific = TRUE), "."
  )
}

results <- do.call(rbind, lapply(KAPPA_CANDIDATES, function(kappa) {
  treatment <- make_treatment_kappa(z1, z2, epsilon, w, kappa)
  data.frame(
    kappa = kappa,
    var_d_latent = var(treatment$d_latent),
    cov_z1_d = cov(z1, treatment$d_latent),
    cov_z2_d = cov(z2, treatment$d_latent),
    cov_u_d = cov(u, treatment$d_latent),
    cor_u_d = cor(u, treatment$d_latent),
    mean_D = mean(treatment$D),
    sd_D = sd(treatment$D),
    stringsAsFactors = FALSE
  )
}))

# Prespecified diagnostic tolerances; these do not alter the DGP.
moment_pass <- abs(results$var_d_latent - 3) <= 0.03 &
  abs(results$cov_z1_d - results$kappa) <= 0.02 &
  abs(results$cov_z2_d - results$kappa) <= 0.02 &
  abs(results$cov_u_d - 0.3) <= 0.02
distribution_pass <- diff(range(results$mean_D)) <= 0.005 &
  diff(range(results$sd_D)) <= 0.005
overall_pass <- identity_pass && all(moment_pass) && distribution_pass

csv_path <- file.path(pilot_dir, "kappa_dgp_validation.csv")
txt_path <- file.path(pilot_dir, "kappa_dgp_validation.txt")
write.csv(results, csv_path, row.names = FALSE)

report <- c(
  "DML-IVQR kappa DGP validation",
  "",
  paste("R version used:", R.version.string),
  paste("Diagnostic seed:", DIAGNOSTIC_SEED),
  paste("Diagnostic sample size:", n_check),
  "The diagnostic seed is separate from the authors' seed and is not the final Monte Carlo seed.",
  "",
  "Primitive draw order:",
  "(u, epsilon), x, z1, z2, v1, v2, then added w",
  "The same primitive data and w are used for every candidate kappa.",
  "",
  "Exact kappa=1 baseline identity test:",
  paste("max(abs(d_new - d_original)) =", format(max_abs_d_identity, scientific = TRUE)),
  paste("max(abs(D_new - D_original)) =", format(max_abs_D_identity, scientific = TRUE)),
  paste("Identity result:", if (identity_pass) "PASS" else "FAIL"),
  "",
  "Theoretical targets:",
  "Var(d_latent) = 2*kappa^2 + 2*(1-kappa^2) + 1 = 3",
  "Cov(z1, d_latent) = kappa",
  "Cov(z2, d_latent) = kappa",
  "Cov(u, d_latent) = 0.3",
  "",
  "Empirical results:",
  capture.output(print(results, row.names = FALSE, digits = 8)),
  "",
  "Diagnostic tolerances:",
  "|Var(d_latent)-3| <= 0.03",
  "|Cov(zj,d_latent)-kappa| <= 0.02",
  "|Cov(u,d_latent)-0.3| <= 0.02",
  "range(mean(D)) <= 0.005 and range(sd(D)) <= 0.005",
  "",
  paste("Moment checks:", if (all(moment_pass)) "PASS" else "FAIL"),
  paste("Treatment-distribution stability:", if (distribution_pass) "PASS" else "FAIL"),
  paste("Overall summary:", if (overall_pass) "PASS" else "FAIL")
)
writeLines(report, txt_path)

cat(paste(report, collapse = "\n"), "\n")
if (!overall_pass) {
  stop("One or more DGP validation checks failed; see ", txt_path, ".")
}
