# Validate the separate deterministic-CV DML W_N(alpha) profile.
GRID_DIAGNOSTIC_SEED <- 20260815L
BASELINE_DIAGNOSTIC_SEED <- 20260814L
CV_SEED <- 2021L
N <- 500L
TAU <- 0.50
KAPPA <- 1
TOLERANCE <- 1e-12

if (getRversion() != "3.4.3") {
  stop("This validation must run under R 3.4.3; found ", R.version.string, ".")
}
required_versions <- c(
  quantreg = "5.34", hdm = "0.2.0", hqreg = "1.4",
  mvtnorm = "1.0-6", doSNOW = "1.0.16"
)
for (package in names(required_versions)) {
  actual <- packageDescription(package, fields = "Version")
  if (is.na(actual) || actual != required_versions[[package]]) {
    stop(package, ": expected version ", required_versions[[package]],
         "; found ", ifelse(is.na(actual), "<not installed>", actual), ".")
  }
}

cat(R.version.string, "\n")
suppressPackageStartupMessages({
  library(quantreg)
  library(hdm)
  library(hqreg)
  library(mvtnorm)
})

root_dir <- getwd()
extension_dir <- file.path(root_dir, "thesis_extension")
pilot_dir <- file.path(extension_dir, "pilot")
source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "src", "wn_profiles.R"))

generate_baseline_dataset <- function(n) {
  sigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)
  epsilon <- rmvnorm(n = n, mean = c(0, 0), sigma = sigma)
  x <- matrix(rnorm(n * 100), ncol = 100)
  X <- matrix(pnorm(x), ncol = 100)
  z <- matrix(cbind(rnorm(n, 0, 1), rnorm(n, 0, 1)), ncol = 2)
  d_original <- z[, 1] + z[, 2] + epsilon[, 2]
  D_original <- pnorm(d_original)
  Z1 <- z[, 1] + rnorm(n, 0, 1) + X[, 2] + X[, 3] + X[, 4]
  Z2 <- z[, 2] + rnorm(n, 0, 1) + X[, 7] + X[, 8] + X[, 9] + X[, 10]
  Z <- matrix(cbind(Z1, Z2), nrow = n)
  w <- rnorm(n, 0, 1)
  treatment <- make_treatment_kappa(z[, 1], z[, 2], epsilon[, 2], w, KAPPA)
  if (!identical(treatment$d_latent, d_original) ||
      !identical(treatment$D, D_original)) {
    stop("kappa=1 baseline treatment identity failed.")
  }
  b <- matrix(c(rep(5, 7), rep(0, 93)))
  y <- c(1 + treatment$D + X %*% b + epsilon[, 1] * treatment$D)
  list(y = y, D = treatment$D, X = X, Z = Z)
}

grid_A <- round(seq(-1, 3, by = 0.1), 1)
grid_B <- round(seq(-3, 3, by = 0.1), 1)
grid_C <- rev(grid_A)
a0 <- 1
grid_lower_then_a0 <- c(round(seq(-3, 0.9, by = 0.1), 1), a0)
grid_higher_then_a0 <- c(round(seq(3, 1.1, by = -0.1), 1), a0)
alpha_key <- function(alpha) sprintf("%.12f", alpha)

set.seed(GRID_DIAGNOSTIC_SEED)
grid_data <- generate_baseline_dataset(N)
initial_rng_state <- .Random.seed

run_fixed <- function(grid) {
  dml_wn_profile_fixedcv(
    grid_data$y, grid_data$D, grid_data$X, grid_data$Z, TAU, grid
  )
}

diagnostic_messages <- character(0)
collect_messages <- function(label, profile) {
  messages <- profile$status_by_alpha$status[
    profile$status_by_alpha$status != "OK"
  ]
  if (length(messages)) {
    diagnostic_messages <<- c(
      diagnostic_messages, paste0(label, ": ", unique(messages))
    )
  }
}

comparison_row <- function(test, left, right,
                           W_a0_alone = NA_real_,
                           W_a0_after_lower = NA_real_,
                           W_a0_after_higher = NA_real_) {
  differences <- abs(left - right)
  complete <- !is.na(differences)
  data.frame(
    estimator = "DML-IVQR fixed CV",
    test = test,
    n_common_alpha = length(left),
    max_abs_W_difference = if (all(complete)) max(differences) else NA_real_,
    mean_abs_W_difference = if (all(complete)) mean(differences) else NA_real_,
    exact_match = all(complete) && identical(as.numeric(left), as.numeric(right)),
    tolerance_match = all(complete) && all(differences <= TOLERANCE),
    W_a0_alone = W_a0_alone,
    W_a0_after_lower = W_a0_after_lower,
    W_a0_after_higher = W_a0_after_higher,
    stringsAsFactors = FALSE
  )
}

.Random.seed <- initial_rng_state
A_first <- run_fixed(grid_A)
.Random.seed <- initial_rng_state
A_second <- run_fixed(grid_A)
.Random.seed <- initial_rng_state
B <- run_fixed(grid_B)
.Random.seed <- initial_rng_state
C <- run_fixed(grid_C)
.Random.seed <- initial_rng_state
alone <- run_fixed(a0)
.Random.seed <- initial_rng_state
after_lower <- run_fixed(grid_lower_then_a0)
.Random.seed <- initial_rng_state
after_higher <- run_fixed(grid_higher_then_a0)
for (item in list(
  A_first = A_first, A_second = A_second, B = B, C = C,
  alone = alone, after_lower = after_lower, after_higher = after_higher
)) {
  collect_messages("fixed-CV profile", item)
}

B_common <- B$W[match(alpha_key(grid_A), alpha_key(B$grid))]
C_common <- C$W[match(alpha_key(grid_A), alpha_key(C$grid))]
W_positions <- c(
  alone$W[1], after_lower$W[length(after_lower$W)],
  after_higher$W[length(after_higher$W)]
)
invariance <- rbind(
  comparison_row("A_same_grid_same_rng", A_first$W, A_second$W),
  comparison_row("B_extra_points_before", A_first$W, B_common),
  comparison_row("C_reversed_grid", A_first$W, C_common),
  comparison_row(
    "D_single_alpha_positions",
    rep(W_positions[1], 2), W_positions[2:3],
    W_a0_alone = W_positions[1],
    W_a0_after_lower = W_positions[2],
    W_a0_after_higher = W_positions[3]
  )
)
rownames(invariance) <- NULL

# Directly verify the seed-generated partition through equivalent CV results.
cv_y <- grid_data$y - a0 * grid_data$D
cv_one <- cv.hqreg(
  grid_data$X, cv_y, method = "quantile", tau = TAU, FUN = "hqreg",
  nfolds = 5, type.measure = "mae", seed = CV_SEED
)
cv_two <- cv.hqreg(
  grid_data$X, cv_y, method = "quantile", tau = TAU, FUN = "hqreg",
  nfolds = 5, type.measure = "mae", seed = CV_SEED
)
set.seed(CV_SEED)
expected_fold_id <- ceiling(sample(seq_len(N)) / N * 5)
cv_explicit_fold <- cv.hqreg(
  grid_data$X, cv_y, method = "quantile", tau = TAU, FUN = "hqreg",
  nfolds = 5, type.measure = "mae", fold.id = expected_fold_id
)
cv_components_equal <- function(left, right) {
  identical(left$cve, right$cve) &&
    identical(left$cvse, right$cvse) &&
    identical(left$lambda, right$lambda) &&
    identical(left$lambda.min, right$lambda.min) &&
    identical(left$lambda.1se, right$lambda.1se) &&
    identical(left$fit$beta, right$fit$beta)
}
seeded_cv_repeat_equal <- cv_components_equal(cv_one, cv_two)
seeded_cv_matches_explicit_fold <- cv_components_equal(cv_one, cv_explicit_fold)
fold_counts <- tabulate(expected_fold_id, nbins = 5)

# Recreate the three prior validation datasets and compare their recorded
# stochastic alpha_hat values against the separate fixed-CV profile.
prior_path <- file.path(pilot_dir, "wn_profile_validation.csv")
if (!file.exists(prior_path)) {
  stop("Missing prior stochastic validation evidence: ", prior_path)
}
prior <- read.csv(prior_path, stringsAsFactors = FALSE)
prior_dml <- prior[prior$estimator == "DML-IVQR", ]
set.seed(BASELINE_DIAGNOSTIC_SEED)
baseline_datasets <- lapply(seq_len(3L), function(index) generate_baseline_dataset(N))
comparison_rows <- vector("list", 9L)
comparison_index <- 0L
for (dataset_id in seq_len(3L)) {
  dat <- baseline_datasets[[dataset_id]]
  for (tau in c(0.10, 0.50, 0.90)) {
    comparison_index <- comparison_index + 1L
    fixed <- dml_wn_profile_fixedcv(dat$y, dat$D, dat$X, dat$Z, tau, grid_A)
    collect_messages(
      paste0("baseline dataset=", dataset_id, ", tau=", tau), fixed
    )
    selected <- prior_dml$dataset_id == dataset_id &
      abs(prior_dml$tau - tau) < .Machine$double.eps^0.5
    if (sum(selected) != 1L) {
      stop("Expected one prior DML validation row for dataset=", dataset_id,
           ", tau=", tau, ".")
    }
    stochastic_hat <- prior_dml$profile_alpha_hat[selected]
    comparison_rows[[comparison_index]] <- data.frame(
      dataset_id = dataset_id,
      n = N,
      tau = tau,
      alpha_hat_stochastic = stochastic_hat,
      alpha_hat_fixedcv = fixed$alpha_hat,
      difference = fixed$alpha_hat - stochastic_hat,
      n_alpha = length(fixed$W),
      n_success = sum(!is.na(fixed$W)),
      n_failure = sum(is.na(fixed$W)),
      stringsAsFactors = FALSE
    )
  }
}
baseline_comparison <- do.call(rbind, comparison_rows)

invariance_pass <- all(invariance$tolerance_match) &&
  all(is.finite(invariance$max_abs_W_difference))
cv_partition_pass <- seeded_cv_repeat_equal && seeded_cv_matches_explicit_fold
baseline_profiles_complete <- all(baseline_comparison$n_failure == 0L) &&
  all(baseline_comparison$n_success == 41L)
overall_pass <- invariance_pass && cv_partition_pass &&
  baseline_profiles_complete && !length(diagnostic_messages)

invariance_path <- file.path(pilot_dir, "fixedcv_grid_invariance.csv")
comparison_path <- file.path(pilot_dir, "fixedcv_baseline_comparison.csv")
report_path <- file.path(pilot_dir, "fixedcv_grid_invariance.txt")
write.csv(invariance, invariance_path, row.names = FALSE)
write.csv(baseline_comparison, comparison_path, row.names = FALSE)

position_table <- invariance[invariance$test == "D_single_alpha_positions",
  c("W_a0_alone", "W_a0_after_lower", "W_a0_after_higher")]
report <- c(
  "Deterministic-CV DML W_N(alpha) validation",
  "",
  paste("R version:", R.version.string),
  paste("Grid diagnostic seed:", GRID_DIAGNOSTIC_SEED),
  paste("Prior baseline dataset seed:", BASELINE_DIAGNOSTIC_SEED),
  paste("cv.hqreg seed:", CV_SEED),
  paste("n:", N),
  paste("kappa:", format(KAPPA, nsmall = 2)),
  paste("Grid-invariance tau:", format(TAU, nsmall = 2)),
  paste("Tolerance:", format(TOLERANCE, scientific = TRUE)),
  "Grid values were canonicalized to one decimal so common alphas are bit-identical.",
  "",
  paste(
    "The DML-IVQR estimator, 5-fold cross-validation procedure, loss function,"
  ),
  "lambda-selection rule, orthogonal score, and GMM criterion are unchanged.",
  "The CV random partition is fixed with seed = 2021 so that W_N(a) at a given",
  "candidate value does not depend on the surrounding numerical alpha grid or",
  "its evaluation order.",
  "",
  "Grid-invariance results:",
  capture.output(print(invariance[, 1:7], row.names = FALSE, digits = 16)),
  "",
  "W_N(1) by evaluation position:",
  capture.output(print(position_table, row.names = FALSE, digits = 16)),
  "",
  "Direct CV partition/equivalent-result verification at alpha=1:",
  paste("Repeated seed=2021 CV result components identical:",
        if (seeded_cv_repeat_equal) "PASS" else "FAIL"),
  paste("Seeded CV matches explicit reconstructed fold.id result:",
        if (seeded_cv_matches_explicit_fold) "PASS" else "FAIL"),
  paste("Reconstructed fold counts:", paste(fold_counts, collapse = ", ")),
  "Compared components: cve, cvse, lambda, lambda.min, lambda.1se, fit$beta.",
  "",
  "Prior stochastic versus fixed-CV DML alpha_hat:",
  capture.output(print(baseline_comparison, row.names = FALSE, digits = 10)),
  "",
  paste("Warnings/errors captured:", length(diagnostic_messages)),
  if (length(diagnostic_messages)) c("Details:", unique(diagnostic_messages)) else
    "Details: none",
  "",
  paste("Overall result:", if (overall_pass) "PASS" else "FAIL")
)
writeLines(report, report_path)
cat(paste(report, collapse = "\n"), "\n")
if (!overall_pass) {
  stop("Fixed-CV validation failed; do not proceed to grid calibration.")
}
