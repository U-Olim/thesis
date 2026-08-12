# Validate the grid-invariant Belloni-Chernozhukov DML W_N(alpha) profile.
DIAGNOSTIC_SEED <- 20260816L
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
  treatment <- make_treatment_kappa(
    z[, 1], z[, 2], epsilon[, 2], w, KAPPA
  )
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
    record_type = "grid_invariance",
    estimator = "DML-IVQR BC",
    test = test,
    dataset_id = NA_integer_,
    n = N,
    tau = TAU,
    kappa = KAPPA,
    n_common_alpha = length(left),
    max_abs_W_difference = if (all(complete)) max(differences) else NA_real_,
    mean_abs_W_difference = if (all(complete)) mean(differences) else NA_real_,
    exact_match = all(complete) && identical(as.numeric(left), as.numeric(right)),
    tolerance_match = all(complete) && all(differences <= TOLERANCE),
    W_a0_alone = W_a0_alone,
    W_a0_after_lower = W_a0_after_lower,
    W_a0_after_higher = W_a0_after_higher,
    alpha_hat = NA_real_,
    min_W = NA_real_,
    n_success = NA_integer_,
    n_failure = NA_integer_,
    lambda_length = NA_integer_,
    lambda_min = NA_real_,
    lambda_median = NA_real_,
    lambda_max = NA_real_,
    stringsAsFactors = FALSE
  )
}

set.seed(DIAGNOSTIC_SEED)
grid_data <- generate_baseline_dataset(N)
lambda_bc <- bc_pivotal_lambda(
  grid_data$X, R = 1000, tau = TAU, c = 2, alpha = 0.1
)
rng_after_lambda <- .Random.seed

run_bc <- function(grid) {
  dml_wn_profile_bc(
    grid_data$y, grid_data$D, grid_data$X, grid_data$Z,
    TAU, grid, lambda_bc
  )
}

A_first <- run_bc(grid_A)
A_second <- run_bc(grid_A)
B <- run_bc(grid_B)
C <- run_bc(grid_C)
alone <- run_bc(a0)
after_lower <- run_bc(grid_lower_then_a0)
after_higher <- run_bc(grid_higher_then_a0)
rng_after_profiles <- .Random.seed

profiles <- list(
  A_first = A_first, A_second = A_second, B = B, C = C,
  alone = alone, after_lower = after_lower, after_higher = after_higher
)
for (label in names(profiles)) {
  collect_messages(label, profiles[[label]])
}

B_common <- B$W[match(alpha_key(grid_A), alpha_key(B$grid))]
C_common <- C$W[match(alpha_key(grid_A), alpha_key(C$grid))]
W_positions <- c(
  alone$W[1], after_lower$W[length(after_lower$W)],
  after_higher$W[length(after_higher$W)]
)
invariance <- rbind(
  comparison_row("A_same_grid", A_first$W, A_second$W),
  comparison_row("B_expanded_grid", A_first$W, B_common),
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

penalty_profile_checks <- vapply(
  profiles, function(profile) identical(profile$lambda_bc, lambda_bc), logical(1)
)
penalty_rng_unchanged <- identical(rng_after_lambda, rng_after_profiles)
penalty_exact_all_profiles <- all(penalty_profile_checks)
penalty_summary <- c(
  length = length(lambda_bc), min = min(lambda_bc),
  median = median(lambda_bc), max = max(lambda_bc)
)

penalty_row <- comparison_row(
  "same_penalty_every_alpha_grid_and_order", 0, 0
)
penalty_row$record_type <- "penalty_invariance"
penalty_row$n_common_alpha <- sum(vapply(profiles, function(x) length(x$grid), integer(1)))
penalty_row$exact_match <- penalty_exact_all_profiles && penalty_rng_unchanged
penalty_row$tolerance_match <- penalty_row$exact_match
penalty_row$lambda_length <- as.integer(penalty_summary["length"])
penalty_row$lambda_min <- penalty_summary["min"]
penalty_row$lambda_median <- penalty_summary["median"]
penalty_row$lambda_max <- penalty_summary["max"]

invariance_pass <- all(invariance$tolerance_match) &&
  all(is.finite(invariance$max_abs_W_difference)) &&
  penalty_row$exact_match
if (!invariance_pass) {
  stop("BC grid/penalty invariance failed; baseline diagnostics were not run.")
}

# Continue the one diagnostic RNG stream for three small baseline datasets.
baseline_rows <- vector("list", 9L)
baseline_index <- 0L
for (dataset_id in seq_len(3L)) {
  dat <- generate_baseline_dataset(N)
  for (tau in c(0.10, 0.50, 0.90)) {
    baseline_index <- baseline_index + 1L
    lambda_tau <- bc_pivotal_lambda(
      dat$X, R = 1000, tau = tau, c = 2, alpha = 0.1
    )
    profile <- dml_wn_profile_bc(
      dat$y, dat$D, dat$X, dat$Z, tau, grid_A, lambda_tau
    )
    collect_messages(
      paste0("baseline dataset=", dataset_id, ", tau=", tau), profile
    )
    baseline_rows[[baseline_index]] <- data.frame(
      record_type = "baseline_alpha_hat",
      estimator = "DML-IVQR BC",
      test = paste0("dataset_", dataset_id, "_tau_", format(tau, nsmall = 2)),
      dataset_id = dataset_id,
      n = N,
      tau = tau,
      kappa = KAPPA,
      n_common_alpha = NA_integer_,
      max_abs_W_difference = NA_real_,
      mean_abs_W_difference = NA_real_,
      exact_match = NA,
      tolerance_match = NA,
      W_a0_alone = NA_real_,
      W_a0_after_lower = NA_real_,
      W_a0_after_higher = NA_real_,
      alpha_hat = profile$alpha_hat,
      min_W = profile$min_W,
      n_success = sum(!is.na(profile$W)),
      n_failure = sum(is.na(profile$W)),
      lambda_length = length(lambda_tau),
      lambda_min = min(lambda_tau),
      lambda_median = median(lambda_tau),
      lambda_max = max(lambda_tau),
      stringsAsFactors = FALSE
    )
  }
}
baseline <- do.call(rbind, baseline_rows)
baseline_complete <- all(baseline$n_success == length(grid_A)) &&
  all(baseline$n_failure == 0L)
overall_pass <- invariance_pass && baseline_complete &&
  !length(diagnostic_messages)

results <- rbind(invariance, penalty_row, baseline)
csv_path <- file.path(pilot_dir, "bc_profile_validation.csv")
report_path <- file.path(pilot_dir, "bc_profile_validation.txt")
write.csv(results, csv_path, row.names = FALSE)

report <- c(
  "Grid-invariant Belloni-Chernozhukov DML W_N(alpha) validation",
  "",
  paste("R version:", R.version.string),
  paste("Diagnostic seed:", DIAGNOSTIC_SEED),
  paste("n:", N),
  paste("kappa:", format(KAPPA, nsmall = 2)),
  paste("Grid-invariance tau:", format(TAU, nsmall = 2)),
  paste("Tolerance:", format(TOLERANCE, scientific = TRUE)),
  "Grid values were canonicalized to one decimal so common alphas are bit-identical.",
  "",
  "The Belloni-Chernozhukov pivotal penalty is computed once for each",
  "dataset and quantile and reused across candidate treatment-effect values.",
  "All numerical constants and the remaining DML-IVQR implementation follow",
  "the authors' Table 3 BC implementation.",
  "",
  "BC penalty implementation:",
  "R = 1000; c = 2; alpha = 0.1; empirical quantile = 0.9",
  "U = matrix(runif(n * R), n)",
  "sigma_j = sqrt(mean(X_j^2)); intercept multiplier = 1",
  "The supplied lambda_bc is not recomputed inside dml_wn_profile_bc().",
  "The interface permits the same lambda_bc to be reused across kappa values",
  "whenever X and tau are unchanged.",
  "",
  "Penalty summary:",
  paste("length(lambda_bc):", penalty_summary["length"]),
  paste("min(lambda_bc):", format(penalty_summary["min"], digits = 16)),
  paste("median(lambda_bc):", format(penalty_summary["median"], digits = 16)),
  paste("max(lambda_bc):", format(penalty_summary["max"], digits = 16)),
  paste("Exact equality in every returned profile:", penalty_exact_all_profiles),
  paste("RNG state unchanged throughout all profile evaluations:", penalty_rng_unchanged),
  paste("Per-profile checks:", paste(names(penalty_profile_checks),
                                      penalty_profile_checks, collapse = "; ")),
  "",
  "Grid-invariance results:",
  capture.output(print(invariance[, c(
    "test", "n_common_alpha", "max_abs_W_difference",
    "mean_abs_W_difference", "exact_match", "tolerance_match"
  )], row.names = FALSE, digits = 16)),
  "",
  "W_N(1) by evaluation position:",
  paste("alone:", format(W_positions[1], digits = 16)),
  paste("after lower alpha values:", format(W_positions[2], digits = 16)),
  paste("after higher alpha values:", format(W_positions[3], digits = 16)),
  "",
  "Baseline BC alpha_hat diagnostics:",
  capture.output(print(baseline[, c(
    "dataset_id", "n", "tau", "kappa", "alpha_hat", "min_W",
    "n_success", "n_failure"
  )], row.names = FALSE, digits = 16)),
  "",
  "Warnings/errors:",
  if (length(diagnostic_messages)) diagnostic_messages else "None.",
  "",
  paste("Overall validation:", if (overall_pass) "PASS" else "FAIL")
)
writeLines(report, report_path)
cat(paste(report, collapse = "\n"), "\n")

if (!overall_pass) {
  stop("BC profile validation failed; see ", report_path, ".")
}
