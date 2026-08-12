# Diagnostic only: test whether W_N(alpha) changes with grid composition/order.
DIAGNOSTIC_SEED <- 20260815L
N <- 500L
TAU <- 0.50
KAPPA <- 1
TOLERANCE <- 1e-12

if (getRversion() != "3.4.3") {
  stop("This diagnostic must run under R 3.4.3; found ", R.version.string, ".")
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

grid_A <- round(seq(-1, 3, by = 0.1), 1)
grid_B <- round(seq(-3, 3, by = 0.1), 1)
grid_C <- rev(grid_A)
a0 <- 1
grid_lower_then_a0 <- c(round(seq(-3, 0.9, by = 0.1), 1), a0)
grid_higher_then_a0 <- c(round(seq(3, 1.1, by = -0.1), 1), a0)

set.seed(DIAGNOSTIC_SEED)
sigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)
epsilon <- rmvnorm(n = N, mean = c(0, 0), sigma = sigma)
x <- matrix(rnorm(N * 100), ncol = 100)
X <- matrix(pnorm(x), ncol = 100)
z <- matrix(cbind(rnorm(N, 0, 1), rnorm(N, 0, 1)), ncol = 2)
d_original <- z[, 1] + z[, 2] + epsilon[, 2]
D_original <- pnorm(d_original)
Z1 <- z[, 1] + rnorm(N, 0, 1) + X[, 2] + X[, 3] + X[, 4]
Z2 <- z[, 2] + rnorm(N, 0, 1) + X[, 7] + X[, 8] + X[, 9] + X[, 10]
Z <- matrix(cbind(Z1, Z2), nrow = N)
w <- rnorm(N, 0, 1)
treatment <- make_treatment_kappa(z[, 1], z[, 2], epsilon[, 2], w, KAPPA)
if (!identical(treatment$d_latent, d_original) ||
    !identical(treatment$D, D_original)) {
  stop("kappa=1 baseline treatment identity failed.")
}
b <- matrix(c(rep(5, 7), rep(0, 93)))
y <- c(1 + treatment$D + X %*% b + epsilon[, 1] * treatment$D)
X10 <- X[, 1:10]
initial_rng_state <- .Random.seed

run_profile <- function(estimator, grid) {
  if (estimator == "Oracle-GMM") {
    oracle_wn_profile(y, treatment$D, X10, Z, TAU, grid)
  } else if (estimator == "Full-GMM") {
    full_wn_profile(y, treatment$D, X, Z, TAU, grid)
  } else {
    dml_wn_profile(y, treatment$D, X, Z, TAU, grid)
  }
}

alpha_key <- function(alpha) sprintf("%.12f", alpha)
diagnostic_messages <- character(0)
collect_messages <- function(estimator, label, profile) {
  messages <- profile$status_by_alpha$status[
    profile$status_by_alpha$status != "OK"
  ]
  if (length(messages)) {
    diagnostic_messages <<- c(
      diagnostic_messages,
      paste0(estimator, " ", label, ": ", unique(messages))
    )
  }
}

comparison_row <- function(estimator, test, left, right,
                           W_a0_alone = NA_real_,
                           W_a0_after_lower = NA_real_,
                           W_a0_after_higher = NA_real_) {
  differences <- abs(left - right)
  complete <- !is.na(differences)
  max_difference <- if (all(complete)) max(differences) else NA_real_
  mean_difference <- if (all(complete)) mean(differences) else NA_real_
  exact <- all(complete) && identical(as.numeric(left), as.numeric(right))
  tolerance <- all(complete) && all(differences <= TOLERANCE)
  data.frame(
    estimator = estimator,
    test = test,
    n_common_alpha = length(left),
    max_abs_W_difference = max_difference,
    mean_abs_W_difference = mean_difference,
    exact_match = exact,
    tolerance_match = tolerance,
    W_a0_alone = W_a0_alone,
    W_a0_after_lower = W_a0_after_lower,
    W_a0_after_higher = W_a0_after_higher,
    stringsAsFactors = FALSE
  )
}

results <- list()
result_index <- 0L
for (estimator in c("Oracle-GMM", "Full-GMM", "DML-IVQR")) {
  .Random.seed <- initial_rng_state
  A_first <- run_profile(estimator, grid_A)
  .Random.seed <- initial_rng_state
  A_second <- run_profile(estimator, grid_A)
  collect_messages(estimator, "grid_A first", A_first)
  collect_messages(estimator, "grid_A second", A_second)
  result_index <- result_index + 1L
  results[[result_index]] <- comparison_row(
    estimator, "A_same_grid_same_rng", A_first$W, A_second$W
  )

  .Random.seed <- initial_rng_state
  B <- run_profile(estimator, grid_B)
  collect_messages(estimator, "grid_B", B)
  B_common <- B$W[match(alpha_key(grid_A), alpha_key(B$grid))]
  result_index <- result_index + 1L
  results[[result_index]] <- comparison_row(
    estimator, "B_extra_points_before", A_first$W, B_common
  )

  .Random.seed <- initial_rng_state
  C <- run_profile(estimator, grid_C)
  collect_messages(estimator, "grid_C", C)
  C_common <- C$W[match(alpha_key(grid_A), alpha_key(C$grid))]
  result_index <- result_index + 1L
  results[[result_index]] <- comparison_row(
    estimator, "C_reversed_grid", A_first$W, C_common
  )

  .Random.seed <- initial_rng_state
  alone <- run_profile(estimator, a0)
  .Random.seed <- initial_rng_state
  after_lower <- run_profile(estimator, grid_lower_then_a0)
  .Random.seed <- initial_rng_state
  after_higher <- run_profile(estimator, grid_higher_then_a0)
  collect_messages(estimator, "a0 alone", alone)
  collect_messages(estimator, "a0 after lower", after_lower)
  collect_messages(estimator, "a0 after higher", after_higher)
  W_positions <- c(
    alone$W[1],
    after_lower$W[length(after_lower$W)],
    after_higher$W[length(after_higher$W)]
  )
  result_index <- result_index + 1L
  results[[result_index]] <- comparison_row(
    estimator,
    "D_single_alpha_positions",
    rep(W_positions[1], 2),
    W_positions[2:3],
    W_a0_alone = W_positions[1],
    W_a0_after_lower = W_positions[2],
    W_a0_after_higher = W_positions[3]
  )
}
validation <- do.call(rbind, results)
rownames(validation) <- NULL

csv_path <- file.path(pilot_dir, "grid_invariance_validation.csv")
txt_path <- file.path(pilot_dir, "grid_invariance_validation.txt")
write.csv(validation, csv_path, row.names = FALSE)

dml_rows <- validation$estimator == "DML-IVQR"
dml_material_dependence <- any(!validation$tolerance_match[dml_rows])
position_table <- validation[validation$test == "D_single_alpha_positions",
  c("estimator", "W_a0_alone", "W_a0_after_lower", "W_a0_after_higher")]
report <- c(
  "W_N(alpha) grid-invariance diagnostic",
  "",
  paste("R version:", R.version.string),
  paste("Diagnostic seed:", DIAGNOSTIC_SEED),
  paste("n:", N),
  paste("kappa:", format(KAPPA, nsmall = 2)),
  paste("tau:", format(TAU, nsmall = 2)),
  paste("Numerical tolerance:", format(TOLERANCE, scientific = TRUE)),
  "One fixed generated dataset was used for every evaluation.",
  "Grid values were canonicalized to one decimal so matched common alphas are",
  "bit-identical; this avoids comparing differing binary encodings from seq().",
  "Every compared profile was started from the same saved .Random.seed.",
  "No seed argument or fixed folds were added to cv.hqreg.",
  "",
  "Grid-invariance comparisons:",
  capture.output(print(validation[, 1:7], row.names = FALSE, digits = 16)),
  "",
  "W_N(1) by evaluation position:",
  capture.output(print(position_table, row.names = FALSE, digits = 16)),
  "",
  "Authors' example/dml_ivqr_hqreg.Rmd inspection:",
  "The example calls cv.hqreg(..., seed = 2021). It also uses example-specific",
  "arguments including nlambda, lambda, and cv_fold. The authors' replicated",
  "Monte Carlo hdm_quantile implementation does not pass a seed argument.",
  "The current dml_wn_profile mirrors the Monte Carlo implementation and also",
  "does not pass a seed argument.",
  "",
  paste("Warnings/errors captured:", length(diagnostic_messages)),
  if (length(diagnostic_messages)) c("Details:", unique(diagnostic_messages)) else
    "Details: none",
  "",
  paste("DML material grid/order dependence at tolerance:",
        if (dml_material_dependence) "YES" else "NO"),
  if (dml_material_dependence) {
    "STOP: do not proceed to the weak-instrument wide-grid confidence-region pilot."
  } else {
    "DML passed the specified grid-invariance checks."
  }
)
writeLines(report, txt_path)
cat(paste(report, collapse = "\n"), "\n")
