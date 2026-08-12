# Exact baseline-equivalence validation for W_N(alpha) profile functions.
DIAGNOSTIC_SEED <- 20260814L
N_DATASETS <- 3L
N <- 500L
TAUS <- c(0.10, 0.50, 0.90)
BASELINE_GRID <- seq(-1, 3, length = 41)

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
source(file.path(root_dir, "simulation", "fun_callback.R"))
source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "src", "wn_profiles.R"))

set.seed(DIAGNOSTIC_SEED)
datasets <- vector("list", N_DATASETS)
sigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)
for (dataset_id in seq_len(N_DATASETS)) {
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
  treatment <- make_treatment_kappa(z[, 1], z[, 2], epsilon[, 2], w, 1)
  if (!identical(treatment$d_latent, d_original) ||
      !identical(treatment$D, D_original)) {
    stop("kappa=1 treatment identity failed for dataset ", dataset_id, ".")
  }
  b <- matrix(c(rep(5, 7), rep(0, 93)))
  y <- 1 + treatment$D + X %*% b + epsilon[, 1] * treatment$D
  datasets[[dataset_id]] <- list(
    y = c(y), D = treatment$D, X10 = X[, 1:10], X100 = X, Z = Z
  )
}

captured_messages <- character(0)
capture_original <- function(call, context) {
  call_warnings <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      call(),
      warning = function(warning_condition) {
        call_warnings <<- c(call_warnings, conditionMessage(warning_condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error_condition) error_condition
  )
  if (length(call_warnings)) {
    captured_messages <<- c(
      captured_messages,
      paste0(context, " WARNING: ", unique(call_warnings))
    )
  }
  if (inherits(value, "error")) {
    captured_messages <<- c(
      captured_messages,
      paste0(context, " ERROR: ", conditionMessage(value))
    )
  }
  value
}

validation_rows <- vector("list", N_DATASETS * length(TAUS) * 3L)
row_index <- 0L
for (dataset_id in seq_len(N_DATASETS)) {
  dat <- datasets[[dataset_id]]
  for (tau in TAUS) {
    for (estimator in c("Oracle-GMM", "Full-GMM", "DML-IVQR")) {
      row_index <- row_index + 1L
      context <- paste0("dataset=", dataset_id, ", tau=", tau,
                        ", estimator=", estimator)

      if (estimator == "Oracle-GMM") {
        original <- capture_original(
          function() gmm_quantile(dat$y, dat$D, dat$X10, dat$Z, tau),
          context
        )
        profile <- oracle_wn_profile(
          dat$y, dat$D, dat$X10, dat$Z, tau, BASELINE_GRID
        )
      } else if (estimator == "Full-GMM") {
        original <- capture_original(
          function() gmm_quantile(dat$y, dat$D, dat$X100, dat$Z, tau),
          context
        )
        profile <- full_wn_profile(
          dat$y, dat$D, dat$X100, dat$Z, tau, BASELINE_GRID
        )
      } else {
        saved_seed <- .Random.seed
        original <- capture_original(
          function() hdm_quantile(dat$y, dat$D, dat$X100, dat$Z, tau),
          context
        )
        .Random.seed <- saved_seed
        profile <- dml_wn_profile(
          dat$y, dat$D, dat$X100, dat$Z, tau, BASELINE_GRID
        )
      }

      profile_messages <- profile$status_by_alpha$status[
        profile$status_by_alpha$status != "OK"
      ]
      if (length(profile_messages)) {
        captured_messages <- c(
          captured_messages,
          paste0(context, " PROFILE: ", unique(profile_messages))
        )
      }
      original_failed <- inherits(original, "error")
      original_alpha <- if (original_failed) NA_real_ else as.numeric(original)
      n_success <- sum(!is.na(profile$W))
      n_failure <- sum(is.na(profile$W))
      grid_unchanged <- identical(profile$grid, BASELINE_GRID)
      finite_success <- all(is.finite(profile$W[!is.na(profile$W)]))
      exact_match <- !original_failed && !is.na(profile$alpha_hat) &&
        identical(as.numeric(profile$alpha_hat), original_alpha)
      row_pass <- exact_match && length(profile$W) == 41L &&
        grid_unchanged && finite_success && n_failure == 0L
      status <- if (row_pass) "PASS" else paste(
        c(
          if (!exact_match) "alpha_hat mismatch",
          if (length(profile$W) != 41L) "wrong W length",
          if (!grid_unchanged) "grid changed",
          if (!finite_success) "non-finite successful W",
          if (n_failure != 0L) paste(n_failure, "alpha failures"),
          if (original_failed) "original estimator failed"
        ),
        collapse = "; "
      )
      validation_rows[[row_index]] <- data.frame(
        dataset_id = dataset_id,
        n = N,
        tau = tau,
        estimator = estimator,
        original_alpha_hat = original_alpha,
        profile_alpha_hat = profile$alpha_hat,
        exact_match = exact_match,
        n_alpha = length(profile$W),
        n_success = n_success,
        n_failure = n_failure,
        min_W = profile$min_W,
        status = status,
        stringsAsFactors = FALSE
      )
    }
  }
}
validation <- do.call(rbind, validation_rows)

# Repeat one complete DML profile from exactly the same RNG state.
repro_seed <- .Random.seed
repro_one <- dml_wn_profile(
  datasets[[1]]$y, datasets[[1]]$D, datasets[[1]]$X100,
  datasets[[1]]$Z, 0.50, BASELINE_GRID
)
.Random.seed <- repro_seed
repro_two <- dml_wn_profile(
  datasets[[1]]$y, datasets[[1]]$D, datasets[[1]]$X100,
  datasets[[1]]$Z, 0.50, BASELINE_GRID
)
dml_reproducible <- identical(repro_one$W, repro_two$W) &&
  identical(repro_one$status_by_alpha, repro_two$status_by_alpha) &&
  identical(repro_one$alpha_hat, repro_two$alpha_hat)

overall_pass <- all(validation$status == "PASS") && dml_reproducible
csv_path <- file.path(pilot_dir, "wn_profile_validation.csv")
txt_path <- file.path(pilot_dir, "wn_profile_validation.txt")
write.csv(validation, csv_path, row.names = FALSE)

report <- c(
  "Exact W_N(alpha) profile baseline-equivalence validation",
  "",
  paste("R version:", R.version.string),
  paste("Diagnostic seed:", DIAGNOSTIC_SEED),
  paste("Datasets:", N_DATASETS),
  paste("n per dataset:", N),
  paste("Quantiles:", paste(format(TAUS, nsmall = 2), collapse = ", ")),
  paste("Grid:", paste(format(BASELINE_GRID, nsmall = 1), collapse = ", ")),
  "kappa: 1.00 only",
  "",
  "Density implementation: dnorm(e, mean(e), var(e)).",
  "The authors' use of var(e) as dnorm's third argument was preserved exactly.",
  "No intercept was added to the M/J instrument residualization.",
  "DML uses the authors' cv.hqreg settings and lambda-selection logic unchanged.",
  "Before each original DML call, the complete R RNG state was saved; that",
  "state was restored before the corresponding new DML profile call.",
  "",
  "Equivalence results:",
  capture.output(print(validation, row.names = FALSE, digits = 10)),
  "",
  paste("Alpha-level failures:", sum(validation$n_failure)),
  paste("DML repeated-profile W-vector reproducibility:",
        if (dml_reproducible) "PASS" else "FAIL"),
  paste("Warnings/errors captured:", length(captured_messages)),
  if (length(captured_messages)) c("Details:", unique(captured_messages)) else
    "Details: none",
  "",
  paste("Overall result:", if (overall_pass) "PASS" else "FAIL")
)
writeLines(report, txt_path)
cat(paste(report, collapse = "\n"), "\n")
if (!overall_pass) {
  stop("W_N profile baseline-equivalence validation failed; no grid calibration performed.")
}
