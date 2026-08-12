# Three-replication full-design software/integration smoke test.
SMOKE_SEED <- 20260817L
N_REPLICATIONS <- 3L
EXPECTED_RAW_ROWS <- 360L
EXPECTED_POWER_ROWS <- 1440L
CR_H <- 0.05

if (getRversion() != "3.4.3") {
  stop("This smoke test must run under R 3.4.3; found ", R.version.string, ".")
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

suppressPackageStartupMessages({
  library(quantreg)
  library(hdm)
  library(mvtnorm)
})

root_dir <- getwd()
extension_dir <- file.path(root_dir, "thesis_extension")
pilot_dir <- file.path(extension_dir, "pilot")
source(file.path(extension_dir, "config", "kappa_candidates.R"))
source(file.path(extension_dir, "config", "inference_config.R"))
source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "src", "wn_profiles.R"))

if (!identical(SAMPLE_SIZES, c(500, 1000)) ||
    !identical(TAUS, c(0.10, 0.25, 0.50, 0.75, 0.90)) ||
    !identical(KAPPA_CANDIDATES, c(1.00, 0.50, 0.25, 0.10)) ||
    !identical(POINT_GRID, seq(-1, 3, by = 0.10)) ||
    !identical(CR_GRID, seq(-1, 3, by = 0.05)) ||
    !identical(POWER_DELTAS, c(-0.50, -0.25, 0.25, 0.50))) {
  stop("One or more frozen design values do not match the smoke-test specification.")
}

# Hexadecimal keys preserve the exact binary alpha representation. Nominally
# similar but binary-distinct candidates are therefore not silently collapsed.
alpha_key <- function(alpha) sprintf("%a", alpha)
make_union <- function(tau) {
  truth <- alpha_true(tau)
  requested <- c(POINT_GRID, CR_GRID, truth, truth + POWER_DELTAS)
  requested[!duplicated(alpha_key(requested))]
}
map_indices <- function(union_grid, requested) {
  indices <- match(alpha_key(requested), alpha_key(union_grid))
  if (anyNA(indices)) stop("Failed to map an exact requested alpha into the union.")
  indices
}

generate_primitives <- function(n) {
  sigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)
  epsilon <- rmvnorm(n = n, mean = c(0, 0), sigma = sigma)
  x <- matrix(rnorm(n * 100), ncol = 100)
  X <- matrix(pnorm(x), ncol = 100)
  z1 <- rnorm(n, 0, 1)
  z2 <- rnorm(n, 0, 1)
  v1 <- rnorm(n, 0, 1)
  v2 <- rnorm(n, 0, 1)
  w <- rnorm(n, 0, 1)
  Z1 <- z1 + v1 + X[, 2] + X[, 3] + X[, 4]
  Z2 <- z2 + v2 + X[, 7] + X[, 8] + X[, 9] + X[, 10]
  Z <- matrix(cbind(Z1, Z2), nrow = n)
  list(
    epsilon = epsilon, X100 = X, X10 = X[, 1:10],
    z1 = z1, z2 = z2, w = w, Z = Z
  )
}

make_dataset <- function(primitives, kappa) {
  treatment <- make_treatment_kappa(
    primitives$z1, primitives$z2, primitives$epsilon[, 2],
    primitives$w, kappa
  )
  b <- matrix(c(rep(5, 7), rep(0, 93)))
  y <- c(
    1 + treatment$D + primitives$X100 %*% b +
      primitives$epsilon[, 1] * treatment$D
  )
  list(y = y, D = treatment$D)
}

run_profiles <- function(dat, primitives, tau, union_grid, lambda_bc) {
  list(
    "Oracle-GMM" = oracle_wn_profile(
      dat$y, dat$D, primitives$X10, primitives$Z, tau, union_grid
    ),
    "Full-GMM" = full_wn_profile(
      dat$y, dat$D, primitives$X100, primitives$Z, tau, union_grid
    ),
    "DML-IVQR" = dml_wn_profile_bc(
      dat$y, dat$D, primitives$X100, primitives$Z, tau,
      union_grid, lambda_bc
    )
  )
}

count_components <- function(accepted) {
  clean <- !is.na(accepted) & accepted
  if (!length(clean)) return(0L)
  as.integer(sum(clean & c(TRUE, !head(clean, -1L))))
}

status_at <- function(profile, index) profile$status_by_alpha$status[index]
combine_status <- function(messages) {
  messages <- unique(messages[messages != "OK"])
  if (!length(messages)) "OK" else paste(messages, collapse = " | ")
}

raw_rows <- vector("list", EXPECTED_RAW_ROWS)
power_rows <- vector("list", EXPECTED_POWER_ROWS)
raw_index <- 0L
power_index <- 0L
dataset_count <- 0L
lambda_generation_count <- 0L
lambda_reuse_checks <- logical(0)
exact_truth_checks <- logical(0)
exact_power_checks <- logical(0)

set.seed(SMOKE_SEED)
for (replication in seq_len(N_REPLICATIONS)) {
  for (n in SAMPLE_SIZES) {
    primitives <- generate_primitives(n)
    dataset_count <- dataset_count + 1L
    datasets <- lapply(KAPPA_CANDIDATES, function(kappa) {
      make_dataset(primitives, kappa)
    })
    names(datasets) <- format(KAPPA_CANDIDATES, nsmall = 2)

    for (tau in TAUS) {
      truth <- alpha_true(tau)
      false_alphas <- truth + POWER_DELTAS
      union_grid <- make_union(tau)
      point_indices <- map_indices(union_grid, POINT_GRID)
      cr_indices <- map_indices(union_grid, CR_GRID)
      truth_index <- map_indices(union_grid, truth)
      false_indices <- map_indices(union_grid, false_alphas)
      exact_truth_checks <- c(
        exact_truth_checks,
        identical(union_grid[truth_index], truth)
      )
      exact_power_checks <- c(
        exact_power_checks,
        identical(union_grid[false_indices], false_alphas)
      )

      lambda_bc <- bc_pivotal_lambda(
        primitives$X100, R = 1000, tau = tau, c = 2, alpha = 0.1
      )
      lambda_generation_count <- lambda_generation_count + 1L

      for (kappa_index in seq_along(KAPPA_CANDIDATES)) {
        kappa <- KAPPA_CANDIDATES[kappa_index]
        dat <- datasets[[kappa_index]]
        profiles <- run_profiles(dat, primitives, tau, union_grid, lambda_bc)
        lambda_reuse_checks <- c(
          lambda_reuse_checks,
          identical(profiles[["DML-IVQR"]]$lambda_bc, lambda_bc)
        )

        for (estimator in names(profiles)) {
          profile <- profiles[[estimator]]
          point_W <- profile$W[point_indices]
          point_success <- is.finite(point_W)
          successful_point_indices <- which(point_success)
          if (length(successful_point_indices)) {
            point_minimum <- successful_point_indices[
              which.min(point_W[successful_point_indices])
            ]
            alpha_hat <- POINT_GRID[point_minimum]
          } else {
            alpha_hat <- NA_real_
          }

          cr_W <- profile$W[cr_indices]
          cr_success <- is.finite(cr_W)
          cr_accepted <- cr_W <= CRITICAL_VALUE
          if (all(cr_success)) {
            cr_indicator <- as.numeric(cr_accepted)
            cr_length <- sum(
              CR_H * (head(cr_indicator, -1L) + tail(cr_indicator, -1L)) / 2
            )
            cr_length_status <- "OK"
          } else {
            cr_length <- NA_real_
            cr_length_status <- "ERROR: missing required CR-grid W evaluation"
          }

          W_true <- profile$W[truth_index]
          true_status <- status_at(profile, truth_index)
          covered <- if (is.finite(W_true)) W_true <= CRITICAL_VALUE else NA
          rejected_true <- if (is.finite(W_true)) W_true > CRITICAL_VALUE else NA
          point_statuses <- profile$status_by_alpha$status[point_indices]
          cr_statuses <- profile$status_by_alpha$status[cr_indices]
          overall_status <- combine_status(c(
            point_statuses, cr_statuses, true_status, cr_length_status
          ))

          raw_index <- raw_index + 1L
          raw_rows[[raw_index]] <- data.frame(
            replication = replication, n = n, tau = tau, kappa = kappa,
            estimator = estimator, alpha_true = truth, alpha_hat = alpha_hat,
            bias = alpha_hat - truth,
            abs_error = abs(alpha_hat - truth),
            squared_error = (alpha_hat - truth)^2,
            W_true = W_true, covered = covered, rejected_true = rejected_true,
            cr_length = cr_length,
            cr_any_accepted = any(cr_accepted, na.rm = TRUE),
            cr_all_accepted = all(!is.na(cr_accepted)) && all(cr_accepted),
            cr_left_boundary_accepted = cr_accepted[1],
            cr_right_boundary_accepted = cr_accepted[length(cr_accepted)],
            cr_n_components_grid = count_components(cr_accepted),
            point_left_boundary = is.finite(alpha_hat) && alpha_hat == PARAMETER_LOWER,
            point_right_boundary = is.finite(alpha_hat) && alpha_hat == PARAMETER_UPPER,
            n_W_success = sum(is.finite(profile$W)),
            n_W_failure = sum(!is.finite(profile$W)),
            status = overall_status,
            stringsAsFactors = FALSE
          )

          for (delta_index in seq_along(POWER_DELTAS)) {
            false_index <- false_indices[delta_index]
            W_false <- profile$W[false_index]
            false_status <- status_at(profile, false_index)
            power_index <- power_index + 1L
            power_rows[[power_index]] <- data.frame(
              replication = replication, n = n, tau = tau, kappa = kappa,
              estimator = estimator, delta = POWER_DELTAS[delta_index],
              alpha_false = false_alphas[delta_index], W_false = W_false,
              rejected_false = if (is.finite(W_false)) {
                W_false > CRITICAL_VALUE
              } else NA,
              status = false_status,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    cat("Completed smoke replication", replication, "n", n, "\n")
  }
}

raw <- do.call(rbind, raw_rows)
power <- do.call(rbind, power_rows)
rownames(raw) <- NULL
rownames(power) <- NULL

truth_targets <- c(
  `0.10` = -0.281551565544601,
  `0.25` = 0.325510249803919,
  `0.50` = 1,
  `0.75` = 1.67448975019608,
  `0.90` = 2.2815515655446
)
computed_truths <- vapply(TAUS, alpha_true, numeric(1))
truth_validation <- max(abs(computed_truths - truth_targets)) <= 1e-12
power_inside <- all(power$alpha_false >= PARAMETER_LOWER &
                      power$alpha_false <= PARAMETER_UPPER)
cr_length_validation <- all(
  is.na(raw$cr_length) |
    (raw$cr_length >= -1e-12 & raw$cr_length <= 4 + 1e-12)
)
expected_penalties <- N_REPLICATIONS * length(SAMPLE_SIZES) * length(TAUS)
validation_pass <- nrow(raw) == EXPECTED_RAW_ROWS &&
  nrow(power) == EXPECTED_POWER_ROWS && dataset_count == 6L &&
  truth_validation && power_inside && cr_length_validation &&
  lambda_generation_count == expected_penalties && all(lambda_reuse_checks) &&
  all(exact_truth_checks) && all(exact_power_checks)

raw_path <- file.path(pilot_dir, "inference_smoke_raw.csv")
power_path <- file.path(pilot_dir, "inference_smoke_power.csv")
report_path <- file.path(pilot_dir, "inference_smoke_report.txt")
write.csv(raw, raw_path, row.names = FALSE)
write.csv(power, power_path, row.names = FALSE)

range_text <- function(x) {
  finite <- x[is.finite(x)]
  if (!length(finite)) return("NA to NA")
  paste(format(min(finite), digits = 16), "to", format(max(finite), digits = 16))
}
estimator_summary <- do.call(rbind, lapply(
  c("Oracle-GMM", "Full-GMM", "DML-IVQR"), function(estimator) {
    selected <- raw[raw$estimator == estimator, ]
    data.frame(
      estimator = estimator,
      alpha_hat_range = range_text(selected$alpha_hat),
      W_true_range = range_text(selected$W_true),
      CR_length_range = range_text(selected$cr_length),
      stringsAsFactors = FALSE
    )
  }
))

report <- c(
  "Full-design inference smoke-test integration report",
  "",
  paste("R version:", R.version.string),
  paste("Pilot seed:", SMOKE_SEED),
  paste("Replications:", N_REPLICATIONS),
  paste("Sample sizes:", paste(SAMPLE_SIZES, collapse = ", ")),
  paste("Taus:", paste(format(TAUS, nsmall = 2), collapse = ", ")),
  paste("Kappas:", paste(format(KAPPA_CANDIDATES, nsmall = 2), collapse = ", ")),
  "Estimators: Oracle-GMM, Full-GMM, DML-IVQR",
  paste("Point grid: [-1,3], h=.10; points:", length(POINT_GRID)),
  paste("CR grid: [-1,3], h=.05; points:", length(CR_GRID)),
  paste("Critical value:", format(CRITICAL_VALUE, digits = 16)),
  paste("Power deltas:", paste(POWER_DELTAS, collapse = ", ")),
  "",
  paste("Expected raw rows:", EXPECTED_RAW_ROWS),
  paste("Actual raw rows:", nrow(raw)),
  paste("Expected power rows:", EXPECTED_POWER_ROWS),
  paste("Actual power rows:", nrow(power)),
  paste("Generated replication-n primitive datasets:", dataset_count),
  "",
  paste("Union-profile numerical failures:", sum(raw$n_W_failure)),
  paste("Raw rows with non-OK status:", sum(raw$status != "OK")),
  paste("Coverage evaluation failures:", sum(is.na(raw$covered))),
  paste("Power evaluation failures:", sum(is.na(power$rejected_false))),
  paste("CR length failures:", sum(is.na(raw$cr_length))),
  paste("CR lengths outside [0,4]:", sum(
    !is.na(raw$cr_length) & (raw$cr_length < 0 | raw$cr_length > 4)
  )),
  "",
  paste("Point-estimate left-boundary hits:", sum(raw$point_left_boundary)),
  paste("Point-estimate right-boundary hits:", sum(raw$point_right_boundary)),
  paste("CR left-boundary contacts:", sum(raw$cr_left_boundary_accepted %in% TRUE)),
  paste("CR right-boundary contacts:", sum(raw$cr_right_boundary_accepted %in% TRUE)),
  paste("CR all-A accepted:", sum(raw$cr_all_accepted)),
  "",
  paste("BC penalties generated:", lambda_generation_count,
        "(expected", expected_penalties, ")"),
  paste("All BC reuse checks:", all(lambda_reuse_checks)),
  paste("All exact truth-union mapping checks:", all(exact_truth_checks)),
  paste("All exact power-union mapping checks:", all(exact_power_checks)),
  paste("Truth-value validation:", truth_validation),
  paste("All power alternatives inside [-1,3]:", power_inside),
  "Coverage was evaluated directly at exact alpha_true, not on CR_GRID.",
  "Power was evaluated directly at exact alpha_true + Delta, not on CR_GRID.",
  "",
  "Broad software-diagnostic ranges by estimator:",
  capture.output(print(estimator_summary, row.names = FALSE, right = FALSE)),
  "",
  paste("Overall validation:", if (validation_pass) "PASS" else "FAIL"),
  "This is a three-replication integration smoke test. Coverage and power",
  "rates are not interpreted as substantive Monte Carlo findings."
)
writeLines(report, report_path)
cat(paste(report, collapse = "\n"), "\n")

if (!validation_pass) {
  stop("Inference smoke-test validation failed; see ", report_path, ".")
}
