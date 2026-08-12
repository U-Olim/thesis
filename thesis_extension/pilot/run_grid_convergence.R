# Final alpha-grid convergence diagnostic on A = [-1,3].
DIAGNOSTIC_SEED <- 20260816L
N <- 500L
KAPPA <- 0.10
TAUS <- c(0.10, 0.50, 0.90)
N_REPLICATIONS <- 10L
GRID_0025 <- seq(-1, 3, by = 0.025)
GRID_005 <- seq(-1, 3, by = 0.05)
GRID_010 <- seq(-1, 3, by = 0.10)
CRITICAL_VALUE <- qchisq(0.95, df = 2)
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

suppressPackageStartupMessages({
  library(quantreg)
  library(hdm)
  library(mvtnorm)
})

root_dir <- getwd()
extension_dir <- file.path(root_dir, "thesis_extension")
pilot_dir <- file.path(extension_dir, "pilot")
prior_path <- file.path(pilot_dir, "grid_resolution_raw.csv")
source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "src", "wn_profiles.R"))

if (!file.exists(prior_path)) stop("Missing prior resolution data: ", prior_path)
if (!identical(GRID_0025, seq(-1, 3, by = 0.025)) ||
    length(GRID_0025) != 161L) {
  stop("Fine grid is not exactly seq(-1,3,by=.025).")
}
if (!identical(KAPPA, 0.10)) stop("Kappa is not exactly 0.10.")

alpha_key <- function(alpha) sprintf("%.12f", alpha)

# Reconstruct the exact evaluated h=.05 vector from the previous runner. Its
# nominal h=.20 entries inherited binary values from seq(-3,5,by=.2).
prior_grid_020 <- seq(-3, 5, by = 0.2)
prior_common_020 <- prior_grid_020[
  alpha_key(prior_grid_020) %in% alpha_key(seq(-1, 3, by = 0.2))
]
prior_common_020 <- prior_common_020[
  match(alpha_key(seq(-1, 3, by = 0.2)), alpha_key(prior_common_020))
]
prior_grid_005_eval <- GRID_005
index_020_in_005 <- match(
  alpha_key(seq(-1, 3, by = 0.2)), alpha_key(prior_grid_005_eval)
)
prior_grid_005_eval[index_020_in_005] <- prior_common_020

INDEX_005 <- match(alpha_key(GRID_005), alpha_key(GRID_0025))
INDEX_010 <- match(alpha_key(GRID_010), alpha_key(GRID_0025))
if (anyNA(INDEX_005) || anyNA(INDEX_010)) {
  stop("The h=.05 and h=.10 grids are not nested fine-grid subsets.")
}
GRID_0025_EVAL <- GRID_0025
GRID_0025_EVAL[INDEX_005] <- prior_grid_005_eval

generate_screening_dataset <- function(n, kappa) {
  sigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)
  epsilon <- rmvnorm(n = n, mean = c(0, 0), sigma = sigma)
  x <- matrix(rnorm(n * 100), ncol = 100)
  X <- matrix(pnorm(x), ncol = 100)
  z1 <- rnorm(n, 0, 1)
  z2 <- rnorm(n, 0, 1)
  v1 <- rnorm(n, 0, 1)
  v2 <- rnorm(n, 0, 1)
  w <- rnorm(n, 0, 1)
  treatment <- make_treatment_kappa(z1, z2, epsilon[, 2], w, kappa)
  Z1 <- z1 + v1 + X[, 2] + X[, 3] + X[, 4]
  Z2 <- z2 + v2 + X[, 7] + X[, 8] + X[, 9] + X[, 10]
  Z <- matrix(cbind(Z1, Z2), nrow = n)
  b <- matrix(c(rep(5, 7), rep(0, 93)))
  y <- c(1 + treatment$D + X %*% b + epsilon[, 1] * treatment$D)
  list(y = y, D = treatment$D, X100 = X, X10 = X[, 1:10], Z = Z)
}

run_profiles <- function(dat, tau, lambda_bc) {
  list(
    "Oracle-GMM" = oracle_wn_profile(
      dat$y, dat$D, dat$X10, dat$Z, tau, GRID_0025_EVAL
    ),
    "Full-GMM" = full_wn_profile(
      dat$y, dat$D, dat$X100, dat$Z, tau, GRID_0025_EVAL
    ),
    "DML-IVQR" = dml_wn_profile_bc(
      dat$y, dat$D, dat$X100, dat$Z, tau, GRID_0025_EVAL, lambda_bc
    )
  )
}

profile_to_raw <- function(replication, tau, estimator, profile) {
  data.frame(
    replication = replication, n = N, kappa = KAPPA, tau = tau,
    estimator = estimator, alpha = profile$grid, W = profile$W,
    critical_value = CRITICAL_VALUE,
    accepted = profile$W <= CRITICAL_VALUE,
    status = profile$status_by_alpha$status,
    stringsAsFactors = FALSE
  )
}

set.seed(DIAGNOSTIC_SEED)
raw_rows <- vector("list", N_REPLICATIONS * length(TAUS) * 3L)
raw_index <- 0L
dataset_count <- 0L
lambda_generation_count <- 0L
lambda_reuse_checks <- logical(0)
for (replication in seq_len(N_REPLICATIONS)) {
  dat <- generate_screening_dataset(N, KAPPA)
  dataset_count <- dataset_count + 1L
  for (tau in TAUS) {
    lambda_bc <- bc_pivotal_lambda(
      dat$X100, R = 1000, tau = tau, c = 2, alpha = 0.1
    )
    lambda_generation_count <- lambda_generation_count + 1L
    profiles <- run_profiles(dat, tau, lambda_bc)
    lambda_reuse_checks <- c(
      lambda_reuse_checks,
      identical(profiles[["DML-IVQR"]]$lambda_bc, lambda_bc)
    )
    for (estimator in names(profiles)) {
      raw_index <- raw_index + 1L
      raw_rows[[raw_index]] <- profile_to_raw(
        replication, tau, estimator, profiles[[estimator]]
      )
    }
  }
  cat("Completed convergence replication", replication, "of", N_REPLICATIONS, "\n")
}
raw <- do.call(rbind, raw_rows)
rownames(raw) <- NULL

# Reproduce all prior h=.05 evaluations before calculating convergence results.
prior <- read.csv(prior_path, stringsAsFactors = FALSE)
current_005 <- raw[alpha_key(raw$alpha) %in% alpha_key(GRID_005), ]
ordering <- function(dat) order(dat$replication, dat$tau, dat$estimator, dat$alpha)
prior <- prior[ordering(prior), ]
current_005 <- current_005[ordering(current_005), ]
key_match <- identical(prior$replication, current_005$replication) &&
  identical(prior$estimator, current_005$estimator) &&
  all(abs(prior$tau - current_005$tau) <= TOLERANCE) &&
  all(abs(prior$alpha - current_005$alpha) <= TOLERANCE)
status_match <- identical(prior$status, current_005$status)
failure_match <- identical(is.na(prior$W), is.na(current_005$W))
valid_pair <- is.finite(prior$W) & is.finite(current_005$W)
W_difference <- abs(prior$W[valid_pair] - current_005$W[valid_pair])
max_reproduction_difference <- if (length(W_difference)) max(W_difference) else 0
reproduction_pass <- nrow(prior) == 7290L && nrow(current_005) == 7290L &&
  key_match && status_match && failure_match &&
  max_reproduction_difference <= TOLERANCE
cat("Common h=.05 max |W difference|:",
    format(max_reproduction_difference, scientific = TRUE, digits = 16), "\n")
if (!reproduction_pass) {
  stop("Previous h=.05 reproduction failed; convergence analysis stopped.")
}

count_components <- function(accepted) {
  clean <- !is.na(accepted) & accepted
  if (!length(clean)) return(0L)
  as.integer(sum(clean & c(TRUE, !head(clean, -1L))))
}

component_ranges <- function(accepted) {
  clean <- !is.na(accepted) & accepted
  list(
    starts = which(clean & c(TRUE, !head(clean, -1L))),
    ends = which(clean & c(!tail(clean, -1L), TRUE))
  )
}

summarize_resolution <- function(dat, h) {
  accepted <- dat$accepted
  success <- is.finite(dat$W)
  accepted_alpha <- dat$alpha[!is.na(accepted) & accepted]
  measure_status <- if (all(success)) "none" else if (
    !success[1] || !success[length(success)]
  ) "missing endpoint evaluation" else "missing internal alpha evaluation"
  measure <- if (all(success)) {
    indicator <- as.numeric(accepted)
    sum(h * (head(indicator, -1L) + tail(indicator, -1L)) / 2)
  } else NA_real_
  successful <- which(success)
  minimizer <- if (length(successful)) successful[which.min(dat$W[successful])] else NA_integer_
  data.frame(
    n_gridpoints = nrow(dat), n_success = sum(success), n_failure = sum(!success),
    n_accepted = sum(accepted, na.rm = TRUE),
    any_accepted = any(accepted, na.rm = TRUE),
    all_accepted = all(!is.na(accepted)) && all(accepted),
    left_boundary_accepted = accepted[1],
    right_boundary_accepted = accepted[length(accepted)],
    n_accepted_components_grid = count_components(accepted),
    minimum_accepted_alpha = if (length(accepted_alpha)) min(accepted_alpha) else NA_real_,
    maximum_accepted_alpha = if (length(accepted_alpha)) max(accepted_alpha) else NA_real_,
    alpha_hat_grid = if (is.na(minimizer)) NA_real_ else dat$alpha[minimizer],
    min_W = if (is.na(minimizer)) NA_real_ else dat$W[minimizer],
    accepted_measure_h = measure,
    accepted_measure_status = measure_status,
    stringsAsFactors = FALSE
  )
}

profile_keys <- unique(raw[c("replication", "tau", "estimator")])
profile_keys <- profile_keys[order(
  profile_keys$replication, profile_keys$tau, profile_keys$estimator
), ]
summary_rows <- vector("list", nrow(profile_keys) * 3L)
comparison_rows <- vector("list", nrow(profile_keys))
summary_index <- 0L

for (profile_index in seq_len(nrow(profile_keys))) {
  key <- profile_keys[profile_index, ]
  fine <- raw[
    raw$replication == key$replication & raw$tau == key$tau &
      raw$estimator == key$estimator,
  ]
  subsets <- list(
    "0.10" = fine[INDEX_010, ],
    "0.05" = fine[INDEX_005, ],
    "0.025" = fine
  )
  summaries <- list()
  for (h_name in names(subsets)) {
    summary_index <- summary_index + 1L
    h <- as.numeric(h_name)
    one <- summarize_resolution(subsets[[h_name]], h)
    summaries[[h_name]] <- one
    summary_rows[[summary_index]] <- cbind(key, h = h, one)
  }

  fine_components <- component_ranges(fine$accepted)
  missed_005 <- 0L
  missed_010 <- 0L
  if (length(fine_components$starts)) {
    for (component_id in seq_along(fine_components$starts)) {
      indices <- fine_components$starts[component_id]:fine_components$ends[component_id]
      accepted_keys <- alpha_key(fine$alpha[indices][fine$accepted[indices] %in% TRUE])
      if (!any(accepted_keys %in% alpha_key(GRID_005))) missed_005 <- missed_005 + 1L
      if (!any(accepted_keys %in% alpha_key(GRID_010))) missed_010 <- missed_010 + 1L
    }
  }

  s010 <- summaries[["0.10"]]
  s005 <- summaries[["0.05"]]
  s0025 <- summaries[["0.025"]]
  comparison_rows[[profile_index]] <- data.frame(
    replication = key$replication, tau = key$tau, estimator = key$estimator,
    accepted_measure_005 = s005$accepted_measure_h,
    accepted_measure_0025 = s0025$accepted_measure_h,
    measure_diff_005_0025 = s005$accepted_measure_h - s0025$accepted_measure_h,
    measure_absdiff_005_0025 = abs(s005$accepted_measure_h - s0025$accepted_measure_h),
    accepted_measure_010 = s010$accepted_measure_h,
    measure_diff_010_0025 = s010$accepted_measure_h - s0025$accepted_measure_h,
    measure_absdiff_010_0025 = abs(s010$accepted_measure_h - s0025$accepted_measure_h),
    component_diff_005_0025 = s005$n_accepted_components_grid - s0025$n_accepted_components_grid,
    component_diff_010_0025 = s010$n_accepted_components_grid - s0025$n_accepted_components_grid,
    n_components_0025 = s0025$n_accepted_components_grid,
    components_missed_by_005 = missed_005,
    components_missed_by_010 = missed_010,
    any_accepted_changed_005 = !identical(s005$any_accepted, s0025$any_accepted),
    all_accepted_changed_005 = !identical(s005$all_accepted, s0025$all_accepted),
    left_endpoint_changed_005 = !identical(s005$left_boundary_accepted, s0025$left_boundary_accepted),
    right_endpoint_changed_005 = !identical(s005$right_boundary_accepted, s0025$right_boundary_accepted),
    any_accepted_changed_010 = !identical(s010$any_accepted, s0025$any_accepted),
    all_accepted_changed_010 = !identical(s010$all_accepted, s0025$all_accepted),
    left_endpoint_changed_010 = !identical(s010$left_boundary_accepted, s0025$left_boundary_accepted),
    right_endpoint_changed_010 = !identical(s010$right_boundary_accepted, s0025$right_boundary_accepted),
    comparison_status = if (all(c(
      s010$accepted_measure_status, s005$accepted_measure_status,
      s0025$accepted_measure_status
    ) == "none")) "complete" else "one or more resolutions contain a missing evaluation",
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, summary_rows)
comparison <- do.call(rbind, comparison_rows)
rownames(summary) <- NULL
rownames(comparison) <- NULL
measure_bounds_pass <- all(
  is.na(summary$accepted_measure_h) |
    (summary$accepted_measure_h >= -TOLERANCE &
       summary$accepted_measure_h <= 4 + TOLERANCE)
)
validation_pass <- nrow(raw) == 14490L && nrow(summary) == 270L &&
  nrow(comparison) == 90L && dataset_count == 10L &&
  lambda_generation_count == 30L && all(lambda_reuse_checks) &&
  reproduction_pass && measure_bounds_pass

raw_path <- file.path(pilot_dir, "grid_convergence_raw.csv")
summary_path <- file.path(pilot_dir, "grid_convergence_summary.csv")
comparison_path <- file.path(pilot_dir, "grid_convergence_comparison.csv")
report_path <- file.path(pilot_dir, "grid_convergence_report.txt")
write.csv(raw, raw_path, row.names = FALSE)
write.csv(summary, summary_path, row.names = FALSE)
write.csv(comparison, comparison_path, row.names = FALSE)

safe_stat <- function(x, fun) {
  valid <- x[is.finite(x)]
  if (!length(valid)) NA_real_ else fun(valid)
}
estimators <- c("Oracle-GMM", "Full-GMM", "DML-IVQR")
group_rows <- vector("list", length(TAUS) * length(estimators))
group_index <- 0L
for (tau in TAUS) {
  for (estimator in estimators) {
    group_index <- group_index + 1L
    cdat <- comparison[comparison$tau == tau & comparison$estimator == estimator, ]
    group_rows[[group_index]] <- data.frame(
      tau = tau, estimator = estimator,
      mean_absdiff_005_0025 = safe_stat(cdat$measure_absdiff_005_0025, mean),
      median_absdiff_005_0025 = safe_stat(cdat$measure_absdiff_005_0025, median),
      max_absdiff_005_0025 = safe_stat(cdat$measure_absdiff_005_0025, max),
      component_count_changed_005 = sum(cdat$component_diff_005_0025 != 0),
      components_missed_005 = sum(cdat$components_missed_by_005),
      all_A_changed_005 = sum(cdat$all_accepted_changed_005),
      any_accepted_changed_005 = sum(cdat$any_accepted_changed_005),
      left_endpoint_changed_005 = sum(cdat$left_endpoint_changed_005),
      right_endpoint_changed_005 = sum(cdat$right_endpoint_changed_005),
      numerical_failures = sum(cdat$comparison_status != "complete"),
      fraction_absdiff_le_005 = mean(cdat$measure_absdiff_005_0025 <= 0.05, na.rm = TRUE),
      fraction_absdiff_le_010 = mean(cdat$measure_absdiff_005_0025 <= 0.10, na.rm = TRUE),
      fraction_identical_all_A = mean(!cdat$all_accepted_changed_005),
      fraction_identical_components = mean(cdat$component_diff_005_0025 == 0),
      stringsAsFactors = FALSE
    )
  }
}
group_summary <- do.call(rbind, group_rows)

report <- c(
  "Final alpha-grid convergence diagnostic",
  "",
  paste("R version:", R.version.string),
  paste("Diagnostic seed:", DIAGNOSTIC_SEED),
  paste("n:", N),
  paste("kappa:", format(KAPPA, nsmall = 2)),
  paste("replications:", N_REPLICATIONS),
  "parameter space: [-1,3]",
  "fine grid: seq(-1,3,by=.025); 161 points",
  "h=.05 and h=.10 were derived only as nested fine-grid subsets.",
  "All shared h=.05 entries inherit the exact evaluated alpha values from",
  "the previous resolution runner, including its h=.20 canonicalization.",
  paste("critical value:", format(CRITICAL_VALUE, digits = 16)),
  "",
  "Previous h=.05 reproduction:",
  paste("rows compared:", nrow(current_005)),
  paste("max absolute W difference:",
        format(max_reproduction_difference, scientific = TRUE, digits = 16)),
  paste("tolerance:", format(TOLERANCE, scientific = TRUE)),
  paste("result:", if (reproduction_pass) "PASS" else "FAIL"),
  "",
  paste("fine-grid evaluations:", nrow(raw)),
  paste("fine-grid numerical failures:", sum(is.na(raw$W))),
  paste("lambda_bc generations:", lambda_generation_count, "(expected 30)"),
  paste("all lambda reuse checks:", all(lambda_reuse_checks)),
  paste("accepted measures within [0,4]:", measure_bounds_pass),
  "",
  "h=.05 versus h=.025 by tau and estimator:",
  capture.output(print(group_summary, row.names = FALSE, right = FALSE, digits = 8)),
  "",
  "Overall totals and numerical decision aids:",
  paste("profiles:", nrow(comparison)),
  paste("mean absolute accepted-measure difference:",
        safe_stat(comparison$measure_absdiff_005_0025, mean)),
  paste("median absolute accepted-measure difference:",
        safe_stat(comparison$measure_absdiff_005_0025, median)),
  paste("maximum absolute accepted-measure difference:",
        safe_stat(comparison$measure_absdiff_005_0025, max)),
  paste("profiles with changed component count:",
        sum(comparison$component_diff_005_0025 != 0)),
  paste("fine components missed by h=.05:",
        sum(comparison$components_missed_by_005)),
  paste("fine components missed by h=.10:",
        sum(comparison$components_missed_by_010)),
  paste("all-A classification changes:",
        sum(comparison$all_accepted_changed_005)),
  paste("any-accepted classification changes:",
        sum(comparison$any_accepted_changed_005)),
  paste("left-endpoint classification changes:",
        sum(comparison$left_endpoint_changed_005)),
  paste("right-endpoint classification changes:",
        sum(comparison$right_endpoint_changed_005)),
  paste("fraction with absolute difference <= .05:",
        mean(comparison$measure_absdiff_005_0025 <= 0.05, na.rm = TRUE)),
  paste("fraction with absolute difference <= .10:",
        mean(comparison$measure_absdiff_005_0025 <= 0.10, na.rm = TRUE)),
  paste("fraction with identical all-A classification:",
        mean(!comparison$all_accepted_changed_005)),
  paste("fraction with identical component count:",
        mean(comparison$component_diff_005_0025 == 0)),
  "",
  paste("Overall validation:", if (validation_pass) "PASS" else "FAIL"),
  "These thresholds are descriptive decision aids only. No grid was selected",
  "automatically, and no smoothing, interpolation, coverage, or power was computed."
)
writeLines(report, report_path)
cat(paste(report, collapse = "\n"), "\n")

if (!validation_pass) {
  stop("Grid-convergence validation failed; see ", report_path, ".")
}
