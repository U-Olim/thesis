# Alpha-grid resolution diagnostic on the prespecified parameter space [-1,3].
DIAGNOSTIC_SEED <- 20260816L
N <- 500L
KAPPA <- 0.10
TAUS <- c(0.10, 0.50, 0.90)
N_REPLICATIONS <- 10L
GRID_005 <- seq(-1, 3, by = 0.05)
GRID_010 <- seq(-1, 3, by = 0.10)
GRID_020 <- seq(-1, 3, by = 0.20)
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
plot_dir <- file.path(extension_dir, "figures", "pilot_resolution")
prior_path <- file.path(pilot_dir, "grid_screening_raw.csv")
source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "src", "wn_profiles.R"))
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(prior_path)) stop("Missing prior screening data: ", prior_path)
if (!identical(GRID_005, seq(-1, 3, by = 0.05)) || length(GRID_005) != 81L) {
  stop("Fine grid is not exactly seq(-1,3,by=0.05).")
}
if (!identical(KAPPA, 0.10)) stop("Kappa is not exactly 0.10.")

alpha_key <- function(alpha) sprintf("%.12f", alpha)
INDEX_010 <- match(alpha_key(GRID_010), alpha_key(GRID_005))
INDEX_020 <- match(alpha_key(GRID_020), alpha_key(GRID_005))
if (anyNA(INDEX_010) || anyNA(INDEX_020)) {
  stop("Coarse grids are not exact nested subsets of the fine grid.")
}
# Preserve the exact binary alpha values produced by the previous pilot's
# seq(-3,5,by=.2) expression at all shared nominal points. Quantile-regression
# solutions can be discontinuous under sub-machine perturbations of alpha.
PRIOR_GRID_EXPRESSION <- seq(-3, 5, by = 0.2)
prior_common <- PRIOR_GRID_EXPRESSION[
  alpha_key(PRIOR_GRID_EXPRESSION) %in% alpha_key(GRID_020)
]
prior_common <- prior_common[match(alpha_key(GRID_020), alpha_key(prior_common))]
GRID_005_EVAL <- GRID_005
GRID_005_EVAL[INDEX_020] <- prior_common

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
      dat$y, dat$D, dat$X10, dat$Z, tau, GRID_005_EVAL
    ),
    "Full-GMM" = full_wn_profile(
      dat$y, dat$D, dat$X100, dat$Z, tau, GRID_005_EVAL
    ),
    "DML-IVQR" = dml_wn_profile_bc(
      dat$y, dat$D, dat$X100, dat$Z, tau, GRID_005_EVAL, lambda_bc
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
  cat("Completed resolution replication", replication, "of", N_REPLICATIONS, "\n")
}
raw <- do.call(rbind, raw_rows)
rownames(raw) <- NULL

# Gate resolution analysis on common h=.20 values from the previous pilot.
prior <- read.csv(prior_path, stringsAsFactors = FALSE)
prior <- prior[prior$alpha >= -1 - TOLERANCE & prior$alpha <= 3 + TOLERANCE, ]
fine_020 <- raw[alpha_key(raw$alpha) %in% alpha_key(GRID_020), ]
comparison_keys <- c("replication", "tau", "estimator", "alpha")
ordering <- function(dat) order(
  dat$replication, dat$tau, dat$estimator, alpha_key(dat$alpha)
)
prior <- prior[ordering(prior), ]
fine_020 <- fine_020[ordering(fine_020), ]
key_match <- identical(prior$replication, fine_020$replication) &&
  identical(prior$estimator, fine_020$estimator) &&
  all(abs(prior$tau - fine_020$tau) <= TOLERANCE) &&
  all(abs(prior$alpha - fine_020$alpha) <= TOLERANCE)
status_match <- identical(prior$status, fine_020$status)
valid_pair <- is.finite(prior$W) & is.finite(fine_020$W)
failure_pair_match <- identical(is.na(prior$W), is.na(fine_020$W))
W_differences <- abs(prior$W[valid_pair] - fine_020$W[valid_pair])
max_reproduction_difference <- if (length(W_differences)) max(W_differences) else 0
reproduction_pass <- nrow(prior) == 1890L && nrow(fine_020) == 1890L &&
  key_match && status_match && failure_pair_match &&
  max_reproduction_difference <= TOLERANCE
cat("Common h=.20 max |W difference|:",
    format(max_reproduction_difference, scientific = TRUE, digits = 16), "\n")
if (!reproduction_pass) {
  stop("Previous common-point reproduction failed; resolution analysis stopped.")
}

count_components <- function(accepted) {
  clean <- !is.na(accepted) & accepted
  if (!length(clean)) return(0L)
  as.integer(sum(clean & c(TRUE, !head(clean, -1L))))
}

component_ranges <- function(accepted) {
  clean <- !is.na(accepted) & accepted
  starts <- which(clean & c(TRUE, !head(clean, -1L)))
  ends <- which(clean & c(!tail(clean, -1L), TRUE))
  list(starts = starts, ends = ends)
}

summarize_resolution <- function(dat, h) {
  accepted <- dat$accepted
  success <- is.finite(dat$W)
  accepted_alpha <- dat$alpha[!is.na(accepted) & accepted]
  missing_reason <- if (all(success)) "none" else {
    if (!success[1] || !success[length(success)]) {
      "missing endpoint evaluation"
    } else {
      "missing internal alpha evaluation"
    }
  }
  accepted_measure <- if (all(success)) {
    indicator <- as.numeric(accepted)
    sum(h * (head(indicator, -1L) + tail(indicator, -1L)) / 2)
  } else {
    NA_real_
  }
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
    accepted_measure_h = accepted_measure,
    accepted_measure_status = missing_reason,
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
    "0.20" = fine[INDEX_020, ],
    "0.10" = fine[INDEX_010, ],
    "0.05" = fine
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
  coarse_keys_010 <- alpha_key(GRID_010)
  coarse_keys_020 <- alpha_key(GRID_020)
  missed_010 <- 0L
  missed_020 <- 0L
  if (length(fine_components$starts)) {
    for (component_id in seq_along(fine_components$starts)) {
      component_indices <- fine_components$starts[component_id]:fine_components$ends[component_id]
      component_alpha_keys <- alpha_key(fine$alpha[component_indices])
      accepted_keys <- component_alpha_keys[fine$accepted[component_indices] %in% TRUE]
      if (!any(accepted_keys %in% coarse_keys_010)) missed_010 <- missed_010 + 1L
      if (!any(accepted_keys %in% coarse_keys_020)) missed_020 <- missed_020 + 1L
    }
  }

  s020 <- summaries[["0.20"]]
  s010 <- summaries[["0.10"]]
  s005 <- summaries[["0.05"]]
  comparison_rows[[profile_index]] <- data.frame(
    replication = key$replication, tau = key$tau, estimator = key$estimator,
    accepted_measure_020 = s020$accepted_measure_h,
    accepted_measure_010 = s010$accepted_measure_h,
    accepted_measure_005 = s005$accepted_measure_h,
    measure_diff_020_005 = s020$accepted_measure_h - s005$accepted_measure_h,
    measure_absdiff_020_005 = abs(s020$accepted_measure_h - s005$accepted_measure_h),
    measure_diff_010_005 = s010$accepted_measure_h - s005$accepted_measure_h,
    measure_absdiff_010_005 = abs(s010$accepted_measure_h - s005$accepted_measure_h),
    component_diff_020_005 = s020$n_accepted_components_grid - s005$n_accepted_components_grid,
    component_diff_010_005 = s010$n_accepted_components_grid - s005$n_accepted_components_grid,
    all_accepted_changed_020 = !identical(s020$all_accepted, s005$all_accepted),
    all_accepted_changed_010 = !identical(s010$all_accepted, s005$all_accepted),
    left_endpoint_changed_020 = !identical(s020$left_boundary_accepted, s005$left_boundary_accepted),
    left_endpoint_changed_010 = !identical(s010$left_boundary_accepted, s005$left_boundary_accepted),
    right_endpoint_changed_020 = !identical(s020$right_boundary_accepted, s005$right_boundary_accepted),
    right_endpoint_changed_010 = !identical(s010$right_boundary_accepted, s005$right_boundary_accepted),
    any_accepted_changed_020 = !identical(s020$any_accepted, s005$any_accepted),
    any_accepted_changed_010 = !identical(s010$any_accepted, s005$any_accepted),
    n_components_005 = s005$n_accepted_components_grid,
    components_missed_by_010 = missed_010,
    components_missed_by_020 = missed_020,
    comparison_status = if (
      all(c(s020$accepted_measure_status, s010$accepted_measure_status,
            s005$accepted_measure_status) == "none")
    ) "complete" else "one or more resolutions contain a missing evaluation",
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
validation_pass <- nrow(raw) == 7290L && nrow(summary) == 270L &&
  nrow(comparison) == 90L && dataset_count == 10L &&
  lambda_generation_count == 30L && all(lambda_reuse_checks) &&
  reproduction_pass && measure_bounds_pass

raw_path <- file.path(pilot_dir, "grid_resolution_raw.csv")
summary_path <- file.path(pilot_dir, "grid_resolution_summary.csv")
comparison_path <- file.path(pilot_dir, "grid_resolution_comparison.csv")
report_path <- file.path(pilot_dir, "grid_resolution_report.txt")
write.csv(raw, raw_path, row.names = FALSE)
write.csv(summary, summary_path, row.names = FALSE)
write.csv(comparison, comparison_path, row.names = FALSE)

# Choose fine-grid examples mechanically: simplest accepted profile, most
# disconnected profile, and profile with the largest accepted measure.
fine_summary <- summary[abs(summary$h - 0.05) <= TOLERANCE, ]
complete_fine <- fine_summary[
  fine_summary$n_failure == 0L & fine_summary$any_accepted,
]
simple_pool <- complete_fine[!complete_fine$all_accepted, ]
if (!nrow(simple_pool)) simple_pool <- complete_fine
simple <- simple_pool[
  order(simple_pool$n_accepted_components_grid,
        abs(simple_pool$accepted_measure_h - 2)),
][1, ]
disconnected <- complete_fine[
  order(-complete_fine$n_accepted_components_grid,
        -complete_fine$accepted_measure_h),
][1, ]
used_keys <- paste(
  c(simple$replication, disconnected$replication),
  c(simple$tau, disconnected$tau),
  c(simple$estimator, disconnected$estimator), sep = "_"
)
candidate_keys <- paste(
  complete_fine$replication, complete_fine$tau, complete_fine$estimator,
  sep = "_"
)
nearly_all_pool <- complete_fine[!candidate_keys %in% used_keys, ]
if (!nrow(nearly_all_pool)) nearly_all_pool <- complete_fine
nearly_all <- nearly_all_pool[order(-nearly_all_pool$accepted_measure_h), ][1, ]
examples <- rbind(simple, disconnected, nearly_all)
examples$example_type <- c("simple", "highly_disconnected", "all_or_nearly_all")

plot_files <- character(0)
for (row in seq_len(nrow(examples))) {
  example <- examples[row, ]
  dat <- raw[
    raw$replication == example$replication & raw$tau == example$tau &
      raw$estimator == example$estimator,
  ]
  safe_estimator <- tolower(gsub("-", "_", example$estimator))
  safe_tau <- gsub("\\.", "p", format(example$tau, nsmall = 2))
  plot_path <- file.path(
    plot_dir,
    paste0(example$example_type, "_rep_", example$replication, "_",
           safe_estimator, "_tau_", safe_tau, ".pdf")
  )
  plot_files <- c(plot_files, plot_path)
  grDevices::pdf(plot_path, width = 8, height = 5.5)
  plot(
    dat$alpha, dat$W, type = "l", lwd = 1,
    xlab = "alpha", ylab = "W_N(alpha)",
    main = paste(example$example_type, "rep", example$replication,
                 example$estimator, "tau", format(example$tau, nsmall = 2))
  )
  points(dat$alpha, dat$W, pch = 16, cex = 0.35, col = "grey45")
  points(dat$alpha[INDEX_010], dat$W[INDEX_010], pch = 1, cex = 0.7, col = "blue")
  points(dat$alpha[INDEX_020], dat$W[INDEX_020], pch = 4, cex = 0.8, col = "red")
  abline(h = CRITICAL_VALUE, lty = 2, lwd = 2)
  legend(
    "topright", c("h=.05", "h=.10", "h=.20", "critical value"),
    col = c("grey45", "blue", "red", "black"),
    pch = c(16, 1, 4, NA), lty = c(NA, NA, NA, 2), bty = "n"
  )
  grDevices::dev.off()
}

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
    s <- summary[summary$tau == tau & summary$estimator == estimator, ]
    cdat <- comparison[comparison$tau == tau & comparison$estimator == estimator, ]
    extract_h <- function(h) s[abs(s$h - h) <= TOLERANCE, ]
    s020 <- extract_h(0.20); s010 <- extract_h(0.10); s005 <- extract_h(0.05)
    group_rows[[group_index]] <- data.frame(
      tau = tau, estimator = estimator,
      mean_measure_020 = safe_stat(s020$accepted_measure_h, mean),
      median_measure_020 = safe_stat(s020$accepted_measure_h, median),
      mean_measure_010 = safe_stat(s010$accepted_measure_h, mean),
      median_measure_010 = safe_stat(s010$accepted_measure_h, median),
      mean_measure_005 = safe_stat(s005$accepted_measure_h, mean),
      median_measure_005 = safe_stat(s005$accepted_measure_h, median),
      mean_absdiff_010_005 = safe_stat(cdat$measure_absdiff_010_005, mean),
      median_absdiff_010_005 = safe_stat(cdat$measure_absdiff_010_005, median),
      max_absdiff_010_005 = safe_stat(cdat$measure_absdiff_010_005, max),
      mean_absdiff_020_005 = safe_stat(cdat$measure_absdiff_020_005, mean),
      median_absdiff_020_005 = safe_stat(cdat$measure_absdiff_020_005, median),
      max_absdiff_020_005 = safe_stat(cdat$measure_absdiff_020_005, max),
      component_count_changed_010 = sum(cdat$component_diff_010_005 != 0),
      component_count_changed_020 = sum(cdat$component_diff_020_005 != 0),
      components_missed_010 = sum(cdat$components_missed_by_010),
      components_missed_020 = sum(cdat$components_missed_by_020),
      all_accepted_020 = sum(s020$all_accepted),
      all_accepted_010 = sum(s010$all_accepted),
      all_accepted_005 = sum(s005$all_accepted),
      left_accepted_020 = sum(s020$left_boundary_accepted %in% TRUE),
      left_accepted_010 = sum(s010$left_boundary_accepted %in% TRUE),
      left_accepted_005 = sum(s005$left_boundary_accepted %in% TRUE),
      right_accepted_020 = sum(s020$right_boundary_accepted %in% TRUE),
      right_accepted_010 = sum(s010$right_boundary_accepted %in% TRUE),
      right_accepted_005 = sum(s005$right_boundary_accepted %in% TRUE),
      numerical_failures_005 = sum(s005$n_failure),
      stringsAsFactors = FALSE
    )
  }
}
group_summary <- do.call(rbind, group_rows)

report <- c(
  "Alpha-grid resolution diagnostic",
  "",
  paste("R version:", R.version.string),
  paste("Diagnostic seed:", DIAGNOSTIC_SEED),
  paste("n:", N),
  paste("kappa:", format(KAPPA, nsmall = 2)),
  paste("replications:", N_REPLICATIONS),
  "parameter space: [-1,3]",
  "fine grid: seq(-1,3,by=.05); 81 points",
  "h=.10 and h=.20 were derived only as nested fine-grid subsets.",
  "Shared h=.20 entries use the exact binary alpha values generated by",
  "the previous pilot's seq(-3,5,by=.2) expression.",
  paste("critical value:", format(CRITICAL_VALUE, digits = 16)),
  "",
  "Previous common h=.20 reproduction:",
  paste("rows compared:", nrow(fine_020)),
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
  "Resolution diagnostics by tau and estimator:",
  capture.output(print(group_summary, row.names = FALSE, right = FALSE, digits = 8)),
  "",
  "Overall totals:",
  paste("profiles:", nrow(comparison)),
  paste("h=.10 component-count changes:",
        sum(comparison$component_diff_010_005 != 0)),
  paste("h=.20 component-count changes:",
        sum(comparison$component_diff_020_005 != 0)),
  paste("fine components missed by h=.10:",
        sum(comparison$components_missed_by_010)),
  paste("fine components missed by h=.20:",
        sum(comparison$components_missed_by_020)),
  paste("h=.10 vs h=.05 maximum accepted-measure absolute difference:",
        safe_stat(comparison$measure_absdiff_010_005, max)),
  paste("h=.20 vs h=.05 maximum accepted-measure absolute difference:",
        safe_stat(comparison$measure_absdiff_020_005, max)),
  "",
  "Representative plotted profiles:",
  capture.output(print(examples[, c(
    "example_type", "replication", "tau", "estimator",
    "n_accepted_components_grid", "accepted_measure_h"
  )], row.names = FALSE, right = FALSE)),
  "",
  "Plot files:",
  plot_files,
  "",
  paste("Overall validation:", if (validation_pass) "PASS" else "FAIL"),
  "accepted_measure_h is a numerical grid-resolution diagnostic only,",
  "not a final confidence-region length. No smoothing or interpolation was used."
)
writeLines(report, report_path)
cat(paste(report, collapse = "\n"), "\n")

if (!validation_pass) {
  stop("Grid-resolution validation failed; see ", report_path, ".")
}
