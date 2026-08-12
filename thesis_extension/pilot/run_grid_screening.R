# Weak-identification wide-grid screening. This is a diagnostic pilot only.
DIAGNOSTIC_SEED <- 20260816L
N <- 500L
KAPPA <- 0.10
TAUS <- c(0.10, 0.50, 0.90)
N_REPLICATIONS <- 10L
ALPHA_GRID <- seq(-3, 5, by = 0.2)
CRITICAL_VALUE <- qchisq(0.95, df = 2)

if (getRversion() != "3.4.3") {
  stop("This screening must run under R 3.4.3; found ", R.version.string, ".")
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
plot_dir <- file.path(extension_dir, "figures", "pilot_grid")
source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "src", "wn_profiles.R"))
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!identical(ALPHA_GRID, seq(-3, 5, by = 0.2)) ||
    length(ALPHA_GRID) != 41L) {
  stop("The screening alpha grid is not exactly seq(-3, 5, by=0.2).")
}
if (!identical(KAPPA, 0.10)) {
  stop("Kappa is not exactly 0.10.")
}

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

count_components <- function(accepted) {
  accepted_clean <- !is.na(accepted) & accepted
  if (!length(accepted_clean)) return(0L)
  as.integer(sum(accepted_clean & c(TRUE, !head(accepted_clean, -1L))))
}

profile_to_raw <- function(replication, tau, estimator, profile) {
  accepted <- profile$W <= CRITICAL_VALUE
  data.frame(
    replication = replication,
    n = N,
    kappa = KAPPA,
    tau = tau,
    estimator = estimator,
    alpha = profile$grid,
    W = profile$W,
    critical_value = CRITICAL_VALUE,
    accepted = accepted,
    status = profile$status_by_alpha$status,
    stringsAsFactors = FALSE
  )
}

profile_to_summary <- function(replication, tau, estimator, profile) {
  accepted <- profile$W <= CRITICAL_VALUE
  data.frame(
    replication = replication,
    tau = tau,
    estimator = estimator,
    alpha_hat_grid = profile$alpha_hat,
    min_W = profile$min_W,
    n_alpha_success = sum(!is.na(profile$W)),
    n_alpha_failure = sum(is.na(profile$W)),
    n_accepted_gridpoints = sum(accepted, na.rm = TRUE),
    left_endpoint_accepted = accepted[1],
    right_endpoint_accepted = accepted[length(accepted)],
    any_accepted = any(accepted, na.rm = TRUE),
    all_accepted = all(!is.na(accepted)) && all(accepted),
    n_accepted_components_grid = count_components(accepted),
    stringsAsFactors = FALSE
  )
}

set.seed(DIAGNOSTIC_SEED)
raw_rows <- vector("list", N_REPLICATIONS * length(TAUS) * 3L)
summary_rows <- vector("list", length(raw_rows))
row_index <- 0L
dataset_count <- 0L
lambda_generation_count <- 0L
lambda_reuse_checks <- logical(0)
order_checks <- data.frame(
  estimator = character(0), max_abs_W_difference = numeric(0),
  exact_match = logical(0), stringsAsFactors = FALSE
)

for (replication in seq_len(N_REPLICATIONS)) {
  dat <- generate_screening_dataset(N, KAPPA)
  dataset_count <- dataset_count + 1L
  for (tau in TAUS) {
    lambda_bc <- bc_pivotal_lambda(
      dat$X100, R = 1000, tau = tau, c = 2, alpha = 0.1
    )
    lambda_generation_count <- lambda_generation_count + 1L
    profiles <- list(
      "Oracle-GMM" = oracle_wn_profile(
        dat$y, dat$D, dat$X10, dat$Z, tau, ALPHA_GRID
      ),
      "Full-GMM" = full_wn_profile(
        dat$y, dat$D, dat$X100, dat$Z, tau, ALPHA_GRID
      ),
      "DML-IVQR" = dml_wn_profile_bc(
        dat$y, dat$D, dat$X100, dat$Z, tau, ALPHA_GRID, lambda_bc
      )
    )
    lambda_reuse_checks <- c(
      lambda_reuse_checks,
      identical(profiles[["DML-IVQR"]]$lambda_bc, lambda_bc)
    )

    for (estimator in names(profiles)) {
      row_index <- row_index + 1L
      profile <- profiles[[estimator]]
      raw_rows[[row_index]] <- profile_to_raw(
        replication, tau, estimator, profile
      )
      summary_rows[[row_index]] <- profile_to_summary(
        replication, tau, estimator, profile
      )
    }

    # Direct order-invariance check on one dataset/quantile, using the same
    # lambda_bc for both DML evaluations. These diagnostic reverse evaluations
    # are not added to the 3,690-row screening data.
    if (replication == 1L && identical(tau, 0.50)) {
      reversed_profiles <- list(
        "Oracle-GMM" = oracle_wn_profile(
          dat$y, dat$D, dat$X10, dat$Z, tau, rev(ALPHA_GRID)
        ),
        "Full-GMM" = full_wn_profile(
          dat$y, dat$D, dat$X100, dat$Z, tau, rev(ALPHA_GRID)
        ),
        "DML-IVQR" = dml_wn_profile_bc(
          dat$y, dat$D, dat$X100, dat$Z, tau, rev(ALPHA_GRID), lambda_bc
        )
      )
      for (estimator in names(profiles)) {
        forward <- profiles[[estimator]]$W
        reversed <- rev(reversed_profiles[[estimator]]$W)
        differences <- abs(forward - reversed)
        order_checks <- rbind(
          order_checks,
          data.frame(
            estimator = estimator,
            max_abs_W_difference = if (anyNA(differences)) NA_real_ else max(differences),
            exact_match = identical(forward, reversed),
            stringsAsFactors = FALSE
          )
        )
      }
      lambda_reuse_checks <- c(
        lambda_reuse_checks,
        identical(reversed_profiles[["DML-IVQR"]]$lambda_bc, lambda_bc)
      )
    }
  }
  cat("Completed screening replication", replication, "of", N_REPLICATIONS, "\n")
}

raw <- do.call(rbind, raw_rows)
summary <- do.call(rbind, summary_rows)
rownames(raw) <- NULL
rownames(summary) <- NULL

expected_rows <- N_REPLICATIONS * length(TAUS) * 3L * length(ALPHA_GRID)
expected_profiles <- N_REPLICATIONS * length(TAUS) * 3L
validation_pass <- nrow(raw) == expected_rows &&
  nrow(summary) == expected_profiles &&
  dataset_count == N_REPLICATIONS &&
  lambda_generation_count == N_REPLICATIONS * length(TAUS) &&
  all(lambda_reuse_checks) &&
  nrow(order_checks) == 3L &&
  all(order_checks$exact_match) &&
  all(order_checks$max_abs_W_difference <= 1e-12)

raw_path <- file.path(pilot_dir, "grid_screening_raw.csv")
profiles_path <- file.path(pilot_dir, "grid_screening_profiles.csv")
summary_path <- file.path(pilot_dir, "grid_screening_summary.csv")
report_path <- file.path(pilot_dir, "grid_screening_report.txt")
write.csv(raw, raw_path, row.names = FALSE)
write.csv(raw, profiles_path, row.names = FALSE)
write.csv(summary, summary_path, row.names = FALSE)

estimators <- c("Oracle-GMM", "Full-GMM", "DML-IVQR")
plot_files <- character(0)
colors <- grDevices::rainbow(N_REPLICATIONS)
for (tau in TAUS) {
  for (estimator in estimators) {
    selected <- raw$tau == tau & raw$estimator == estimator
    plot_data <- raw[selected, ]
    safe_estimator <- tolower(gsub("-", "_", estimator))
    safe_tau <- gsub("\\.", "p", format(tau, nsmall = 2))
    plot_path <- file.path(
      plot_dir, paste0(safe_estimator, "_tau_", safe_tau, ".pdf")
    )
    plot_files <- c(plot_files, plot_path)
    finite_W <- plot_data$W[is.finite(plot_data$W)]
    y_limits <- range(c(finite_W, CRITICAL_VALUE), finite = TRUE)
    grDevices::pdf(plot_path, width = 8, height = 5.5)
    plot(
      range(ALPHA_GRID), y_limits, type = "n",
      xlab = "alpha", ylab = "W_N(alpha)",
      main = paste(estimator, "tau =", format(tau, nsmall = 2),
                   "kappa = 0.10, n = 500")
    )
    for (replication in seq_len(N_REPLICATIONS)) {
      one <- plot_data[plot_data$replication == replication, ]
      lines(one$alpha, one$W, col = colors[replication], lwd = 1)
      points(one$alpha, one$W, col = colors[replication], pch = 16, cex = 0.3)
    }
    abline(h = CRITICAL_VALUE, lty = 2, lwd = 2)
    legend(
      "topright", legend = paste("rep", seq_len(N_REPLICATIONS)),
      col = colors, lty = 1, cex = 0.65, ncol = 2, bty = "n"
    )
    grDevices::dev.off()
  }
}

range_text <- function(x) {
  finite <- x[is.finite(x)]
  if (!length(finite)) return("NA to NA")
  paste(format(min(finite), digits = 16), "to", format(max(finite), digits = 16))
}

group_rows <- vector("list", length(TAUS) * length(estimators))
group_index <- 0L
for (tau in TAUS) {
  for (estimator in estimators) {
    group_index <- group_index + 1L
    s <- summary[summary$tau == tau & summary$estimator == estimator, ]
    group_rows[[group_index]] <- data.frame(
      tau = tau,
      estimator = estimator,
      n_profiles = nrow(s),
      n_alpha_failures = sum(s$n_alpha_failure),
      n_no_accepted = sum(!s$any_accepted),
      n_all_accepted = sum(s$all_accepted),
      n_left_boundary = sum(s$left_endpoint_accepted %in% TRUE),
      n_right_boundary = sum(s$right_endpoint_accepted %in% TRUE),
      n_multiple_components = sum(s$n_accepted_components_grid > 1L),
      alpha_hat_range = range_text(s$alpha_hat_grid),
      min_W_range = range_text(s$min_W),
      stringsAsFactors = FALSE
    )
  }
}
group_summary <- do.call(rbind, group_rows)

report <- c(
  "Weak-identification wide-grid screening report",
  "",
  paste("R version:", R.version.string),
  paste("Diagnostic seed:", DIAGNOSTIC_SEED),
  paste("n:", N),
  paste("kappa:", format(KAPPA, nsmall = 2)),
  paste("taus:", paste(format(TAUS, nsmall = 2), collapse = ", ")),
  paste("estimators:", paste(estimators, collapse = ", ")),
  paste("replications/datasets:", N_REPLICATIONS),
  paste("alpha grid: seq(-3, 5, by=0.2); points:", length(ALPHA_GRID)),
  paste("chi-square critical value:", format(CRITICAL_VALUE, digits = 16)),
  paste("number of profiles:", nrow(summary)),
  paste("raw alpha-level rows:", nrow(raw)),
  paste("number of alpha-level failures:", sum(summary$n_alpha_failure)),
  "",
  "DML pivotal-penalty validation:",
  paste("lambda_bc generations:", lambda_generation_count,
        "(expected", N_REPLICATIONS * length(TAUS), ")"),
  paste("all supplied-lambda exact checks:", all(lambda_reuse_checks)),
  "One lambda_bc was generated per replication-tau and reused for all 41 alpha values.",
  "The same interface permits reuse across kappa when X and tau are unchanged.",
  "",
  "Direct reversed-grid validation (replication 1, tau=0.50):",
  capture.output(print(order_checks, row.names = FALSE, digits = 16)),
  "",
  "Counts and observed ranges by tau and estimator:",
  capture.output(print(group_summary, row.names = FALSE, right = FALSE)),
  "",
  paste("Profiles with no accepted points:", sum(!summary$any_accepted)),
  paste("Profiles with all points accepted:", sum(summary$all_accepted)),
  paste("Profiles touching left boundary:",
        sum(summary$left_endpoint_accepted %in% TRUE)),
  paste("Profiles touching right boundary:",
        sum(summary$right_endpoint_accepted %in% TRUE)),
  paste("Profiles with multiple accepted components:",
        sum(summary$n_accepted_components_grid > 1L)),
  "",
  "Plot files:",
  plot_files,
  "",
  paste("Validation result:", if (validation_pass) "PASS" else "FAIL"),
  "The acceptance indicator is a preliminary screening diagnostic only;",
  "no confidence-region length, interpolation, or automatic expansion was computed."
)
writeLines(report, report_path)
cat(paste(report, collapse = "\n"), "\n")

if (!validation_pass) {
  stop("Wide-grid screening validation failed; see ", report_path, ".")
}
