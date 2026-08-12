# Descriptive instrument-strength calibration only. This script does not run
# any IVQR estimator, GMM objective, alpha-grid search, or inference procedure.
PILOT_SEED <- 20260813L
N_REPLICATIONS <- 100L
SAMPLE_SIZES <- c(500L, 1000L)
P <- 100L
Q_EXCLUDED <- 2L
K_UNRESTRICTED <- 13L

if (getRversion() != "3.4.3") {
  stop("This pilot must run under R 3.4.3; found ", R.version.string, ".")
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
  stop("Run this pilot from the thesis repository root.")
}
source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "config", "kappa_candidates.R"))

set.seed(PILOT_SEED)
sigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)
warning_log <- character(0)
raw_rows <- vector("list", N_REPLICATIONS * length(SAMPLE_SIZES) *
                    length(KAPPA_CANDIDATES))
row_index <- 0L
identity_checks <- 0L
common_draw_checks <- 0L

empty_row <- function(replication, n, kappa, status) {
  data.frame(
    replication = replication,
    n = n,
    kappa = kappa,
    first_stage_F = NA_real_,
    partial_R2 = NA_real_,
    beta_Z1 = NA_real_,
    beta_Z2 = NA_real_,
    R2_restricted = NA_real_,
    R2_unrestricted = NA_real_,
    RSS_restricted = NA_real_,
    RSS_unrestricted = NA_real_,
    nested_F = NA_real_,
    partial_R2_from_model_R2 = NA_real_,
    cov_z1_d_latent = NA_real_,
    cov_z2_d_latent = NA_real_,
    status = status,
    stringsAsFactors = FALSE
  )
}

for (n in SAMPLE_SIZES) {
  for (replication in seq_len(N_REPLICATIONS)) {
    # Authors' primitive draw order: (u, epsilon), x, z1, z2, v1, v2.
    error_draws <- mvtnorm::rmvnorm(
      n = n,
      mean = c(0, 0),
      sigma = sigma
    )
    u <- error_draws[, 1]
    epsilon <- error_draws[, 2]
    x <- matrix(rnorm(n * P), ncol = P)
    X <- matrix(pnorm(x), ncol = P)
    z <- matrix(cbind(rnorm(n, 0, 1), rnorm(n, 0, 1)), ncol = 2)
    z1 <- z[, 1]
    z2 <- z[, 2]
    Z1 <- z1 + rnorm(n, 0, 1) + X[, 2] + X[, 3] + X[, 4]
    Z2 <- z2 + rnorm(n, 0, 1) + X[, 7] + X[, 8] + X[, 9] + X[, 10]

    # The added primitive is drawn only after every authors' primitive.
    w <- rnorm(n, 0, 1)
    X10 <- X[, 1:10, drop = FALSE]
    colnames(X10) <- paste0("X", seq_len(10L))

    # References make the common-random-number invariants explicit.
    X10_reference <- X10
    Z1_reference <- Z1
    Z2_reference <- Z2
    w_reference <- w

    d_original <- z1 + z2 + epsilon
    D_original <- pnorm(d_original)
    treatment_one <- make_treatment_kappa(z1, z2, epsilon, w, 1)
    if (!identical(treatment_one$d_latent, d_original) ||
        !identical(treatment_one$D, D_original)) {
      stop("kappa=1 identity failed at n=", n,
           ", replication=", replication, ".")
    }
    identity_checks <- identity_checks + 1L

    for (kappa in KAPPA_CANDIDATES) {
      row_index <- row_index + 1L

      if (!identical(X10, X10_reference) ||
          !identical(Z1, Z1_reference) ||
          !identical(Z2, Z2_reference) ||
          !identical(w, w_reference)) {
        stop("Common-draw invariant failed at n=", n,
             ", replication=", replication, ", kappa=", kappa, ".")
      }
      common_draw_checks <- common_draw_checks + 1L

      row_result <- tryCatch(
        withCallingHandlers({
          treatment <- make_treatment_kappa(z1, z2, epsilon, w, kappa)
          restricted_data <- data.frame(D = treatment$D, X10)
          unrestricted_data <- data.frame(
            D = treatment$D,
            X10,
            Z1 = Z1,
            Z2 = Z2
          )

          restricted_fit <- lm(D ~ ., data = restricted_data)
          unrestricted_fit <- lm(D ~ ., data = unrestricted_data)
          if (restricted_fit$rank != 11L || unrestricted_fit$rank != K_UNRESTRICTED) {
            stop("Unexpected regression rank: restricted=", restricted_fit$rank,
                 ", unrestricted=", unrestricted_fit$rank, ".")
          }

          RSS_R <- sum(residuals(restricted_fit)^2)
          RSS_U <- sum(residuals(unrestricted_fit)^2)
          TSS <- sum((treatment$D - mean(treatment$D))^2)
          R2_R <- 1 - RSS_R / TSS
          R2_U <- 1 - RSS_U / TSS

          F_manual <- ((RSS_R - RSS_U) / Q_EXCLUDED) /
            (RSS_U / (n - K_UNRESTRICTED))
          F_nested <- anova(restricted_fit, unrestricted_fit)$F[2]
          partial_R2_rss <- (RSS_R - RSS_U) / RSS_R
          partial_R2_models <- (R2_U - R2_R) / (1 - R2_R)

          f_tolerance <- sqrt(.Machine$double.eps) *
            max(1, abs(F_manual), abs(F_nested))
          r2_tolerance <- sqrt(.Machine$double.eps) *
            max(1, abs(partial_R2_rss), abs(partial_R2_models))
          if (!is.finite(F_manual) || !is.finite(F_nested) ||
              abs(F_manual - F_nested) > f_tolerance) {
            stop("Manual and nested-model F-statistics disagree: manual=",
                 format(F_manual, digits = 16), ", nested=",
                 format(F_nested, digits = 16), ".")
          }
          if (!is.finite(partial_R2_rss) || !is.finite(partial_R2_models) ||
              abs(partial_R2_rss - partial_R2_models) > r2_tolerance) {
            stop("Partial R-squared formulas disagree: RSS=",
                 format(partial_R2_rss, digits = 16), ", model R2=",
                 format(partial_R2_models, digits = 16), ".")
          }

          coefficients_U <- coef(unrestricted_fit)
          data.frame(
            replication = replication,
            n = n,
            kappa = kappa,
            first_stage_F = F_manual,
            partial_R2 = partial_R2_rss,
            beta_Z1 = unname(coefficients_U["Z1"]),
            beta_Z2 = unname(coefficients_U["Z2"]),
            R2_restricted = R2_R,
            R2_unrestricted = R2_U,
            RSS_restricted = RSS_R,
            RSS_unrestricted = RSS_U,
            nested_F = F_nested,
            partial_R2_from_model_R2 = partial_R2_models,
            cov_z1_d_latent = cov(z1, treatment$d_latent),
            cov_z2_d_latent = cov(z2, treatment$d_latent),
            status = "OK",
            stringsAsFactors = FALSE
          )
        }, warning = function(warning_condition) {
          warning_log <<- c(
            warning_log,
            paste0("n=", n, ", replication=", replication,
                   ", kappa=", kappa, ": ", conditionMessage(warning_condition))
          )
          invokeRestart("muffleWarning")
        }),
        error = function(error_condition) {
          empty_row(
            replication,
            n,
            kappa,
            paste0("ERROR: ", conditionMessage(error_condition))
          )
        }
      )
      raw_rows[[row_index]] <- row_result
    }
  }
}

raw_results <- do.call(rbind, raw_rows)
expected_rows <- N_REPLICATIONS * length(SAMPLE_SIZES) *
  length(KAPPA_CANDIDATES)
if (nrow(raw_results) != expected_rows) {
  stop("Expected ", expected_rows, " design rows; found ",
       nrow(raw_results), ".")
}

safe_stats <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(c(mean = NA_real_, median = NA_real_, sd = NA_real_,
             p10 = NA_real_, p90 = NA_real_, min = NA_real_, max = NA_real_))
  }
  c(
    mean = mean(values),
    median = median(values),
    sd = sd(values),
    p10 = unname(quantile(values, 0.10)),
    p90 = unname(quantile(values, 0.90)),
    min = min(values),
    max = max(values)
  )
}

summary_rows <- vector("list", length(SAMPLE_SIZES) * length(KAPPA_CANDIDATES))
summary_index <- 0L
for (n in SAMPLE_SIZES) {
  for (kappa in KAPPA_CANDIDATES) {
    summary_index <- summary_index + 1L
    selected <- raw_results$n == n & raw_results$kappa == kappa
    group <- raw_results[selected, , drop = FALSE]
    successful <- group$status == "OK"
    f_stats <- safe_stats(group$first_stage_F[successful])
    r2_stats <- safe_stats(group$partial_R2[successful])
    summary_rows[[summary_index]] <- data.frame(
      n = n,
      kappa = kappa,
      F_mean = f_stats["mean"],
      F_median = f_stats["median"],
      F_sd = f_stats["sd"],
      F_p10 = f_stats["p10"],
      F_p90 = f_stats["p90"],
      F_min = f_stats["min"],
      F_max = f_stats["max"],
      partial_R2_mean = r2_stats["mean"],
      partial_R2_median = r2_stats["median"],
      partial_R2_sd = r2_stats["sd"],
      partial_R2_p10 = r2_stats["p10"],
      partial_R2_p90 = r2_stats["p90"],
      partial_R2_min = r2_stats["min"],
      partial_R2_max = r2_stats["max"],
      n_success = sum(successful),
      n_failed = sum(!successful),
      stringsAsFactors = FALSE
    )
  }
}
summary_results <- do.call(rbind, summary_rows)
rownames(summary_results) <- NULL

raw_path <- file.path(pilot_dir, "kappa_strength_raw.csv")
summary_path <- file.path(pilot_dir, "kappa_strength_summary.csv")
report_path <- file.path(pilot_dir, "kappa_strength_report.txt")
write.csv(raw_results, raw_path, row.names = FALSE)
write.csv(summary_results, summary_path, row.names = FALSE)

f_report <- summary_results[, c("n", "kappa", "F_mean", "F_median")]
r2_report <- summary_results[, c(
  "n", "kappa", "partial_R2_mean", "partial_R2_median"
)]
n_failed <- sum(raw_results$status != "OK")
report <- c(
  "DML-IVQR descriptive instrument-strength calibration pilot",
  "",
  paste("R version:", R.version.string),
  paste("Pilot seed:", PILOT_SEED),
  paste("Replications per sample size:", N_REPLICATIONS),
  paste("Sample sizes:", paste(SAMPLE_SIZES, collapse = ", ")),
  paste("Candidate kappa values:",
        paste(format(KAPPA_CANDIDATES, nsmall = 2), collapse = ", ")),
  paste("Expected and observed design rows:", expected_rows),
  "",
  "Restricted regression: D(kappa) ~ 1 + X1 + ... + X10",
  "Unrestricted regression: D(kappa) ~ 1 + X1 + ... + X10 + Z1 + Z2",
  "Z1 and Z2 are the observed instruments, not latent z1 and z2.",
  "X1 through X10 are the controls used for this calibration.",
  "Primitive variables are generated once per replication and sample size.",
  "All kappa designs use the same primitives, observed instruments, controls, and w.",
  "",
  paste("kappa=1 exact baseline identity checks passed:", identity_checks),
  paste("Common-draw invariant checks passed:", common_draw_checks),
  "Manual and base-R nested-model F-statistics agreed for every successful row.",
  "Both partial-R2 formulas agreed for every successful row.",
  "",
  "Mean and median joint first-stage F-statistics:",
  capture.output(print(f_report, row.names = FALSE, digits = 8)),
  "",
  "Mean and median partial R-squared:",
  capture.output(print(r2_report, row.names = FALSE, digits = 8)),
  "",
  paste("Failed design rows:", n_failed),
  paste("Warnings captured:", length(warning_log)),
  if (length(warning_log)) c("Warning details:", unique(warning_log)) else
    "Warning details: none",
  "",
  "These statistics are used only to descriptively calibrate the simulation",
  "design. They are not treated as formal IVQR identification statistics."
)
writeLines(report, report_path)
cat(paste(report, collapse = "\n"), "\n")
