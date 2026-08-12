# Wide-tail expansion diagnostic. This remains a pilot, not a final CR rule.
DIAGNOSTIC_SEED <- 20260816L
N <- 500L
KAPPA <- 0.10
TAUS <- c(0.10, 0.50, 0.90)
N_REPLICATIONS <- 10L
ORIGINAL_GRID <- seq(-3, 5, by = 0.2)
TAIL_GRID <- seq(-40, 40, by = 0.2)
WINDOWS <- list(
  "[-3,5]" = c(-3, 5),
  "[-5,5]" = c(-5, 5),
  "[-10,10]" = c(-10, 10),
  "[-20,20]" = c(-20, 20),
  "[-40,40]" = c(-40, 40)
)
CRITICAL_VALUE <- qchisq(0.95, df = 2)
REPRODUCTION_TOLERANCE <- 1e-12

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
plot_dir <- file.path(extension_dir, "figures", "pilot_tail")
prior_path <- file.path(pilot_dir, "grid_screening_raw.csv")
source(file.path(extension_dir, "src", "dgp_kappa.R"))
source(file.path(extension_dir, "src", "wn_profiles.R"))
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(prior_path)) {
  stop("Missing prior screening data: ", prior_path)
}
if (!identical(ORIGINAL_GRID, seq(-3, 5, by = 0.2)) ||
    length(ORIGINAL_GRID) != 41L) {
  stop("Original grid is not exactly seq(-3,5,by=0.2).")
}
if (!identical(TAIL_GRID, seq(-40, 40, by = 0.2)) ||
    length(TAIL_GRID) != 401L) {
  stop("Tail grid is not exactly seq(-40,40,by=0.2).")
}
if (!identical(KAPPA, 0.10)) stop("Kappa is not exactly 0.10.")

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

run_profiles <- function(dat, tau, grid, lambda_bc) {
  list(
    "Oracle-GMM" = oracle_wn_profile(
      dat$y, dat$D, dat$X10, dat$Z, tau, grid
    ),
    "Full-GMM" = full_wn_profile(
      dat$y, dat$D, dat$X100, dat$Z, tau, grid
    ),
    "DML-IVQR" = dml_wn_profile_bc(
      dat$y, dat$D, dat$X100, dat$Z, tau, grid, lambda_bc
    )
  )
}

profile_to_raw <- function(replication, tau, estimator, profile) {
  data.frame(
    replication = replication,
    n = N,
    kappa = KAPPA,
    tau = tau,
    estimator = estimator,
    alpha = profile$grid,
    W = profile$W,
    critical_value = CRITICAL_VALUE,
    accepted = profile$W <= CRITICAL_VALUE,
    status = profile$status_by_alpha$status,
    stringsAsFactors = FALSE
  )
}

alpha_key <- function(alpha) sprintf("%.12f", alpha)
lambda_key <- function(replication, tau) {
  paste(replication, sprintf("%.2f", tau), sep = "_")
}

# Stage 1: recreate the previous pilot in the exact generation order. Keep the
# ten datasets and 30 lambda vectors so Stage 2 does not generate new ones.
set.seed(DIAGNOSTIC_SEED)
stored_datasets <- vector("list", N_REPLICATIONS)
stored_lambdas <- list()
reproduction_rows <- vector("list", N_REPLICATIONS * length(TAUS) * 3L)
reproduction_index <- 0L
lambda_generation_count <- 0L
lambda_reuse_checks <- logical(0)

for (replication in seq_len(N_REPLICATIONS)) {
  dat <- generate_screening_dataset(N, KAPPA)
  stored_datasets[[replication]] <- dat
  for (tau in TAUS) {
    lambda_bc <- bc_pivotal_lambda(
      dat$X100, R = 1000, tau = tau, c = 2, alpha = 0.1
    )
    lambda_generation_count <- lambda_generation_count + 1L
    stored_lambdas[[lambda_key(replication, tau)]] <- lambda_bc
    profiles <- run_profiles(dat, tau, ORIGINAL_GRID, lambda_bc)
    lambda_reuse_checks <- c(
      lambda_reuse_checks,
      identical(profiles[["DML-IVQR"]]$lambda_bc, lambda_bc)
    )
    for (estimator in names(profiles)) {
      reproduction_index <- reproduction_index + 1L
      reproduction_rows[[reproduction_index]] <- profile_to_raw(
        replication, tau, estimator, profiles[[estimator]]
      )
    }
  }
  cat("Completed reproduction replication", replication, "of", N_REPLICATIONS, "\n")
}

reproduced <- do.call(rbind, reproduction_rows)
rownames(reproduced) <- NULL
prior <- read.csv(prior_path, stringsAsFactors = FALSE)
key_columns <- c("replication", "n", "kappa", "tau", "estimator", "alpha")
key_match <- identical(reproduced$replication, prior$replication) &&
  identical(reproduced$n, prior$n) &&
  identical(reproduced$estimator, prior$estimator) &&
  all(abs(reproduced$kappa - prior$kappa) <= REPRODUCTION_TOLERANCE) &&
  all(abs(reproduced$tau - prior$tau) <= REPRODUCTION_TOLERANCE) &&
  all(abs(reproduced$alpha - prior$alpha) <= REPRODUCTION_TOLERANCE)
status_match <- identical(reproduced$status, prior$status)
accepted_match <- identical(reproduced$accepted, prior$accepted)
critical_differences <- abs(reproduced$critical_value - prior$critical_value)
W_differences <- abs(reproduced$W - prior$W)
reproduction_exact_W <- identical(reproduced$W, prior$W)
reproduction_max_abs_W <- if (anyNA(W_differences)) NA_real_ else max(W_differences)
reproduction_pass <- nrow(reproduced) == 3690L && nrow(prior) == 3690L &&
  key_match && status_match && accepted_match &&
  all(critical_differences <= REPRODUCTION_TOLERANCE) &&
  is.finite(reproduction_max_abs_W) &&
  reproduction_max_abs_W <= REPRODUCTION_TOLERANCE &&
  lambda_generation_count == 30L && all(lambda_reuse_checks)

cat("Prior screening reproduction max |W difference|:",
    format(reproduction_max_abs_W, scientific = TRUE, digits = 16), "\n")
cat("Reproduction checks: keys=", key_match,
    "; status=", status_match,
    "; accepted=", accepted_match,
    "; exact_W=", reproduction_exact_W, "\n", sep = "")
if (!reproduction_pass) {
  stop("Prior [-3,5] screening reproduction failed; expanded profiles not run.")
}

# Stage 2: use only the stored datasets and stored lambda vectors.
tail_rows <- vector("list", N_REPLICATIONS * length(TAUS) * 3L)
tail_profiles <- vector("list", length(tail_rows))
tail_index <- 0L
for (replication in seq_len(N_REPLICATIONS)) {
  dat <- stored_datasets[[replication]]
  for (tau in TAUS) {
    lambda_bc <- stored_lambdas[[lambda_key(replication, tau)]]
    profiles <- run_profiles(dat, tau, TAIL_GRID, lambda_bc)
    lambda_reuse_checks <- c(
      lambda_reuse_checks,
      identical(profiles[["DML-IVQR"]]$lambda_bc, lambda_bc)
    )
    for (estimator in names(profiles)) {
      tail_index <- tail_index + 1L
      tail_rows[[tail_index]] <- profile_to_raw(
        replication, tau, estimator, profiles[[estimator]]
      )
      tail_profiles[[tail_index]] <- list(
        replication = replication, tau = tau,
        estimator = estimator, profile = profiles[[estimator]]
      )
    }
  }
  cat("Completed tail replication", replication, "of", N_REPLICATIONS, "\n")
}

tail_raw <- do.call(rbind, tail_rows)
rownames(tail_raw) <- NULL

count_components <- function(accepted) {
  accepted_clean <- !is.na(accepted) & accepted
  if (!length(accepted_clean)) return(0L)
  as.integer(sum(accepted_clean & c(TRUE, !head(accepted_clean, -1L))))
}

value_at <- function(profile, alpha) {
  profile$W[match(alpha_key(alpha), alpha_key(profile$grid))]
}

outer_diagnostics <- function(profile) {
  left <- profile$grid >= -40 - 1e-12 & profile$grid <= -30 + 1e-12
  right <- profile$grid >= 30 - 1e-12 & profile$grid <= 40 + 1e-12
  accepted <- profile$W <= CRITICAL_VALUE
  data.frame(
    any_accept_left_outer = any(accepted[left], na.rm = TRUE),
    any_accept_right_outer = any(accepted[right], na.rm = TRUE),
    all_accept_left_outer = all(!is.na(accepted[left])) && all(accepted[left]),
    all_accept_right_outer = all(!is.na(accepted[right])) && all(accepted[right]),
    W_minus40 = value_at(profile, -40),
    W_minus30 = value_at(profile, -30),
    W_30 = value_at(profile, 30),
    W_40 = value_at(profile, 40),
    stringsAsFactors = FALSE
  )
}

nested_rows <- vector("list", length(tail_profiles) * length(WINDOWS))
nested_index <- 0L
components_rows <- list()
components_index <- 0L
for (item in tail_profiles) {
  profile <- item$profile
  accepted_full <- profile$W <= CRITICAL_VALUE
  outer <- outer_diagnostics(profile)
  for (window_name in names(WINDOWS)) {
    limits <- WINDOWS[[window_name]]
    selected <- profile$grid >= limits[1] - 1e-12 &
      profile$grid <= limits[2] + 1e-12
    alpha <- profile$grid[selected]
    accepted <- accepted_full[selected]
    accepted_alpha <- alpha[!is.na(accepted) & accepted]
    nested_index <- nested_index + 1L
    nested_rows[[nested_index]] <- data.frame(
      replication = item$replication,
      tau = item$tau,
      estimator = item$estimator,
      window = window_name,
      window_left = limits[1],
      window_right = limits[2],
      any_accepted = any(accepted, na.rm = TRUE),
      all_accepted = all(!is.na(accepted)) && all(accepted),
      left_endpoint_accepted = accepted[1],
      right_endpoint_accepted = accepted[length(accepted)],
      n_accepted_gridpoints = sum(accepted, na.rm = TRUE),
      n_accepted_components_grid = count_components(accepted),
      min_accepted_alpha = if (length(accepted_alpha)) min(accepted_alpha) else NA_real_,
      max_accepted_alpha = if (length(accepted_alpha)) max(accepted_alpha) else NA_real_,
      alpha_hat_grid = profile$alpha_hat,
      min_W = profile$min_W,
      n_alpha_success_full = sum(!is.na(profile$W)),
      n_alpha_failure_full = sum(is.na(profile$W)),
      any_accept_left_outer = outer$any_accept_left_outer,
      any_accept_right_outer = outer$any_accept_right_outer,
      all_accept_left_outer = outer$all_accept_left_outer,
      all_accept_right_outer = outer$all_accept_right_outer,
      W_minus40 = outer$W_minus40,
      W_minus30 = outer$W_minus30,
      W_30 = outer$W_30,
      W_40 = outer$W_40,
      left_tail_label = if (outer$any_accept_left_outer) {
        "tail unresolved within search range"
      } else {
        "no acceptance near search boundary"
      },
      right_tail_label = if (outer$any_accept_right_outer) {
        "tail unresolved within search range"
      } else {
        "no acceptance near search boundary"
      },
      stringsAsFactors = FALSE
    )
  }

  accepted_clean <- !is.na(accepted_full) & accepted_full
  starts <- which(accepted_clean & c(TRUE, !head(accepted_clean, -1L)))
  ends <- which(accepted_clean & c(!tail(accepted_clean, -1L), TRUE))
  if (length(starts)) {
    for (component_id in seq_along(starts)) {
      components_index <- components_index + 1L
      components_rows[[components_index]] <- data.frame(
        replication = item$replication,
        tau = item$tau,
        estimator = item$estimator,
        component_id = component_id,
        first_accepted_alpha = profile$grid[starts[component_id]],
        last_accepted_alpha = profile$grid[ends[component_id]],
        n_gridpoints = ends[component_id] - starts[component_id] + 1L,
        stringsAsFactors = FALSE
      )
    }
  }
}

nested_summary <- do.call(rbind, nested_rows)
rownames(nested_summary) <- NULL
if (length(components_rows)) {
  components <- do.call(rbind, components_rows)
  rownames(components) <- NULL
} else {
  components <- data.frame(
    replication = integer(0), tau = numeric(0), estimator = character(0),
    component_id = integer(0), first_accepted_alpha = numeric(0),
    last_accepted_alpha = numeric(0), n_gridpoints = integer(0),
    stringsAsFactors = FALSE
  )
}

# A direct small reversed-grid check, separate from saved tail evaluations.
check_dat <- stored_datasets[[1L]]
check_tau <- 0.50
check_lambda <- stored_lambdas[[lambda_key(1L, check_tau)]]
check_grid <- seq(-2, 2, by = 0.2)
check_forward <- run_profiles(check_dat, check_tau, check_grid, check_lambda)
check_reverse <- run_profiles(check_dat, check_tau, rev(check_grid), check_lambda)
order_checks <- do.call(rbind, lapply(names(check_forward), function(estimator) {
  forward <- check_forward[[estimator]]$W
  reverse <- rev(check_reverse[[estimator]]$W)
  difference <- abs(forward - reverse)
  data.frame(
    estimator = estimator,
    max_abs_W_difference = if (anyNA(difference)) NA_real_ else max(difference),
    exact_match = identical(forward, reverse),
    stringsAsFactors = FALSE
  )
}))
lambda_reuse_checks <- c(
  lambda_reuse_checks,
  identical(check_forward[["DML-IVQR"]]$lambda_bc, check_lambda),
  identical(check_reverse[["DML-IVQR"]]$lambda_bc, check_lambda)
)

expected_tail_rows <- N_REPLICATIONS * length(TAUS) * 3L * length(TAIL_GRID)
validation_pass <- reproduction_pass && nrow(tail_raw) == expected_tail_rows &&
  nrow(nested_summary) == 90L * length(WINDOWS) &&
  lambda_generation_count == 30L && all(lambda_reuse_checks) &&
  all(order_checks$exact_match) &&
  all(order_checks$max_abs_W_difference <= 1e-12)

raw_path <- file.path(pilot_dir, "tail_expansion_raw.csv")
components_path <- file.path(pilot_dir, "tail_expansion_components.csv")
summary_path <- file.path(pilot_dir, "tail_expansion_summary.csv")
report_path <- file.path(pilot_dir, "tail_expansion_report.txt")
write.csv(tail_raw, raw_path, row.names = FALSE)
write.csv(components, components_path, row.names = FALSE)
write.csv(nested_summary, summary_path, row.names = FALSE)

estimators <- c("Oracle-GMM", "Full-GMM", "DML-IVQR")
plot_files <- character(0)
colors <- grDevices::rainbow(N_REPLICATIONS)
for (tau in TAUS) {
  for (estimator in estimators) {
    selected <- tail_raw$tau == tau & tail_raw$estimator == estimator
    plot_data <- tail_raw[selected, ]
    safe_estimator <- tolower(gsub("-", "_", estimator))
    safe_tau <- gsub("\\.", "p", format(tau, nsmall = 2))
    for (view in c("full", "zoom")) {
      x_limits <- if (view == "full") c(-40, 40) else c(-10, 10)
      view_data <- plot_data[
        plot_data$alpha >= x_limits[1] & plot_data$alpha <= x_limits[2],
      ]
      finite_W <- view_data$W[is.finite(view_data$W)]
      y_limits <- range(c(finite_W, CRITICAL_VALUE), finite = TRUE)
      plot_path <- file.path(
        plot_dir,
        paste0(safe_estimator, "_tau_", safe_tau, "_", view, ".pdf")
      )
      plot_files <- c(plot_files, plot_path)
      grDevices::pdf(plot_path, width = 8, height = 5.5)
      plot(
        x_limits, y_limits, type = "n", xlab = "alpha", ylab = "W_N(alpha)",
        main = paste(estimator, "tau =", format(tau, nsmall = 2),
                     if (view == "full") "[-40,40]" else "[-10,10]")
      )
      for (replication in seq_len(N_REPLICATIONS)) {
        one <- view_data[view_data$replication == replication, ]
        lines(one$alpha, one$W, col = colors[replication], lwd = 1)
      }
      abline(h = CRITICAL_VALUE, lty = 2, lwd = 2)
      legend(
        "topright", legend = paste("rep", seq_len(N_REPLICATIONS)),
        col = colors, lty = 1, cex = 0.65, ncol = 2, bty = "n"
      )
      grDevices::dev.off()
    }
  }
}

range_text <- function(x) {
  finite <- x[is.finite(x)]
  if (!length(finite)) return("NA to NA")
  paste(format(min(finite), digits = 16), "to", format(max(finite), digits = 16))
}

full_summary <- nested_summary[nested_summary$window == "[-40,40]", ]
group_rows <- vector("list", length(TAUS) * length(estimators))
group_index <- 0L
for (tau in TAUS) {
  for (estimator in estimators) {
    group_index <- group_index + 1L
    s <- full_summary[
      full_summary$tau == tau & full_summary$estimator == estimator,
    ]
    group_rows[[group_index]] <- data.frame(
      tau = tau,
      estimator = estimator,
      n_profiles = nrow(s),
      n_left_endpoint_accepted = sum(s$left_endpoint_accepted %in% TRUE),
      n_right_endpoint_accepted = sum(s$right_endpoint_accepted %in% TRUE),
      n_any_left_outer = sum(s$any_accept_left_outer),
      n_any_right_outer = sum(s$any_accept_right_outer),
      n_all_full_range = sum(s$all_accepted),
      n_multiple_components = sum(s$n_accepted_components_grid > 1L),
      max_components = max(s$n_accepted_components_grid),
      alpha_hat_range = range_text(s$alpha_hat_grid),
      min_W_range = range_text(s$min_W),
      stringsAsFactors = FALSE
    )
  }
}
group_summary <- do.call(rbind, group_rows)

boundary_rows <- vector("list", length(TAUS) * length(estimators) * length(WINDOWS))
boundary_index <- 0L
for (tau in TAUS) {
  for (estimator in estimators) {
    for (window_name in names(WINDOWS)) {
      boundary_index <- boundary_index + 1L
      s <- nested_summary[
        nested_summary$tau == tau &
          nested_summary$estimator == estimator &
          nested_summary$window == window_name,
      ]
      boundary_rows[[boundary_index]] <- data.frame(
        tau = tau, estimator = estimator, window = window_name,
        left_boundary_count = sum(s$left_endpoint_accepted %in% TRUE),
        right_boundary_count = sum(s$right_endpoint_accepted %in% TRUE),
        all_accepted_count = sum(s$all_accepted),
        stringsAsFactors = FALSE
      )
    }
  }
}
boundary_summary <- do.call(rbind, boundary_rows)

report <- c(
  "Wide-tail expansion diagnostic",
  "",
  paste("R version:", R.version.string),
  paste("Diagnostic seed:", DIAGNOSTIC_SEED),
  paste("n:", N),
  paste("kappa:", format(KAPPA, nsmall = 2)),
  paste("taus:", paste(format(TAUS, nsmall = 2), collapse = ", ")),
  paste("replications:", N_REPLICATIONS),
  paste("tail grid: seq(-40,40,by=0.2); points:", length(TAIL_GRID)),
  paste("critical value:", format(CRITICAL_VALUE, digits = 16)),
  "",
  "Previous [-3,5] pilot reproduction:",
  paste("rows compared:", nrow(reproduced)),
  paste("key columns exact:", key_match),
  paste("statuses exact:", status_match),
  paste("acceptance indicators exact:", accepted_match),
  paste("W vector exact after CSV round trip:", reproduction_exact_W),
  paste("max absolute W difference:",
        format(reproduction_max_abs_W, scientific = TRUE, digits = 16)),
  paste("reproduction tolerance:", format(REPRODUCTION_TOLERANCE, scientific = TRUE)),
  paste("reproduction result:", if (reproduction_pass) "PASS" else "FAIL"),
  "",
  paste("tail evaluations:", nrow(tail_raw)),
  paste("alpha-level failures:", sum(is.na(tail_raw$W))),
  paste("lambda_bc generations:", lambda_generation_count, "(expected 30)"),
  paste("all lambda reuse checks:", all(lambda_reuse_checks)),
  "",
  "Small reversed-grid check:",
  capture.output(print(order_checks, row.names = FALSE, digits = 16)),
  "",
  "Full-range and outer-tail diagnostics by tau and estimator:",
  capture.output(print(group_summary, row.names = FALSE, right = FALSE)),
  "",
  "Boundary contact across nested windows:",
  capture.output(print(boundary_summary, row.names = FALSE, right = FALSE)),
  "",
  paste("Total discrete accepted components:", nrow(components)),
  paste("Profiles with multiple components:",
        sum(full_summary$n_accepted_components_grid > 1L)),
  paste("Maximum components in one profile:",
        max(full_summary$n_accepted_components_grid)),
  "",
  "Interpretation labels used:",
  "no acceptance near search boundary",
  "tail unresolved within search range",
  "Rejected endpoints alone are not treated as proof of mathematical boundedness.",
  "",
  "Plot files:",
  plot_files,
  "",
  paste("Overall validation:", if (validation_pass) "PASS" else "FAIL"),
  "No adaptive expansion, interpolation, refinement, final CR length, coverage,",
  "or power calculation was performed."
)
writeLines(report, report_path)
cat(paste(report, collapse = "\n"), "\n")

if (!validation_pass) {
  stop("Tail expansion validation failed; see ", report_path, ".")
}
