# Exact W_N(alpha) profiles for the authors' replicated Monte Carlo criteria.
# Successful evaluations intentionally preserve the original implementation,
# including dnorm(e, mean(e), var(e)).

.finish_wn_profile <- function(grid, W, status_by_alpha) {
  successful <- which(!is.na(W))
  if (length(successful)) {
    minimizer <- successful[which.min(W[successful])]
    alpha_hat <- grid[minimizer]
    min_W <- W[minimizer]
  } else {
    alpha_hat <- NA_real_
    min_W <- NA_real_
  }
  list(
    grid = grid,
    W = W,
    alpha_hat = alpha_hat,
    min_W = min_W,
    status_by_alpha = data.frame(
      alpha = grid,
      status = status_by_alpha,
      stringsAsFactors = FALSE
    )
  )
}

.evaluate_profile_alpha <- function(alpha, evaluator) {
  alpha_warnings <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      evaluator(alpha),
      warning = function(warning_condition) {
        alpha_warnings <<- c(alpha_warnings, conditionMessage(warning_condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error_condition) error_condition
  )
  if (inherits(value, "error")) {
    return(list(W = NA_real_, status = paste0("ERROR: ", conditionMessage(value))))
  }
  value <- as.numeric(value)
  if (length(value) != 1L || !is.finite(value)) {
    return(list(W = NA_real_, status = "ERROR: criterion was not one finite value"))
  }
  status <- if (length(alpha_warnings)) {
    paste0("OK; WARNING: ", paste(unique(alpha_warnings), collapse = " | "))
  } else {
    "OK"
  }
  list(W = value, status = status)
}

.gmm_wn_profile <- function(y, D, X, Z, tau, grid) {
  W <- rep(NA_real_, length(grid))
  status <- rep("NOT RUN", length(grid))
  for (i in seq_along(grid)) {
    evaluated <- .evaluate_profile_alpha(grid[i], function(alpha) {
      beta <- rq(y - (alpha * D) ~ X, tau = tau)
      beta <- matrix(beta$coefficients, nrow = 1)
      e <- y - alpha * D - cbind(1, X) %*% t(beta)
      distribition <- c(dnorm(e, mean(e), var(e)))
      distribition <- diag(distribition)
      M <- t(Z) %*% distribition %*% X
      J <- t(X) %*% distribition %*% X
      delta <- M %*% solve(J)
      psi <- t(Z) - delta %*% t(X)
      indicator <- ifelse(e <= 0, 1, 0)
      g <- psi %*% (tau - indicator)
      invsigma <- solve(
        psi %*% diag(diag((tau - indicator) %*% t(tau - indicator))) %*% t(psi)
      )
      t(g) %*% invsigma %*% g
    })
    W[i] <- evaluated$W
    status[i] <- evaluated$status
  }
  .finish_wn_profile(grid, W, status)
}

oracle_wn_profile <- function(y, D, X10, Z, tau, grid) {
  .gmm_wn_profile(y, D, X10, Z, tau, grid)
}

full_wn_profile <- function(y, D, X100, Z, tau, grid) {
  .gmm_wn_profile(y, D, X100, Z, tau, grid)
}

dml_wn_profile <- function(y, D, X100, Z, tau, grid) {
  W <- rep(NA_real_, length(grid))
  status <- rep("NOT RUN", length(grid))
  for (i in seq_along(grid)) {
    evaluated <- .evaluate_profile_alpha(grid[i], function(alpha) {
      lasso <- cv.hqreg(
        X100,
        y - alpha * D,
        method = c("quantile"),
        tau = tau,
        FUN = c("hqreg"),
        nfolds = 5,
        type.measure = c("mae")
      )
      cv.beta <- as.matrix(lasso$fit$beta)
      kfold <- which(lasso$lambda == lasso$lambda.min, arr.ind = TRUE)
      kfold.beta <- cv.beta[, kfold]
      beta <- matrix(kfold.beta, nrow = 1)
      e <- y - alpha * D - cbind(1, X100) %*% t(beta)
      distribition <- c(dnorm(e, mean(e), var(e)))
      distribition <- diag(distribition)
      distribition <- sqrt(distribition)
      psi <- matrix(0, nrow = length(Z[1, ]), ncol = length(Z[, 1]))
      for (j in seq_len(length(Z[1, ]))) {
        delta <- rlasso(
          distribition %*% Z[, j] ~ distribition %*% X100,
          post = FALSE
        )
        delta <- matrix(delta$coefficients, ncol = 1)
        delta <- Z[, j] - cbind(1, X100) %*% delta
        psi[j, ] <- t(delta)
      }
      indicator <- ifelse(e <= 0, 1, 0)
      g <- psi %*% (tau - indicator)
      invsigma <- solve(
        psi %*% diag(diag((tau - indicator) %*% t(tau - indicator))) %*% t(psi)
      )
      t(g) %*% invsigma %*% g
    })
    W[i] <- evaluated$W
    status[i] <- evaluated$status
  }
  .finish_wn_profile(grid, W, status)
}

# Deterministic-CV variant. This is identical to dml_wn_profile() except that
# cv.hqreg receives seed = 2021, following the authors' inference/example code.
dml_wn_profile_fixedcv <- function(y, D, X100, Z, tau, grid) {
  W <- rep(NA_real_, length(grid))
  status <- rep("NOT RUN", length(grid))
  for (i in seq_along(grid)) {
    evaluated <- .evaluate_profile_alpha(grid[i], function(alpha) {
      lasso <- cv.hqreg(
        X100,
        y - alpha * D,
        method = c("quantile"),
        tau = tau,
        FUN = c("hqreg"),
        nfolds = 5,
        type.measure = c("mae"),
        seed = 2021
      )
      cv.beta <- as.matrix(lasso$fit$beta)
      kfold <- which(lasso$lambda == lasso$lambda.min, arr.ind = TRUE)
      kfold.beta <- cv.beta[, kfold]
      beta <- matrix(kfold.beta, nrow = 1)
      e <- y - alpha * D - cbind(1, X100) %*% t(beta)
      distribition <- c(dnorm(e, mean(e), var(e)))
      distribition <- diag(distribition)
      distribition <- sqrt(distribition)
      psi <- matrix(0, nrow = length(Z[1, ]), ncol = length(Z[, 1]))
      for (j in seq_len(length(Z[1, ]))) {
        delta <- rlasso(
          distribition %*% Z[, j] ~ distribition %*% X100,
          post = FALSE
        )
        delta <- matrix(delta$coefficients, ncol = 1)
        delta <- Z[, j] - cbind(1, X100) %*% delta
        psi[j, ] <- t(delta)
      }
      indicator <- ifelse(e <= 0, 1, 0)
      g <- psi %*% (tau - indicator)
      invsigma <- solve(
        psi %*% diag(diag((tau - indicator) %*% t(tau - indicator))) %*% t(psi)
      )
      t(g) %*% invsigma %*% g
    })
    W[i] <- evaluated$W
    status[i] <- evaluated$status
  }
  .finish_wn_profile(grid, W, status)
}

# Belloni-Chernozhukov pivotal penalty used by the authors' Table 3 code.
# The returned vector includes the intercept penalty followed by the slope
# penalties. For a fixed dataset and tau, callers should compute this once and
# pass the same vector to every profile evaluation (and every kappa value that
# shares X and tau).
bc_pivotal_lambda <- function(X, R = 1000, tau = 0.5, c = 2, alpha = 0.1) {
  norm2n <- function(z) sqrt(mean(z^2))
  n <- nrow(X)
  sigs <- apply(X, 2, norm2n)
  U <- matrix(runif(n * R), n)
  R <- (t(X) %*% (tau - (U < tau))) /
    (sigs * sqrt(tau * (1 - tau)))
  r <- apply(abs(R), 2, max)
  c * quantile(r, 1 - alpha) * sqrt(tau * (1 - tau)) * c(1, sigs)
}

# Grid-invariant Table 3 Belloni-Chernozhukov DML profile. The pivotal penalty
# is supplied by the caller and is never recomputed inside the alpha loop.
dml_wn_profile_bc <- function(y, D, X100, Z, tau, grid, lambda_bc) {
  if (!is.numeric(lambda_bc) || length(lambda_bc) != ncol(X100) + 1L ||
      any(!is.finite(lambda_bc))) {
    stop("lambda_bc must be a finite numeric vector of length ncol(X100) + 1.")
  }
  W <- rep(NA_real_, length(grid))
  status <- rep("NOT RUN", length(grid))
  for (i in seq_along(grid)) {
    evaluated <- .evaluate_profile_alpha(grid[i], function(alpha) {
      lasso <- quantreg::rq(
        y - alpha * D ~ X100,
        tau = tau,
        method = "lasso",
        lambda = lambda_bc
      )
      beta <- matrix(lasso$coefficients, ncol = 1)
      e <- y - alpha * D - cbind(1, X100) %*% beta
      distribition <- c(dnorm(e, mean(e), var(e)))
      distribition <- diag(distribition)
      distribition <- sqrt(distribition)
      psi <- matrix(0, nrow = length(Z[1, ]), ncol = length(Z[, 1]))
      for (j in seq_len(length(Z[1, ]))) {
        delta <- hdm::rlasso(
          distribition %*% Z[, j] ~ distribition %*% X100,
          post = FALSE
        )
        delta <- matrix(delta$coefficients, ncol = 1)
        delta <- Z[, j] - cbind(1, X100) %*% delta
        psi[j, ] <- t(delta)
      }
      indicator <- ifelse(e <= 0, 1, 0)
      g <- psi %*% (tau - indicator)
      invsigma <- solve(
        psi %*% diag(diag((tau - indicator) %*% t(tau - indicator))) %*% t(psi)
      )
      t(g) %*% invsigma %*% g
    })
    W[i] <- evaluated$W
    status[i] <- evaluated$status
  }
  result <- .finish_wn_profile(grid, W, status)
  result$lambda_bc <- lambda_bc
  result
}
