# Frozen thesis simulation and inference specification.
SAMPLE_SIZES <- c(500, 1000)

TAUS <- c(
  0.10,
  0.25,
  0.50,
  0.75,
  0.90
)

PARAMETER_LOWER <- -1
PARAMETER_UPPER <- 3

POINT_GRID <- seq(-1, 3, by = 0.10)
CR_GRID <- seq(-1, 3, by = 0.05)
CRITICAL_VALUE <- qchisq(0.95, df = 2)
POWER_DELTAS <- c(-0.50, -0.25, 0.25, 0.50)

alpha_true <- function(tau) {
  1 + qnorm(tau)
}

# Frozen performance definitions:
#   bias          = estimate - true_value
#   MAE term      = abs(estimate - true_value)
#   squared error = (estimate - true_value)^2
#   RMSE later    = sqrt(mean(squared_error)) across replications
#   coverage      = W_N(alpha_true) <= CRITICAL_VALUE
#   rejection     = W_N(alpha_true) > CRITICAL_VALUE
#   power(Delta)  = W_N(alpha_true + Delta) > CRITICAL_VALUE
# Coverage and power are evaluated directly at their exact requested alpha,
# never approximated with CR_GRID.
#
# Point estimates minimize W_N over POINT_GRID. R's which.min() chooses the
# first/smallest grid alpha when minima are exactly tied.
# The numerical CR is the accepted subset of CR_GRID; it is not expanded,
# smoothed, interpolated, or locally refined.

truths <- vapply(TAUS, alpha_true, numeric(1))
power_alternatives <- as.vector(outer(truths, POWER_DELTAS, "+"))
if (any(power_alternatives < PARAMETER_LOWER |
        power_alternatives > PARAMETER_UPPER)) {
  stop("One or more frozen power alternatives lie outside [-1,3].")
}
rm(truths, power_alternatives)
