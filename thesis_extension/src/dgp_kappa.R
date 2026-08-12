# Construct only the treatment component of the thesis DGP extension.
make_treatment_kappa <- function(z1, z2, epsilon, w, kappa) {
  if (!is.numeric(kappa) || length(kappa) != 1L || is.na(kappa) ||
      !is.finite(kappa) || kappa < 0 || kappa > 1) {
    stop("kappa must be one finite numeric value in [0, 1].")
  }

  input_lengths <- c(length(z1), length(z2), length(epsilon), length(w))
  if (length(unique(input_lengths)) != 1L) {
    stop("z1, z2, epsilon, and w must have equal lengths.")
  }

  d_latent <- kappa * (z1 + z2) +
    sqrt(2 * (1 - kappa^2)) * w +
    epsilon

  list(
    d_latent = d_latent,
    D = pnorm(d_latent)
  )
}
