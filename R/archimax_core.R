# -----------------------------------------------------------------------------
# Clayton-Archimax building blocks
#
#   C(u_1,...,u_d) = psi( sum_j phi(u_j) * A( phi(u_1)/S, ..., phi(u_d)/S ) ),
#   S = sum_j phi(u_j),  psi(t) = (1 + theta t)^(-1/theta),  phi = psi^{-1}.
#
# Key identity (Chatelain et al., 2020, Sec. 3): for fixed omega in the unit
# simplex, xi(omega) = min_j phi(U_j) / omega_j satisfies
#   P(xi(omega) > x) = psi(x A(omega)),
# i.e. A(omega) is a pure scale parameter of xi(omega).  Every estimator in
# this package is built on that identity.
# -----------------------------------------------------------------------------

#' Clayton generator and its inverse
#'
#' @param t Non-negative numeric vector.
#' @param u Numeric vector in (0, 1].
#' @param theta Clayton parameter, `theta > 0`.
#' @return Numeric vector.
#' @name clayton
NULL

#' @rdname clayton
#' @export
clayton_psi <- function(t, theta) {
  check_theta(theta)
  (1 + theta * t)^(-1 / theta)
}

#' @rdname clayton
#' @export
clayton_phi <- function(u, theta) {
  check_theta(theta)
  (u^(-theta) - 1) / theta
}

#' Moments of the radial variable Z with survival function psi_theta
#'
#' `E(Z) = 1 / (1 - theta)` (finite only for `theta < 1`) and
#' `E(log Z) = digamma(1) - digamma(1/theta) - log(theta)`.
#'
#' @inheritParams clayton
#' @return A single numeric value.
#' @name clayton_moments
NULL

#' @rdname clayton_moments
#' @export
clayton_EZ <- function(theta) {
  check_theta(theta)
  if (theta >= 1) stop("E(Z) is finite only for theta < 1.", call. = FALSE)
  1 / (1 - theta)
}

#' @rdname clayton_moments
#' @export
clayton_ElogZ <- function(theta) {
  check_theta(theta)
  digamma(1) - digamma(1 / theta) - log(theta)
}

#' Pseudo-sample xi_i(omega) = min_j phi(U_ij) / omega_j
#'
#' Coordinates with `omega_j = 0` are dropped from the minimum, so a
#' bivariate face of the d-simplex reproduces the bivariate estimator.
#'
#' @param U n x d matrix of pseudo-observations in (0, 1).
#' @param omega Numeric vector of length d on the unit simplex.
#' @inheritParams clayton
#' @return Numeric vector of length n.
#' @export
xi_transform <- function(U, omega, theta) {
  U <- as.matrix(U)
  if (length(omega) != ncol(U)) {
    stop("length(omega) must equal ncol(U).", call. = FALSE)
  }
  if (any(omega < 0) || abs(sum(omega) - 1) > 1e-8) {
    stop("omega must lie on the unit simplex.", call. = FALSE)
  }
  out <- rep(Inf, nrow(U))
  for (j in which(omega > 0)) {
    out <- pmin(out, clayton_phi(U[, j], theta) / omega[j])
  }
  out
}

#' Logistic (Gumbel) Pickands dependence function
#'
#' `A(w) = (w^rho + (1 - w)^rho)^(1/rho)`, `rho >= 1`.
#'
#' @param w Numeric vector in \[0, 1\].
#' @param rho Dependence parameter (`rho = 1`: independence,
#'   `rho -> Inf`: comonotonicity).
#' @export
pickands_logistic <- function(w, rho) {
  stopifnot(rho >= 1)
  (w^rho + (1 - w)^rho)^(1 / rho)
}

#' Simulate from a bivariate Clayton-Archimax copula with logistic A
#'
#' Exact sampling by numerical inversion of the conditional distribution
#' `C(v | u) = dC(u, v) / du`.
#'
#' @param n Sample size.
#' @inheritParams clayton
#' @param rho Logistic Pickands parameter, `rho >= 1`.
#' @return A data frame with columns `u`, `v`.
#' @export
rclayton_archimax <- function(n, theta, rho) {
  check_theta(theta)
  stopifnot(rho >= 1, n >= 1)
  u <- stats::runif(n)
  p <- stats::runif(n)
  v <- numeric(n)
  eps <- 1e-12
  for (i in seq_len(n)) {
    x1 <- clayton_phi(u[i], theta)
    cond_cdf <- function(vv) {
      x2 <- clayton_phi(vv, theta)
      ell <- (x1^rho + x2^rho)^(1 / rho)
      (1 + theta * ell)^(-1 / theta - 1) * x1^(rho - 1) * ell^(1 - rho) *
        u[i]^(-theta - 1)
    }
    v[i] <- stats::uniroot(function(vv) cond_cdf(vv) - p[i],
                           interval = c(eps, 1 - eps), tol = 1e-10)$root
  }
  data.frame(u = u, v = v)
}

#' Rank-based pseudo-observations, rank / (n + 1), column-wise
#'
#' @param X Numeric matrix or data frame (no missing values).
#' @export
pseudo_obs <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  U <- apply(X, 2, function(x) rank(x, ties.method = "average") / (n + 1))
  matrix(as.numeric(U), nrow = n, dimnames = dimnames(X))
}

#' Pseudo-observations re-ranked within groups
#'
#' Conditional copula inference given a discrete covariate: ranks are
#' computed separately inside each group so that every subgroup is a
#' proper rank sample (marginal shifts between groups do not leak into
#' the dependence estimate).
#'
#' @param X Numeric matrix (no missing values).
#' @param group Vector of group labels, length `nrow(X)`.
#' @export
pseudo_obs_by_group <- function(X, group) {
  X <- as.matrix(X)
  U <- matrix(NA_real_, nrow(X), ncol(X), dimnames = dimnames(X))
  for (g in unique(group)) {
    i <- which(group == g)
    U[i, ] <- pseudo_obs(X[i, , drop = FALSE])
  }
  U
}

check_theta <- function(theta) {
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta) ||
      theta <= 0) {
    stop("theta must be a single positive finite number.", call. = FALSE)
  }
  invisible(TRUE)
}
