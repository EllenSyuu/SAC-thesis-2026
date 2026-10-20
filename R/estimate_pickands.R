# -----------------------------------------------------------------------------
# Estimators of the Pickands dependence function A(omega) of a
# Clayton-Archimax copula, given theta.
#
# Since xi(omega) =d Z / A(omega), with Z ~ psi_theta:
#   Pickands-type : A = E(Z) / mean(xi)
#   CFG-type      : A = exp( E(log Z) - mean(log xi) )
#   MLE (new)     : A = argmax_a  n log a - (1/theta + 1) sum log(1 + theta a xi)
#                   over a in [max(omega), 1]  -> always within the bounds.
# -----------------------------------------------------------------------------

negloglik_a <- function(a, xi, theta) {
  term <- 1 + theta * a * xi
  if (any(!is.finite(term)) || any(term <= 0)) return(1e12)
  -(length(xi) * log(a) - (1 / theta + 1) * sum(log(term)))
}

mle_scale <- function(xi, theta, lower) {
  if (lower >= 1) return(1)
  stats::optimize(negloglik_a, interval = c(lower, 1), xi = xi,
                  theta = theta, tol = 1e-10)$minimum
}

#' Estimate A(omega) at a single point of the unit simplex
#'
#' Works in any dimension d; `omega` with a single non-zero coordinate
#' returns 1 by definition.
#'
#' @param U n x d matrix of pseudo-observations.
#' @param omega Point on the unit simplex (length d).
#' @inheritParams clayton
#' @param method One of `"mle"`, `"cfg"`, `"pickands"`.
#' @return A single numeric value.
#' @export
pickands_at <- function(U, omega, theta,
                        method = c("mle", "cfg", "pickands")) {
  method <- match.arg(method)
  if (sum(omega > 0) == 1) return(1)
  U <- as.matrix(U)
  U <- U[stats::complete.cases(U), , drop = FALSE]
  xi <- xi_transform(U, omega, theta)
  switch(method,
    mle      = mle_scale(xi, theta, lower = max(omega)),
    cfg      = exp(clayton_ElogZ(theta) - mean(log(xi))),
    pickands = clayton_EZ(theta) / mean(xi)
  )
}

#' Bivariate Pickands dependence function on a grid
#'
#' @param U n x 2 matrix (or data frame) of pseudo-observations.
#' @param w Grid in \[0, 1\]; `A(w)` is evaluated at `omega = (w, 1 - w)`.
#' @inheritParams pickands_at
#' @param project Clip estimates into the admissible band
#'   `[max(w, 1 - w), 1]` (only relevant for CFG / Pickands).
#' @return Numeric vector, same length as `w`.
#' @export
pickands_curve <- function(U, theta, w = seq(0.01, 0.99, by = 0.01),
                           method = c("mle", "cfg", "pickands"),
                           project = FALSE) {
  method <- match.arg(method)
  U <- as.matrix(U)[, 1:2]
  A <- vapply(w, function(wk) pickands_at(U, c(wk, 1 - wk), theta, method),
              numeric(1))
  if (project) A <- pmin(pmax(A, pmax(w, 1 - w)), 1)
  A
}

#' Covariate-dependent MLE of the Pickands dependence function
#'
#' At each grid point `w`, fits `A(w | z_i) = g(z_i' beta(w))` by maximising
#' the Clayton-Archimax pseudo-likelihood
#' `sum_i log a_i - (1/theta + 1) log(1 + theta a_i xi_i(w))`.
#'
#' * `link = "identity"`: `a_i = z_i' beta`, maximised subject to
#'   `max(w, 1-w) <= a_i <= 1` for every observed covariate profile
#'   (multi-start Nelder-Mead plus a warm start from group-wise MLEs).
#' * `link = "logit"`:
#'   `a_i = m + (1 - m) plogis(z_i' beta)`, `m = max(w, 1-w)`; the bounds
#'   hold automatically and the problem is unconstrained (BFGS).
#'
#' @param U n x 2 matrix of pseudo-observations.
#' @param X n x p design matrix whose first column is the intercept.
#' @inheritParams clayton
#' @param w Grid in (0, 1).
#' @param link `"identity"` or `"logit"`.
#' @return A list with `beta` (length(w) x p), `A` (length(w) x G fitted
#'   values for the G unique covariate profiles), `profiles`, `w`.
#' @export
pickands_cov_mle <- function(U, X, theta, w = seq(0.01, 0.99, by = 0.01),
                             link = c("identity", "logit")) {
  link <- match.arg(link)
  check_theta(theta)
  U <- as.matrix(U)[, 1:2]
  X <- as.matrix(X)
  ok <- stats::complete.cases(U, X)
  U <- U[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]
  X_u <- unique(X)
  p <- ncol(X)

  beta <- matrix(NA_real_, length(w), p, dimnames = list(NULL, colnames(X)))
  prev <- NULL
  for (k in seq_along(w)) {
    xi <- xi_transform(U, c(w[k], 1 - w[k]), theta)
    lo <- max(w[k], 1 - w[k])
    beta[k, ] <- if (link == "identity") {
      fit_identity(xi, X, X_u, theta, lo, init = prev)
    } else {
      fit_logit(xi, X, theta, lo, init = prev)
    }
    prev <- beta[k, ]
  }
  A <- matrix(vapply(seq_along(w), function(k) {
    eta <- as.numeric(X_u %*% beta[k, ])
    if (link == "identity") return(eta)
    lo <- max(w[k], 1 - w[k])
    lo + (1 - lo) * stats::plogis(eta)
  }, numeric(nrow(X_u))), nrow = length(w), byrow = TRUE)
  colnames(A) <- apply(X_u, 1, paste, collapse = ",")
  list(beta = beta, A = A, profiles = as.data.frame(X_u), w = w, link = link)
}

fit_identity <- function(xi, X, X_u, theta, lo, init = NULL) {
  p <- ncol(X)
  nll <- function(b) {
    a <- as.numeric(X %*% b)
    a_u <- as.numeric(X_u %*% b)
    if (any(!is.finite(a_u)) || any(a_u < lo) || any(a_u > 1)) return(1e12)
    term <- 1 + theta * a * xi
    if (any(term <= 0)) return(1e12)
    -sum(log(a) - (1 / theta + 1) * log(term))
  }
  starts <- lapply(seq(lo + 0.01, 0.99, length.out = 5),
                   function(b0) c(b0, rep(0, p - 1)))
  # warm start: group-wise MLE mapped back to beta by least squares
  key <- apply(X, 1, paste, collapse = "_")
  key_u <- apply(X_u, 1, paste, collapse = "_")
  a_grp <- vapply(key_u, function(g) mle_scale(xi[key == g], theta, lo),
                  numeric(1))
  b_ls <- tryCatch(as.numeric(qr.solve(X_u, a_grp)), error = function(e) NULL)
  if (!is.null(b_ls) && all(is.finite(b_ls))) starts <- c(starts, list(b_ls))
  if (!is.null(init) && all(is.finite(init))) starts <- c(starts, list(init))

  best <- list(value = Inf, par = starts[[1]])
  for (s in starts) {
    opt <- tryCatch(stats::optim(s, nll, method = "Nelder-Mead",
                                 control = list(maxit = 5000)),
                    error = function(e) NULL)
    if (!is.null(opt) && opt$value < best$value) best <- opt
  }
  best$par
}

fit_logit <- function(xi, X, theta, lo, init = NULL) {
  p <- ncol(X)
  span <- 1 - lo
  nll <- function(b) {
    a <- lo + span * stats::plogis(as.numeric(X %*% b))
    term <- 1 + theta * a * xi
    if (any(!is.finite(a)) || any(a <= 0) || any(term <= 0)) return(1e12)
    -sum(log(a) - (1 / theta + 1) * log(term))
  }
  starts <- list(rep(0, p), c(-2, rep(0, p - 1)), c(2, rep(0, p - 1)))
  if (p > 1) {
    starts <- c(starts, list(c(0, -2, rep(0, p - 2)), c(0, 2, rep(0, p - 2))))
  }
  if (!is.null(init) && all(is.finite(init))) starts <- c(starts, list(init))
  best <- list(value = Inf, par = starts[[1]])
  for (s in starts) {
    opt <- tryCatch(stats::optim(s, nll, method = "BFGS",
                                 control = list(maxit = 2000)),
                    error = function(e) NULL)
    if (!is.null(opt) && opt$value < best$value) best <- opt
  }
  best$par
}

#' Points on the segment from a bivariate face to the simplex centroid
#'
#' `omega(t) = (1 - t) omega_jk + t omega_0`, where `omega_jk` puts mass 1/2
#' on coordinates j and k and `omega_0 = (1/d, ..., 1/d)`.
#'
#' @param j,k Indices of the starting face.
#' @param t Position in \[0, 1\].
#' @param d Dimension.
#' @export
simplex_path <- function(j, k, t, d) {
  start <- numeric(d)
  start[c(j, k)] <- 0.5
  (1 - t) * start + t * rep(1 / d, d)
}

#' Extremal coefficient eta = d * A(omega_0) at the simplex centroid
#'
#' `eta = 1` under complete dependence and `eta = d` under independence.
#'
#' @inheritParams pickands_at
#' @export
extremal_coefficient <- function(U, theta, method = "mle") {
  d <- ncol(U)
  d * pickands_at(U, rep(1 / d, d), theta, method)
}
