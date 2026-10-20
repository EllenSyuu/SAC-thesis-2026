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
  if (!all(is.finite(term)) || any(term <= 0)) return(1e12)
  -(length(xi) * log(a) - (1 / theta + 1) * sum(log(term)))
}

mle_scale <- function(xi, theta, lower) {
  if (lower >= 1) return(1)
  opt <- stats::optimize(negloglik_a, interval = c(lower, 1), xi = xi,
                         theta = theta, tol = 1e-10)
  if (!is.finite(opt$objective) || opt$objective >= 1e12) return(NA_real_)
  opt$minimum
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
#' @return A single numeric value; `NA` (with a warning) if the estimate is
#'   not finite, e.g. because `U` contains 0 or 1.
#' @export
pickands_at <- function(U, omega, theta,
                        method = c("mle", "cfg", "pickands")) {
  method <- match.arg(method)
  if (sum(omega > 0) == 1) return(1)
  U <- as.matrix(U)
  U <- U[stats::complete.cases(U), , drop = FALSE]
  xi <- xi_transform(U, omega, theta)
  out <- switch(method,
    mle      = mle_scale(xi, theta, lower = max(omega)),
    cfg      = exp(clayton_ElogZ(theta) - mean(log(xi))),
    pickands = clayton_EZ(theta) / mean(xi)
  )
  if (!is.finite(out)) {
    warning("pickands_at: non-finite estimate (check that U lies strictly in (0, 1)).",
            call. = FALSE)
    out <- NA_real_
  }
  out
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
#' Failure is never silent: if no start reaches a feasible finite optimum at
#' some `w`, that row of `beta` / `A` is `NA`, `diagnostics$status` is
#' `"failed"` and a single warning summarises the affected grid points.
#'
#' @param U n x 2 matrix of pseudo-observations.
#' @param X n x p design matrix whose first column is the intercept.
#' @inheritParams clayton
#' @param w Grid in (0, 1).
#' @param link `"identity"` or `"logit"`.
#' @return A list with `beta` (length(w) x p), `A` (length(w) x G fitted
#'   values for the G unique covariate profiles), `profiles`, `w`, `link`
#'   and `diagnostics`: a data frame with one row per `w` giving `status`
#'   (`"ok"`, `"maxit"` or `"failed"`), the minimised negative
#'   log-likelihood `nll`, the number of starts that returned a feasible
#'   value `n_feasible_starts`, and `at_bound` (some fitted profile lies on
#'   the boundary of `[max(w, 1-w), 1]`).
#' @export
pickands_cov_mle <- function(U, X, theta, w = seq(0.01, 0.99, by = 0.01),
                             link = c("identity", "logit")) {
  link <- match.arg(link)
  check_theta(theta)
  U <- as.matrix(U)[, 1:2, drop = FALSE]
  X <- as.matrix(X)
  if (nrow(U) != nrow(X)) stop("U and X must have the same number of rows.", call. = FALSE)
  ok <- stats::complete.cases(U, X)
  U <- U[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]
  if (nrow(U) < 2) stop("Need at least 2 complete observations.", call. = FALSE)
  X_u <- unique(X)
  p <- ncol(X)

  beta <- matrix(NA_real_, length(w), p, dimnames = list(NULL, colnames(X)))
  diag <- data.frame(w = w, status = NA_character_, nll = NA_real_,
                     n_feasible_starts = NA_integer_, at_bound = NA)
  prev <- NULL
  for (k in seq_along(w)) {
    xi <- xi_transform(U, c(w[k], 1 - w[k]), theta)
    lo <- max(w[k], 1 - w[k])
    fit <- if (link == "identity") {
      fit_identity(xi, X, X_u, theta, lo, init = prev)
    } else {
      fit_logit(xi, X, theta, lo, init = prev)
    }
    beta[k, ] <- fit$par
    diag$status[k] <- fit$status
    diag$nll[k] <- fit$nll
    diag$n_feasible_starts[k] <- fit$n_feasible
    if (fit$status != "failed") prev <- fit$par
  }
  A <- matrix(vapply(seq_along(w), function(k) {
    eta <- as.numeric(X_u %*% beta[k, ])
    if (link == "identity") return(eta)
    lo <- max(w[k], 1 - w[k])
    lo + (1 - lo) * stats::plogis(eta)
  }, numeric(nrow(X_u))), nrow = length(w), byrow = TRUE)
  colnames(A) <- apply(X_u, 1, paste, collapse = ",")
  lo_w <- pmax(w, 1 - w)
  diag$at_bound <- apply(abs(A - lo_w) < 1e-6 | abs(A - 1) < 1e-6, 1, any)

  failed <- diag$status == "failed"
  if (any(failed)) {
    warning(sprintf("pickands_cov_mle: no feasible optimum at %d of %d grid points (w = %s); returned NA there.",
                    sum(failed), length(w),
                    paste(utils::head(signif(w[failed], 3), 5), collapse = ", ")),
            call. = FALSE)
  }
  list(beta = beta, A = A, profiles = as.data.frame(X_u), w = w, link = link,
       diagnostics = diag)
}

PENALTY <- 1e12

# Run optim from every start and keep the best feasible solution. Returns
# status "failed" (and NA parameters) if no start reaches a finite value
# below the infeasibility penalty.
multistart <- function(starts, nll, method, maxit, bounds_1d) {
  one_d <- length(starts[[1]]) == 1
  if (one_d) starts <- starts[1]    # Brent is deterministic on an interval
  best <- NULL
  n_feasible <- 0L
  for (s in starts) {
    opt <- tryCatch(
      if (one_d) {
        stats::optim(s, nll, method = "Brent", lower = bounds_1d[1],
                     upper = bounds_1d[2])
      } else {
        stats::optim(s, nll, method = method, control = list(maxit = maxit))
      },
      error = function(e) NULL)
    if (is.null(opt) || !is.finite(opt$value) || opt$value >= PENALTY) next
    n_feasible <- n_feasible + 1L
    if (is.null(best) || opt$value < best$value) best <- opt
  }
  if (is.null(best)) {
    return(list(par = rep(NA_real_, length(starts[[1]])), status = "failed",
                nll = NA_real_, n_feasible = 0L))
  }
  list(par = best$par, status = if (best$convergence == 0) "ok" else "maxit",
       nll = best$value, n_feasible = n_feasible)
}

fit_identity <- function(xi, X, X_u, theta, lo, init = NULL) {
  p <- ncol(X)
  nll <- function(b) {
    a <- as.numeric(X %*% b)
    a_u <- as.numeric(X_u %*% b)
    if (any(!is.finite(a_u)) || any(a_u < lo) || any(a_u > 1)) return(PENALTY)
    term <- 1 + theta * a * xi
    if (!all(is.finite(term)) || any(term <= 0)) return(PENALTY)
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
  multistart(starts, nll, "Nelder-Mead", 5000, bounds_1d = c(lo, 1))
}

fit_logit <- function(xi, X, theta, lo, init = NULL) {
  p <- ncol(X)
  span <- 1 - lo
  nll <- function(b) {
    a <- lo + span * stats::plogis(as.numeric(X %*% b))
    term <- 1 + theta * a * xi
    if (!all(is.finite(a)) || any(a <= 0) || !all(is.finite(term)) ||
        any(term <= 0)) return(PENALTY)
    -sum(log(a) - (1 / theta + 1) * log(term))
  }
  starts <- list(rep(0, p), c(-2, rep(0, p - 1)), c(2, rep(0, p - 1)))
  if (p > 1) {
    starts <- c(starts, list(c(0, -2, rep(0, p - 2)), c(0, 2, rep(0, p - 2))))
  }
  if (!is.null(init) && all(is.finite(init))) starts <- c(starts, list(init))
  multistart(starts, nll, "BFGS", 2000, bounds_1d = c(-30, 30))
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
