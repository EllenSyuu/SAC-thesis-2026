#' Model-implied tail-dependence coefficients of a Clayton-Archimax pair
#'
#' For the bivariate Clayton-Archimax copula (Chatelain et al., 2020,
#' Remark 2): `lambda_L = (2 A(1/2))^(-1/theta)` and
#' `lambda_U = 2 - 2 A(1/2)`.
#'
#' @param A_half Value(s) of the Pickands function at `w = 1/2`.
#' @inheritParams clayton
#' @return A data frame with columns `lambda_L`, `lambda_U`.
#' @export
tail_dependence_archimax <- function(A_half, theta) {
  check_theta(theta)
  data.frame(lambda_L = (2 * A_half)^(-1 / theta),
             lambda_U = 2 - 2 * A_half)
}

#' Empirical upper / lower tail-dependence coefficients at threshold q
#'
#' `chi_U(q) = P(V > q | U > q)` and `chi_L(q) = P(V < 1 - q | U < 1 - q)`,
#' computed with [extRemes::taildep()] on rank-transformed data.
#'
#' @param x,y Numeric vectors (missing values dropped pairwise).
#' @param q Threshold(s) in (0, 1).
#' @return A data frame with `q`, `lambda_U`, `lambda_L`.
#' @export
tail_dependence_empirical <- function(x, y, q = 0.8) {
  ok <- stats::complete.cases(x, y)
  U <- pseudo_obs(cbind(x[ok], y[ok]))
  chi <- function(a, b) {
    vapply(q, function(qq) {
      unname(extRemes::taildep(a, b, u = qq, type = "chi"))
    }, numeric(1))
  }
  data.frame(q = q, lambda_U = chi(U[, 1], U[, 2]),
             lambda_L = chi(-U[, 1], -U[, 2]))
}

#' Integrated squared error and integrated relative absolute error
#'
#' Trapezoidal approximations of
#' `ISE = int (A_hat - A)^2` and `IRAE = int |A_hat - A| / A`.
#'
#' @param w Grid.
#' @param A_hat,A_true Estimated and true values on `w`.
#' @name curve_errors
NULL

#' @rdname curve_errors
#' @export
ise <- function(w, A_hat, A_true) trapz(w, (A_hat - A_true)^2)

#' @rdname curve_errors
#' @export
irae <- function(w, A_hat, A_true) trapz(w, abs(A_hat - A_true) / A_true)

trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
