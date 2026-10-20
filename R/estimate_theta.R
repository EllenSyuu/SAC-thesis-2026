# -----------------------------------------------------------------------------
# Moment-based estimation of the Clayton generator parameter theta
# (Chatelain, Fougeres & Neslehova, 2020, Example 7.1, Eq. 7.4), with a
# delete-one jackknife standard error.
#
# Kendall-distribution moments m_k = E[W^k] satisfy
#   (1 - 2 m_1) / (1 - 3 m_2) = (theta + 3) / (2 (theta + 2)),
# which inverts to theta_n = S_n / R_n with
#   S_n = 8 m_{n,1} - 9 m_{n,2} - 1,   R_n = 1 - 4 m_{n,1} + 3 m_{n,2},
# where m_{n,1}, m_{n,2} are U-statistics (Ben Ghorbal et al., 2009).
# -----------------------------------------------------------------------------

concordance_matrix <- function(u, v) {
  M <- outer(u, u, "<=") & outer(v, v, "<=")
  diag(M) <- FALSE
  M * 1
}

#' Moment estimator of the Clayton-Archimax generator parameter
#'
#' @param u,v Pseudo-observations of one station pair (missing values are
#'   dropped pairwise).
#' @param alpha Significance level for the Wald interval.
#' @param se Logical; compute a delete-one jackknife standard error?
#' @return A list with `theta`, `se`, `ci` (length 2), `n`, `m1`, `m2`.
#' @export
estimate_theta <- function(u, v, alpha = 0.10, se = TRUE) {
  ok <- stats::complete.cases(u, v)
  u <- u[ok]
  v <- v[ok]
  n <- length(u)
  if (n < 4) stop("Need at least 4 complete pairs.", call. = FALSE)
  if (any(u <= 0 | u >= 1) || any(v <= 0 | v >= 1)) {
    stop("u and v must be pseudo-observations strictly inside (0, 1); ",
         "use pseudo_obs() on raw data first.", call. = FALSE)
  }

  M  <- concordance_matrix(u, v)
  A  <- rowSums(M)                    # #{j : X_i <= X_j}
  B  <- colSums(M)                    # #{i : X_i <= X_j}
  m1 <- sum(M) / (n * (n - 1))
  m2 <- sum(B * (B - 1)) / (n * (n - 1) * (n - 2))
  S  <- 8 * m1 - 9 * m2 - 1
  R  <- 1 - 4 * m1 + 3 * m2
  theta <- S / R

  out <- list(theta = theta, se = NA_real_, ci = c(NA_real_, NA_real_),
              n = n, m1 = m1, m2 = m2)
  if (!se) return(out)

  # Delete-one jackknife, O(n^2) via closed-form updates of the U-statistics
  MB <- as.numeric(M %*% B)
  tot1 <- sum(M)
  tot2 <- sum(B * (B - 1))
  m1_l <- (tot1 - A - B) / ((n - 1) * (n - 2))
  m2_l <- (tot2 - B * (B - 1) - 2 * (MB - A)) / ((n - 1) * (n - 2) * (n - 3))
  th_l <- (8 * m1_l - 9 * m2_l - 1) / (1 - 4 * m1_l + 3 * m2_l)
  out$se <- sqrt((n - 1) / n * sum((th_l - mean(th_l))^2))
  out$ci <- theta + c(-1, 1) * stats::qnorm(1 - alpha / 2) * out$se
  out
}

#' Pairwise and global theta estimates for a d-variate sample
#'
#' The global estimate is the average of the `choose(d, 2)` pairwise
#' estimates, as recommended by Chatelain et al. (2020, Sec. 7).
#'
#' @param U n x d matrix of pseudo-observations (may contain `NA`).
#' @inheritParams estimate_theta
#' @return A list: `pairs` (data frame) and `global` (numeric).
#' @export
estimate_theta_pairs <- function(U, alpha = 0.10) {
  U <- as.matrix(U)
  nm <- colnames(U)
  if (is.null(nm)) nm <- paste0("X", seq_len(ncol(U)))
  idx <- utils::combn(ncol(U), 2)
  rows <- lapply(seq_len(ncol(idx)), function(k) {
    i <- idx[1, k]
    j <- idx[2, k]
    e <- estimate_theta(U[, i], U[, j], alpha = alpha)
    data.frame(station_1 = nm[i], station_2 = nm[j], n = e$n,
               theta = e$theta, se = e$se, ci_lower = e$ci[1],
               ci_upper = e$ci[2])
  })
  pairs <- do.call(rbind, rows)
  list(pairs = pairs, global = mean(pairs$theta))
}
