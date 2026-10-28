# -----------------------------------------------------------------------------
# Year-block bootstrap
#
# Monthly maxima within a year share large-scale climate conditions, and the
# six stations are observed simultaneously. Resampling whole calendar years
# (with replacement) keeps both the cross-station dependence and the
# within-year seasonal structure intact.
# -----------------------------------------------------------------------------

#' Year-block bootstrap of an arbitrary statistic
#'
#' Replicate `b` resamples the distinct values of `data[[year_col]]` with
#' replacement and evaluates `statistic` on the stacked rows. A per-replicate
#' seed is drawn up front from `seed`, so results are identical for any
#' number of cores. Replicates that error are kept as `NULL`, counted, and
#' reported in a warning - never silently dropped.
#'
#' @param data Data frame with a year column.
#' @param statistic Function of a resampled data frame; returns any object.
#' @param B Number of replicates.
#' @param seed Integer seed.
#' @param year_col Name of the block column.
#' @param cores Number of forked workers (`parallel::mclapply`; 1 on Windows).
#' @return A list with `replicates` (length-B list), `n_failed`, `errors`,
#'   `seeds`, `B`.
#' @examples
#' set.seed(1)
#' d <- data.frame(year = rep(1:30, each = 4), x = stats::rnorm(120))
#' b <- boot_years(d, function(dat) mean(dat$x), B = 50, seed = 1)
#' b$n_failed
#' boot_ci(unlist(b$replicates))
#' @export
boot_years <- function(data, statistic, B = 200, seed = 1, year_col = "year",
                       cores = 1) {
  if (!year_col %in% names(data)) stop("year_col not found in data.", call. = FALSE)
  if (B < 1) stop("B must be >= 1.", call. = FALSE)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  rows_by_year <- split(seq_len(nrow(data)), data[[year_col]])
  years <- names(rows_by_year)
  one <- function(b) {
    set.seed(seeds[b])
    pick <- sample(years, length(years), replace = TRUE)
    rb <- data[unlist(rows_by_year[pick], use.names = FALSE), , drop = FALSE]
    tryCatch(statistic(rb), error = function(e) {
      structure(conditionMessage(e), class = "boot_error")
    })
  }
  reps <- if (cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_len(B), one, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(B), one)
  }
  bad <- vapply(reps, function(x) inherits(x, "boot_error") || inherits(x, "try-error"),
                logical(1))
  errors <- unique(vapply(reps[bad], as.character, character(1)))
  if (any(bad)) {
    warning(sprintf("boot_years: %d of %d replicates failed (first error: %s)",
                    sum(bad), B, errors[1]), call. = FALSE)
    reps[bad] <- list(NULL)
  }
  list(replicates = reps, n_failed = sum(bad), errors = errors, seeds = seeds, B = B)
}

#' Percentile bootstrap interval
#'
#' @param x Numeric vector of bootstrap replicates (`NA` allowed).
#' @param level Confidence level.
#' @return Named vector `lower`, `upper`, `n_ok`.
#' @examples
#' boot_ci(stats::rnorm(1000))
#' @export
boot_ci <- function(x, level = 0.90) {
  a <- (1 - level) / 2
  ok <- is.finite(x)
  if (sum(ok) < 2) return(c(lower = NA_real_, upper = NA_real_, n_ok = sum(ok)))
  q <- stats::quantile(x[ok], c(a, 1 - a), names = FALSE, type = 7)
  c(lower = q[1], upper = q[2], n_ok = sum(ok))
}

#' Rank columns into pseudo-observations, keeping missing values
#'
#' Column-wise `rank / (n_j + 1)` where `n_j` counts the non-missing values
#' of column j. Unlike [pseudo_obs()] this accepts `NA`.
#'
#' @param X Numeric matrix or data frame.
#' @return Numeric matrix with the dimensions of `X`; `NA` stays `NA`.
#' @examples
#' rank_columns(cbind(c(1, NA, 3), c(2, 1, 3)))
#' @export
rank_columns <- function(X) {
  X <- as.matrix(X)
  U <- apply(X, 2, function(x) {
    rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  })
  matrix(U, nrow = nrow(X), dimnames = dimnames(X))
}
