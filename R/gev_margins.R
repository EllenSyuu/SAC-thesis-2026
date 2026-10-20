# -----------------------------------------------------------------------------
# Data preparation and non-stationary GEV margins
# -----------------------------------------------------------------------------

#' Quebec seasons used in the analysis
#'
#' AM = Apr-May, JJ = Jun-Jul, AS = Aug-Sep, ON = Oct-Nov; other months
#' (snow-dominated) map to `NA`.
#'
#' @param month Integer month 1-12.
#' @return Factor with levels `AM`, `JJ`, `AS`, `ON`.
#' @export
qc_season <- function(month) {
  lab <- c(NA, NA, NA, "AM", "AM", "JJ", "JJ", "AS", "AS", "ON", "ON", NA)
  factor(lab[as.integer(month)], levels = c("AM", "JJ", "AS", "ON"))
}

#' Pre / post climate-period indicator
#'
#' @param year Integer year.
#' @param cutoff Last year of the reference period (default 1975).
#' @return Factor with levels `pre`, `post`.
#' @export
qc_period <- function(year, cutoff = 1975) {
  factor(ifelse(year <= cutoff, "pre", "post"), levels = c("pre", "post"))
}

#' Monthly block maxima from daily records
#'
#' @param daily Data frame with `station_name`, `year`, `month`,
#'   `total_precip`.
#' @return One row per station-month with at least one valid day:
#'   `station_name`, `year`, `month`, `max_precip`, `n_days`.
#' @export
monthly_maxima <- function(daily) {
  daily <- daily[!is.na(daily$total_precip),
                 c("station_name", "year", "month", "total_precip")]
  f <- total_precip ~ station_name + year + month
  mx <- stats::aggregate(f, data = daily, FUN = max)
  nd <- stats::aggregate(f, data = daily, FUN = length)
  names(mx)[4] <- "max_precip"
  names(nd)[4] <- "n_days"
  out <- merge(mx, nd, by = c("station_name", "year", "month"))
  out[order(out$station_name, out$year, out$month), ]
}

#' Fit a non-stationary GEV with covariate-dependent location / scale
#'
#' Thin wrapper around [extRemes::fevd()] with constant shape.
#'
#' @param data Data frame containing `max_precip` and the covariates.
#' @param location,scale One-sided formulas.
#' @return An `fevd` object.
#' @export
fit_gev <- function(data, location = ~ period + season, scale = ~ season) {
  extRemes::fevd(x = data$max_precip, data = data, type = "GEV",
                 method = "MLE", location.fun = location,
                 scale.fun = scale, shape.fun = ~ 1)
}

#' Coefficient table (estimate, SE, z, p, 95\% Wald CI) for an fevd fit
#'
#' @param fit An `fevd` object.
#' @export
gev_coef_table <- function(fit) {
  est <- fit$results$par
  se <- sqrt(diag(extRemes::parcov.fevd(fit)))
  z <- est / se
  data.frame(term = names(est), estimate = unname(est), se = unname(se),
             z = unname(z), p_value = 2 * stats::pnorm(-abs(unname(z))),
             ci_lower = unname(est - 1.96 * se),
             ci_upper = unname(est + 1.96 * se))
}

#' Information criteria for an fevd fit
#'
#' @param fit An `fevd` object.
#' @export
gev_ic <- function(fit) {
  nll <- fit$results$value
  k <- length(fit$results$par)
  n <- length(fit$x)
  data.frame(loglik = -nll, n_par = k, AIC = 2 * nll + 2 * k,
             BIC = 2 * nll + k * log(n))
}

#' Observation-wise GEV parameters of a fitted non-stationary model
#'
#' @param fit An `fevd` object.
#' @return Data frame with `mu`, `sigma`, `xi`, one row per observation.
#' @export
gev_params <- function(fit) {
  p <- extRemes::findpars(fit)
  data.frame(mu = p$location, sigma = p$scale, xi = p$shape)
}

#' Transform observations to the standard Gumbel scale
#'
#' `z = log(1 + xi (x - mu) / sigma) / xi`, so that
#' `exp(-exp(-z))` is approximately Uniform(0, 1).
#'
#' @param x Observations.
#' @param mu,sigma,xi GEV parameters (recycled).
#' @export
gev_to_gumbel <- function(x, mu, sigma, xi) {
  n <- max(length(x), length(mu), length(sigma), length(xi))
  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  xi <- rep_len(xi, n)
  s <- (x - mu) / sigma
  ifelse(abs(xi) < 1e-10, s, log1p(xi * s) / xi)
}

#' Mean of a GEV distribution (xi < 1)
#'
#' @inheritParams gev_to_gumbel
#' @export
gev_mean <- function(mu, sigma, xi) {
  n <- max(length(mu), length(sigma), length(xi))
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  xi <- rep_len(xi, n)
  ifelse(abs(xi) < 1e-10, mu - sigma * digamma(1),
         mu + sigma * (gamma(1 - xi) - 1) / xi)
}

#' Backward elimination for location / scale terms by likelihood-ratio tests
#'
#' Removes, one at a time, the least significant *removable* term (respecting
#' marginality: a main effect cannot be dropped while an interaction
#' containing it remains) until every remaining candidate has `p < alpha`.
#'
#' @param data Station data frame.
#' @param terms Character vector of candidate terms (shared start for both
#'   parameters).
#' @param order `"scale_first"` (thesis) or `"location_first"`.
#' @param alpha LRT level.
#' @return A list with final `location`, `scale` terms, the `fit`, and a
#'   `path` data frame recording each step.
#' @export
gev_backward <- function(data,
                         terms = c("anomaly", "period", "season",
                                   "anomaly:period", "anomaly:season",
                                   "period:season"),
                         order = c("scale_first", "location_first"),
                         alpha = 0.05) {
  order <- match.arg(order)
  loc <- terms
  sca <- terms
  path <- list()
  phases <- if (order == "scale_first") c("scale", "location") else
    c("location", "scale")
  for (ph in phases) {
    repeat {
      cur <- safe_fit(data, loc, sca)
      target <- if (ph == "scale") sca else loc
      cand <- removable_terms(target)
      if (length(cand) == 0 || is.null(cur)) break
      tests <- lapply(cand, function(tm) {
        new <- setdiff(target, tm)
        red <- if (ph == "scale") safe_fit(data, loc, new) else
          safe_fit(data, new, sca)
        if (is.null(red)) return(c(LR = NA, df = NA, p = NA))
        lr <- 2 * (red$results$value - cur$results$value)
        df <- length(cur$results$par) - length(red$results$par)
        c(LR = lr, df = df, p = stats::pchisq(lr, df, lower.tail = FALSE))
      })
      pv <- vapply(tests, function(t) t[["p"]], numeric(1))
      if (all(is.na(pv))) break
      k <- which.max(pv)
      step <- data.frame(parameter = ph, term = cand[k],
                         LR = tests[[k]][["LR"]], df = tests[[k]][["df"]],
                         p_value = pv[k], removed = pv[k] > alpha)
      path[[length(path) + 1]] <- step
      if (pv[k] <= alpha) break
      if (ph == "scale") sca <- setdiff(sca, cand[k]) else
        loc <- setdiff(loc, cand[k])
    }
  }
  list(location = loc, scale = sca, fit = safe_fit(data, loc, sca),
       path = do.call(rbind, path))
}

removable_terms <- function(tl) {
  parts <- strsplit(tl, ":", fixed = TRUE)
  keep <- vapply(seq_along(tl), function(i) {
    !any(vapply(seq_along(tl), function(j) {
      j != i && length(parts[[i]]) < length(parts[[j]]) &&
        all(parts[[i]] %in% parts[[j]])
    }, logical(1)))
  }, logical(1))
  tl[keep]
}

as_formula <- function(tl) {
  if (length(tl) == 0) ~ 1 else
    stats::as.formula(paste("~", paste(tl, collapse = " + ")))
}

safe_fit <- function(data, loc, sca) {
  tryCatch(suppressWarnings(fit_gev(data, as_formula(loc), as_formula(sca))),
           error = function(e) NULL)
}
