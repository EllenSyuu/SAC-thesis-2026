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

#' Coefficient table for an fevd fit, with optional cluster-robust SEs
#'
#' Wald table based on the inverse observed information. When `cluster` is
#' given, a Godambe (sandwich) standard error is added:
#' `V = H^{-1} J H^{-1}`, `J = sum_c s_c s_c'`, where `s_c` is the score
#' summed over the observations of cluster `c` (e.g. all stations in one
#' month). This corrects the composite-likelihood SEs of a pooled model for
#' within-cluster dependence.
#'
#' @param fit An `fevd` object.
#' @param level Confidence level of the Wald intervals.
#' @param cluster Optional vector (one entry per observation) defining
#'   independent clusters.
#' @export
gev_coef_table <- function(fit, level = 0.95, cluster = NULL) {
  est <- fit$results$par
  V <- extRemes::parcov.fevd(fit)
  se <- sqrt(diag(V))
  zc <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(term = names(est), covariate = gev_term_labels(fit),
                    estimate = unname(est), se = unname(se),
                    z = unname(est / se),
                    p_value = 2 * stats::pnorm(-abs(unname(est / se))),
                    ci_lower = unname(est - zc * se),
                    ci_upper = unname(est + zc * se))
  if (!is.null(cluster)) {
    se_s <- sqrt(diag(gev_sandwich_vcov(fit, cluster, V)))
    out$se_sandwich <- se_s
    out$p_sandwich <- 2 * stats::pnorm(-abs(out$estimate / se_s))
    out$ci_lower_sandwich <- out$estimate - zc * se_s
    out$ci_upper_sandwich <- out$estimate + zc * se_s
  }
  out
}

#' Human-readable labels for fevd parameters (e.g. mu8 -> "location: seasonAS")
#'
#' @param fit An `fevd` object.
#' @export
gev_term_labels <- function(fit) {
  d <- fit$cov.data
  lab <- function(f, nm) {
    cn <- colnames(stats::model.matrix(f, d))
    paste0(nm, ": ", sub("^\\(Intercept\\)$", "intercept", cn))
  }
  c(lab(fit$par.models$location, "location"), lab(fit$par.models$scale, "scale"),
    "shape: intercept")
}

#' Per-observation GEV log-likelihood of a fitted non-stationary model
#'
#' @param par Parameter vector in `fevd` order (location, scale, shape).
#' @param fit An `fevd` object with linear location / scale models and a
#'   constant shape.
#' @keywords internal
gev_loglik_obs <- function(par, fit) {
  d <- fit$cov.data
  Xl <- stats::model.matrix(fit$par.models$location, d)
  Xs <- stats::model.matrix(fit$par.models$scale, d)
  pl <- ncol(Xl)
  ps <- ncol(Xs)
  mu <- as.numeric(Xl %*% par[seq_len(pl)])
  sg <- as.numeric(Xs %*% par[pl + seq_len(ps)])
  xi <- par[pl + ps + 1]
  t <- 1 + xi * (fit$x - mu) / sg
  -log(sg) - (1 + 1 / xi) * log(t) - t^(-1 / xi)
}

#' Godambe sandwich covariance for an fevd fit
#'
#' @param fit An `fevd` object.
#' @param cluster Cluster labels, one per observation.
#' @param Hinv Inverse observed information (defaults to `parcov.fevd`).
#' @export
gev_sandwich_vcov <- function(fit, cluster, Hinv = extRemes::parcov.fevd(fit)) {
  par <- fit$results$par
  if (length(cluster) != length(fit$x)) {
    stop("cluster must have one entry per observation.", call. = FALSE)
  }
  ll0 <- gev_loglik_obs(par, fit)
  if (abs(sum(ll0) + fit$results$value) > 1e-4 * abs(fit$results$value)) {
    stop("Could not reproduce the fevd log-likelihood; sandwich SE unavailable ",
         "for this model specification.", call. = FALSE)
  }
  scores <- vapply(seq_along(par), function(k) {
    h <- 1e-5 * max(1, abs(par[k]))
    up <- par
    dn <- par
    up[k] <- up[k] + h
    dn[k] <- dn[k] - h
    (gev_loglik_obs(up, fit) - gev_loglik_obs(dn, fit)) / (2 * h)
  }, numeric(length(ll0)))
  S <- rowsum(scores, group = cluster)
  J <- crossprod(S)
  V <- Hinv %*% J %*% Hinv
  dimnames(V) <- list(names(par), names(par))
  V
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
#' Failed fits are never dropped silently: each failure raises a warning and
#' is recorded in `path` (column `note`). A candidate whose reduced model
#' cannot be fitted is not eligible for removal at that step.
#'
#' @param data Station data frame.
#' @param terms Character vector of candidate terms (shared start for both
#'   parameters).
#' @param order `"scale_first"` (thesis) or `"location_first"`.
#' @param alpha LRT level.
#' @return A list with final `location`, `scale` terms, the `fit`, and a
#'   `path` data frame recording each step (`parameter`, `term`, `LR`,
#'   `df`, `p_value`, `removed`, `note`).
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
  record <- function(ph, term, LR = NA_real_, df = NA_real_, p = NA_real_,
                     removed = FALSE, note = "") {
    path[[length(path) + 1]] <<- data.frame(
      parameter = ph, term = term, LR = LR, df = df, p_value = p,
      removed = removed, note = note, stringsAsFactors = FALSE)
  }
  phases <- if (order == "scale_first") c("scale", "location") else
    c("location", "scale")
  for (ph in phases) {
    repeat {
      cur <- safe_fit(data, loc, sca)
      if (fit_failed(cur)) {
        record(ph, "<current model>", note = paste("fit failed:", cur$error))
        break
      }
      target <- if (ph == "scale") sca else loc
      cand <- removable_terms(target)
      if (length(cand) == 0) break
      tests <- lapply(cand, function(tm) {
        new <- setdiff(target, tm)
        red <- if (ph == "scale") safe_fit(data, loc, new) else
          safe_fit(data, new, sca)
        if (fit_failed(red)) {
          record(ph, tm, note = paste("reduced fit failed:", red$error))
          return(c(LR = NA, df = NA, p = NA))
        }
        lr <- 2 * (red$results$value - cur$results$value)
        df <- length(cur$results$par) - length(red$results$par)
        c(LR = lr, df = df, p = stats::pchisq(lr, df, lower.tail = FALSE))
      })
      pv <- vapply(tests, function(t) t[["p"]], numeric(1))
      if (all(is.na(pv))) {
        warning("gev_backward: every candidate reduction failed in the ", ph,
                " phase; stopping this phase.", call. = FALSE)
        break
      }
      k <- which.max(pv)
      record(ph, cand[k], tests[[k]][["LR"]], tests[[k]][["df"]], pv[k],
             removed = pv[k] > alpha)
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

#' Fit a GEV and surface (rather than swallow) failures
#'
#' Returns the `fevd` fit or, on error, an object of class `gev_fit_error`
#' (a list with element `error`) after emitting a warning that names both
#' formulas and the error message.
#'
#' @param data Data frame.
#' @param loc,sca Character vectors of location / scale terms.
#' @keywords internal
safe_fit <- function(data, loc, sca) {
  tryCatch(suppressWarnings(fit_gev(data, as_formula(loc), as_formula(sca))),
           error = function(e) {
             msg <- conditionMessage(e)
             f <- function(tl) if (length(tl)) paste(tl, collapse = " + ") else "1"
             warning(sprintf("GEV fit failed [location ~ %s | scale ~ %s]: %s",
                             f(loc), f(sca), msg), call. = FALSE)
             structure(list(error = msg), class = "gev_fit_error")
           })
}

fit_failed <- function(x) inherits(x, "gev_fit_error")
