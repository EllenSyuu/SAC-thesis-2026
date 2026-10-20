# Every failure must be visible: NA + diagnostics + warning, never a
# plausible-looking number.

test_that("covariate MLE reports failure instead of returning a start value", {
  set.seed(3)
  U <- pseudo_obs(matrix(runif(200), ncol = 2))
  U[1, ] <- 0                       # xi = Inf -> every start infeasible
  X <- cbind(1, rep(0:1, 50))
  for (lk in c("identity", "logit")) {
    expect_warning(fit <- pickands_cov_mle(U, X, 0.4, w = c(0.3, 0.5), link = lk),
                   "no feasible optimum at 2 of 2")
    expect_true(all(is.na(fit$beta)))
    expect_true(all(is.na(fit$A)))
    expect_equal(fit$diagnostics$status, c("failed", "failed"))
    expect_equal(fit$diagnostics$n_feasible_starts, c(0L, 0L))
  }
})

test_that("successful covariate fits carry clean diagnostics", {
  set.seed(4)
  U <- pseudo_obs(rclayton_archimax(300, 0.4, 2))
  X <- cbind(1, rep(0:1, 150))
  fit <- expect_silent(pickands_cov_mle(U, X, 0.4, w = c(0.25, 0.5, 0.75)))
  expect_true(all(fit$diagnostics$status %in% c("ok", "maxit")))
  expect_true(all(fit$diagnostics$n_feasible_starts >= 1))
  expect_true(all(is.finite(fit$diagnostics$nll)))
})

test_that("a single covariate group (intercept only) is handled", {
  set.seed(5)
  U <- pseudo_obs(rclayton_archimax(200, 0.4, 2))
  fit <- pickands_cov_mle(U, matrix(1, 200, 1), 0.4, w = 0.5)
  expect_equal(dim(fit$A), c(1L, 1L))
  # identical to the homogeneous MLE
  expect_equal(unname(fit$A[1, 1]), pickands_at(U, c(0.5, 0.5), 0.4), tolerance = 1e-4)
})

test_that("pickands_at returns NA with a warning on degenerate input", {
  U <- cbind(c(0, 0.5, 0.7), c(0, 0.4, 0.9))
  expect_warning(a <- pickands_at(U, c(0.5, 0.5), 0.4), "non-finite")
  expect_true(is.na(a))
})

test_that("estimate_theta refuses raw (non-pseudo) data", {
  expect_error(estimate_theta(c(10, 20, 30, 40), c(1, 2, 3, 4)), "pseudo-observations")
})

gev_toy <- function(n = 240, seed = 11) {
  set.seed(seed)
  season <- factor(rep(c("AM", "JJ", "AS", "ON"), length.out = n),
                   levels = c("AM", "JJ", "AS", "ON"))
  mu <- 17 + c(0, 4, 6, 2)[as.integer(season)]
  data.frame(season = season,
             period = factor(rep(c("pre", "post"), each = n / 2), c("pre", "post")),
             max_precip = extRemes::revd(n, loc = mu, scale = 7, shape = 0.05))
}

test_that("safe_fit warns and returns a classed error object", {
  d <- gev_toy()
  expect_warning(res <- safe_fit(d, "not_a_column", character(0)),
                 "GEV fit failed \\[location ~ not_a_column \\| scale ~ 1\\]")
  expect_s3_class(res, "gev_fit_error")
  expect_match(res$error, "not_a_column")
})

test_that("gev_backward records failed fits in the selection path", {
  d <- gev_toy()
  warns <- capture_warnings(sel <- gev_backward(d, terms = c("season", "ghost")))
  expect_true(length(warns) >= 1)
  expect_true(all(grepl("GEV fit failed|every candidate", warns)))
  expect_true(any(grepl("fit failed", sel$path$note)))
})

test_that("gev_backward selects season on data generated with a season effect", {
  d <- gev_toy(n = 480)
  sel <- gev_backward(d, terms = c("period", "season"))
  expect_true("season" %in% sel$location)
  expect_s3_class(sel$fit, "fevd")
  expect_true(all(c("parameter", "term", "p_value", "removed", "note") %in%
                    names(sel$path)))
})

test_that("sandwich SE equals model SE when every observation is its own cluster", {
  d <- gev_toy(n = 400)
  fit <- fit_gev(d, ~ season, ~ 1)
  tab <- gev_coef_table(fit, cluster = seq_len(nrow(d)))
  # J ~ H for a correctly specified model, so the two should be close
  expect_lt(max(abs(log(tab$se_sandwich / tab$se))), 0.35)
  expect_equal(nrow(tab), length(fit$results$par))
})
