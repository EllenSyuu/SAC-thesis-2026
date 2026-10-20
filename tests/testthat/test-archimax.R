test_that("generator and inverse are inverses", {
  u <- c(0.01, 0.3, 0.7, 0.99)
  expect_equal(clayton_psi(clayton_phi(u, 0.4), 0.4), u)
})

test_that("closed-form moments of Z match numerical integration", {
  th <- 0.35
  ez <- integrate(function(x) clayton_psi(x, th), 0, Inf)$value
  expect_equal(clayton_EZ(th), ez, tolerance = 1e-6)
  elog <- integrate(function(y) y * exp(y) * (1 + th * exp(y))^(-1 / th - 1),
                    -40, 60)$value
  expect_equal(clayton_ElogZ(th), elog, tolerance = 1e-6)
})

test_that("xi_transform drops zero-weight coordinates", {
  U <- matrix(c(0.2, 0.5, 0.9, 0.4, 0.6, 0.1), ncol = 3)
  expect_equal(xi_transform(U, c(0.5, 0.5, 0), 0.4),
               xi_transform(U[, 1:2], c(0.5, 0.5), 0.4))
  expect_error(xi_transform(U, c(0.5, 0.6, 0), 0.4), "simplex")
})

test_that("MLE estimate always respects the Pickands bounds", {
  set.seed(1)
  U <- pseudo_obs(matrix(runif(400), ncol = 2))
  w <- seq(0.05, 0.95, by = 0.05)
  A <- pickands_curve(U, 0.4, w, method = "mle")
  expect_true(all(A >= pmax(w, 1 - w) - 1e-8 & A <= 1 + 1e-8))
  expect_equal(pickands_at(U, c(1, 0), 0.4), 1)
})

test_that("theta and A are recovered from a simulated Clayton-Archimax sample", {
  skip_on_cran()
  set.seed(2026)
  dat <- rclayton_archimax(800, theta = 0.4, rho = 3)
  U <- pseudo_obs(dat)
  th <- estimate_theta(U[, 1], U[, 2])
  expect_lt(abs(th$theta - 0.4), 3 * th$se)
  w <- seq(0.05, 0.95, by = 0.05)
  A_hat <- pickands_curve(U, th$theta, w, method = "mle")
  expect_lt(ise(w, A_hat, pickands_logistic(w, 3)), 0.005)
})

test_that("covariate MLE separates two dependence regimes", {
  skip_on_cran()
  set.seed(7)
  d0 <- rclayton_archimax(500, 0.4, 1.5)
  d1 <- rclayton_archimax(500, 0.4, 4.5)
  z <- rep(0:1, each = 500)
  U <- pseudo_obs(rbind(d0, d1))
  X <- cbind(intercept = 1, z = z)
  for (lk in c("identity", "logit")) {
    fit <- pickands_cov_mle(U, X, theta = 0.4, w = 0.5, link = lk)
    expect_gt(fit$A[1, 1], fit$A[1, 2])  # z = 1 has stronger dependence
    expect_true(all(fit$A >= 0.5 & fit$A <= 1))
  }
})

test_that("tail-dependence formulas match closed forms", {
  td <- tail_dependence_archimax(0.75, 0.35)
  expect_equal(td$lambda_U, 0.5)
  expect_equal(td$lambda_L, 1.5^(-1 / 0.35))
})

test_that("Gumbel transform inverts the GEV cdf", {
  x <- c(10, 25, 60)
  z <- gev_to_gumbel(x, mu = 17, sigma = 7, xi = 0.06)
  F <- extRemes::pevd(x, loc = 17, scale = 7, shape = 0.06)
  expect_equal(exp(-exp(-z)), F)
})

test_that("season / period helpers", {
  expect_equal(as.character(qc_season(c(3, 4, 7, 9, 11, 12))),
               c(NA, "AM", "JJ", "AS", "ON", NA))
  expect_equal(as.character(qc_period(c(1975, 1976))), c("pre", "post"))
})

test_that("marginality is respected in backward elimination", {
  expect_setequal(removable_terms(c("a", "b", "a:b", "c")), c("a:b", "c"))
})
