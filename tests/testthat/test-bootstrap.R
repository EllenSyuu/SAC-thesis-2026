toy_years <- function() {
  data.frame(year = rep(2001:2010, each = 3), month = rep(1:3, 10),
             x = seq_len(30) + rnorm(30))
}

test_that("whole years are resampled as blocks", {
  set.seed(1)
  d <- toy_years()
  bt <- boot_years(d, function(rb) rb, B = 5, seed = 9)
  for (rb in bt$replicates) {
    expect_equal(nrow(rb), 30)
    expect_true(all(table(rb$year) %% 3 == 0))     # every drawn year kept whole
    expect_true(all(rb$year %in% d$year))
  }
})

test_that("results are identical for any number of cores", {
  skip_on_os("windows")
  d <- toy_years()
  stat <- function(rb) mean(rb$x) + stats::rnorm(1)   # uses the RNG inside
  a <- boot_years(d, stat, B = 8, seed = 42, cores = 1)
  b <- boot_years(d, stat, B = 8, seed = 42, cores = 2)
  expect_identical(unlist(a$replicates), unlist(b$replicates))
  expect_identical(a$seeds, b$seeds)
})

test_that("failing replicates are counted and reported, not dropped silently", {
  d <- toy_years()
  stat <- function(rb) if (2001 %in% rb$year) stop("boom") else 1
  expect_warning(bt <- boot_years(d, stat, B = 20, seed = 3), "replicates failed .*boom")
  expect_gt(bt$n_failed, 0)
  expect_equal(sum(vapply(bt$replicates, is.null, logical(1))), bt$n_failed)
  expect_equal(bt$errors, "boom")
})

test_that("boot_ci returns percentile limits and ignores non-finite draws", {
  ci <- boot_ci(c(1:99, NA, Inf), level = 0.90)
  expect_equal(unname(ci["n_ok"]), 99)
  expect_equal(unname(ci[c("lower", "upper")]),
               stats::quantile(1:99, c(0.05, 0.95), names = FALSE))
  expect_true(all(is.na(boot_ci(c(NA, 1))[1:2])))
})

test_that("rank_columns keeps NA and ranks within each column", {
  X <- cbind(a = c(3, 1, NA, 2), b = c(10, 40, 20, 30))
  U <- rank_columns(X)
  expect_true(is.na(U[3, "a"]))
  expect_equal(U[, "a"][-3], c(3, 1, 2) / 4)
  expect_equal(U[, "b"], c(1, 4, 2, 3) / 5)
})

test_that("the single-start shortcut used in the bootstrap matches the full search", {
  set.seed(8)
  d0 <- rclayton_archimax(300, 0.4, 1.5)
  d1 <- rclayton_archimax(300, 0.4, 4)
  U <- pseudo_obs(rbind(d0, d1))
  X <- cbind(1, rep(0:1, each = 300))
  full <- pickands_cov_mle(U, X, 0.4, w = c(0.3, 0.5, 0.7), n_grid_starts = 5)
  fast <- pickands_cov_mle(U, X, 0.4, w = c(0.3, 0.5, 0.7), n_grid_starts = 1)
  expect_equal(fast$beta, full$beta, tolerance = 1e-3)
})
