read_fixture <- function() {
  utils::read.csv(test_path("fixtures", "mini_daily.csv"), stringsAsFactors = FALSE)
}

test_that("a clean fixture passes and is summarised", {
  s <- validate_daily(read_fixture())
  expect_equal(s$n_rows, 182)
  expect_equal(s$n_stations, 2)
  expect_equal(s$n_missing_precip, 31 + 20)
  expect_equal(format(s$first_date), "2000-04-01")
})

test_that("each class of bad input fails with an actionable message", {
  d <- read_fixture()
  expect_error(validate_daily(d[0, ]), "zero rows")
  expect_error(validate_daily(d[, -6]), "missing required column\\(s\\): total_precip")
  bad <- d; bad$total_precip <- as.character(bad$total_precip)
  expect_error(validate_daily(bad), "'total_precip' must be numeric")
  bad <- d; bad$date[3] <- "2000-13-45"
  expect_error(validate_daily(bad), "cannot be parsed")
  bad <- d; bad$station_name[5] <- NA
  expect_error(validate_daily(bad), "missing station_name or date; first at row 5")
  bad <- d; bad$day[7] <- 30
  expect_error(validate_daily(bad), "disagree with date; first at row 7")
  bad <- d; bad$total_precip[9] <- -1
  expect_error(validate_daily(bad), "negative precipitation .* row 9")
  expect_error(validate_daily(rbind(d, d[10, ])), "duplicated \\(station_name, date\\)")
})

test_that("implausibly large values are flagged, not rejected", {
  d <- read_fixture()
  d$total_precip[1] <- 400
  expect_equal(validate_daily(d)$n_above_plausible, 1)
})

test_that("validate_monthly catches structural problems", {
  mm <- data.frame(station_name = "A", year = 2000, month = 4:5, max_precip = c(10, 12),
                   n_days = c(30, 31), season = "AM", period = "post", anomaly = 0.1)
  expect_equal(validate_monthly(mm)$n_station_months, 2)
  expect_error(validate_monthly(mm[0, ]), "no station-months")
  expect_error(validate_monthly(rbind(mm, mm[1, ])), "more than one row")
  bad <- mm; bad$anomaly[1] <- NA
  expect_error(validate_monthly(bad), "missing value\\(s\\) in 'anomaly'")
})

test_that("monthly_maxima computes maxima and valid-day counts per station-month", {
  d <- read_fixture()
  mm <- monthly_maxima(d)
  # station B, May is entirely missing -> no block maximum at all
  expect_false(any(mm$station_name == "STATION B" & mm$month == 5))
  expect_equal(nrow(mm), 5)
  a_jun <- mm[mm$station_name == "STATION A" & mm$month == 6, ]
  expect_equal(a_jun$n_days, 10)
  ref <- max(d$total_precip[d$station_name == "STATION A" & d$month == 6], na.rm = TRUE)
  expect_equal(a_jun$max_precip, ref)
  expect_true(all(mm$n_days <= 31))
})
