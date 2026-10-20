# -----------------------------------------------------------------------------
# Input validation. Every check fails loudly with an actionable message; on
# success a one-row summary is returned so the pipeline can log it.
# -----------------------------------------------------------------------------

#' Validate daily precipitation records
#'
#' Checks, in order: non-empty input; required columns; column types;
#' parseable dates consistent with `year` / `month` / `day`; no missing
#' station or date; no negative precipitation; no duplicated
#' `(station_name, date)`. Missing `total_precip` values are allowed (they
#' are common in ECCC records) and are counted in the summary.
#'
#' @param daily Data frame of daily records.
#' @param max_plausible Values above this (mm/day) are reported in the
#'   summary as `n_above_plausible`, not rejected.
#' @return Invisibly, a one-row data frame summarising the input.
#' @export
validate_daily <- function(daily, max_plausible = 250) {
  fail <- function(...) stop("validate_daily(): ", ..., call. = FALSE)
  if (!is.data.frame(daily)) fail("input must be a data frame.")
  if (nrow(daily) == 0) fail("input has zero rows.")

  req <- c("station_name", "date", "year", "month", "day", "total_precip")
  miss <- setdiff(req, names(daily))
  if (length(miss)) {
    fail("missing required column(s): ", paste(miss, collapse = ", "),
         ". Found: ", paste(names(daily), collapse = ", "), ".")
  }
  for (v in c("year", "month", "day", "total_precip")) {
    if (!is.numeric(daily[[v]])) {
      fail("column '", v, "' must be numeric, not ", class(daily[[v]])[1],
           " (coerce with as.numeric() after checking for non-numeric codes).")
    }
  }

  d <- tryCatch(as.Date(daily$date), error = function(e) NULL)
  if (is.null(d) || anyNA(d[!is.na(daily$date)])) {
    fail("column 'date' contains values that cannot be parsed as dates ",
         "(expected YYYY-MM-DD).")
  }
  key_na <- is.na(daily$station_name) | is.na(d)
  if (any(key_na)) {
    fail(sum(key_na), " row(s) have a missing station_name or date; first at row ",
         which(key_na)[1], ".")
  }
  ymd_bad <- as.integer(format(d, "%Y")) != daily$year |
    as.integer(format(d, "%m")) != daily$month |
    as.integer(format(d, "%d")) != daily$day
  if (any(ymd_bad, na.rm = TRUE)) {
    i <- which(ymd_bad)[1]
    fail(sum(ymd_bad, na.rm = TRUE), " row(s) where year/month/day disagree with ",
         "date; first at row ", i, " (", format(d[i]), ").")
  }
  neg <- which(daily$total_precip < 0)
  if (length(neg)) {
    fail(length(neg), " negative precipitation value(s); first at row ", neg[1],
         " (", daily$station_name[neg[1]], ", ", format(d[neg[1]]), ").")
  }
  dup <- duplicated(data.frame(daily$station_name, d))
  if (any(dup)) {
    i <- which(dup)[1]
    fail(sum(dup), " duplicated (station_name, date) row(s); first: ",
         daily$station_name[i], " ", format(d[i]),
         ". De-duplicate the download before running the pipeline.")
  }

  invisible(data.frame(
    n_rows = nrow(daily),
    n_stations = length(unique(daily$station_name)),
    first_date = min(d), last_date = max(d),
    n_missing_precip = sum(is.na(daily$total_precip)),
    share_missing_precip = mean(is.na(daily$total_precip)),
    n_above_plausible = sum(daily$total_precip > max_plausible, na.rm = TRUE),
    max_precip = max(daily$total_precip, na.rm = TRUE)
  ))
}

#' Validate the monthly-maxima table used for modelling
#'
#' @param mm Data frame with `station_name`, `year`, `month`, `max_precip`,
#'   `n_days`, `season`, `period`, `anomaly`.
#' @return Invisibly, a one-row summary.
#' @export
validate_monthly <- function(mm) {
  fail <- function(...) stop("validate_monthly(): ", ..., call. = FALSE)
  if (nrow(mm) == 0) fail("no station-months left (check YEARS / MIN_DAYS).")
  req <- c("station_name", "year", "month", "max_precip", "n_days",
           "season", "period", "anomaly")
  miss <- setdiff(req, names(mm))
  if (length(miss)) fail("missing column(s): ", paste(miss, collapse = ", "), ".")
  if (any(duplicated(mm[c("station_name", "year", "month")]))) {
    fail("more than one row per station-month.")
  }
  for (v in c("max_precip", "season", "period", "anomaly")) {
    if (anyNA(mm[[v]])) fail(sum(is.na(mm[[v]])), " missing value(s) in '", v, "'.")
  }
  if (any(mm$n_days < 1 | mm$n_days > 31)) fail("n_days outside 1..31.")
  invisible(data.frame(n_station_months = nrow(mm),
                       n_stations = length(unique(mm$station_name)),
                       min_n_days = min(mm$n_days),
                       n_incomplete_lt20 = sum(mm$n_days < 20)))
}
