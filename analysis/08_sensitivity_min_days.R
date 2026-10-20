# Step 8 - sensitivity to incomplete months (MIN_DAYS)
# The thesis keeps every station-month with at least one valid day (default,
# MIN_DAYS = 0). A month observed on only a few days has a maximum that is
# biased low; dropping such months removes that bias but also discards data
# and can create gaps. Rather than choosing a threshold, this step refits the
# key quantities for a range of thresholds and reports how much they move.
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("08 MIN_DAYS sensitivity")

thresholds <- c(0, 5, 10, 15, 20, 25)
mm_all <- readr::read_csv(file.path(PROC, "monthly_maxima_all.csv"),
                          show_col_types = FALSE) |>
  mutate(season = factor(season, SEASONS), period = factor(period, c("pre", "post")),
         station = factor(station_name, STATIONS$key))
pairs <- pair_index()

one_threshold <- function(m) {
  mm <- filter(mm_all, n_days >= m)
  pooled <- fit_gev(as.data.frame(mm), ~ period * station + season, ~ season)
  ct <- gev_coef_table(pooled)
  pick <- function(lbl) ct$estimate[ct$covariate == lbl]
  # raw-rank pseudo-observations, as in the thesis
  r <- mm |>
    group_by(station) |>
    mutate(u = rank(max_precip) / (n() + 1)) |>
    ungroup() |>
    mutate(col = U_COLS[as.integer(station)]) |>
    select(date, season, col, u) |>
    pivot_wider(names_from = col, values_from = u)
  theta <- estimate_theta_pairs(as.matrix(r[U_COLS]))$global
  cc <- r[stats::complete.cases(r[U_COLS]), ]
  eta <- vapply(SEASONS, function(s) {
    extremal_coefficient(pseudo_obs(cc[cc$season == s, U_COLS]), theta)
  }, numeric(1))
  data.frame(min_days = m,
             station_months = nrow(mm),
             dropped = nrow(mm_all) - nrow(mm),
             complete_6d_months = nrow(cc),
             loc_AS_vs_AM = pick("location: seasonAS"),
             loc_post_drummondville = pick("location: periodpost"),
             shape_xi = pick("shape: intercept"),
             theta_global = theta,
             eta_AM = eta[["AM"]], eta_JJ = eta[["JJ"]],
             eta_AS = eta[["AS"]], eta_ON = eta[["ON"]])
}

sens <- cached("sensitivity_min_days",
               key = list(thresholds, digest_key(mm_all)),
               do.call(rbind, lapply(thresholds, one_threshold)))
save_tab(sens, "tab_18_sensitivity_min_days")

long <- sens |>
  select(-station_months, -complete_6d_months) |>
  pivot_longer(-min_days) |>
  mutate(name = factor(name, c("dropped", "loc_AS_vs_AM", "loc_post_drummondville",
                               "shape_xi", "theta_global", "eta_AM", "eta_JJ",
                               "eta_AS", "eta_ON"),
                       c("station-months dropped", "GEV location: AS vs AM (mm)",
                         "GEV location: post-1975 (mm)", "GEV shape xi",
                         "Clayton theta", "eta AM", "eta JJ", "eta AS", "eta ON")))
p <- ggplot(long, aes(min_days, value)) +
  geom_line(colour = "#0072B2") + geom_point(colour = "#0072B2") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey50") +
  facet_wrap(~ name, scales = "free_y", nrow = 3) +
  labs(title = "Sensitivity to the minimum number of valid days per month",
       subtitle = "MIN_DAYS = 0 (dotted) is the thesis default; no threshold is preferred here",
       x = "MIN_DAYS", y = NULL)
save_fig(p, "fig_13_sensitivity_min_days", 9, 6.5)
finish_step("08")
