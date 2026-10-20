# Step 1 - daily ECCC records -> monthly block maxima with covariates
#   * corrects one verified transcription error (St-Charles, 1996-07-02)
#   * keeps Apr-Nov (rain season) for 1930-2017
#   * joins NOAA North-America land temperature anomalies
#   * exploratory tables / figures for Chapter 3.1
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("01 prepare data")

daily <- readr::read_csv(file.path(RAW, "eccc_daily_precip.csv.gz"),
                         show_col_types = FALSE,
                         col_types = readr::cols(climate_id = "c")) |>
  mutate(year = as.integer(year), month = as.integer(month),
         day = as.integer(day))

# 147.0 mm on 1996-07-02 at St-Charles is a decimal-shift error: no nearby
# station (St-Alexis-des-Monts: 16.2 mm) or news report supports it.
fix <- with(daily, station_name == "ST CHARLES DE MANDEVILLE" & year == 1996 &
              month == 7 & day == 2)
stopifnot(sum(fix) == 1, daily$total_precip[fix] == 147)
daily$total_precip[fix] <- 14.7

# ---- missingness by calendar month (Fig. 3.3) -------------------------------
miss <- daily |>
  filter(year >= YEARS[1], year <= 2019) |>
  group_by(station_name, year, month) |>
  summarise(any_obs = any(!is.na(total_precip)), .groups = "drop") |>
  complete(station_name, year = YEARS[1]:2019, month = 1:12,
           fill = list(any_obs = FALSE)) |>
  group_by(station_name, month) |>
  summarise(missing_rate = mean(!any_obs), .groups = "drop") |>
  mutate(station = STATIONS$label[match(station_name, STATIONS$key)])

p_miss <- ggplot(miss, aes(factor(month, labels = month.abb), station,
                           fill = missing_rate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = scales::percent(missing_rate, 1)), size = 2.8) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C",
                      labels = scales::percent, name = "Missing months") +
  labs(title = "Share of station-months without any record, 1930-2019",
       x = NULL, y = NULL)
save_fig(p_miss, "fig_03_missingness", 9, 3.6)

# ---- monthly maxima + covariates --------------------------------------------
anom <- readr::read_csv(file.path(RAW, "noaa_na_land_temp_anomaly.csv"),
                        show_col_types = FALSE)

mm <- monthly_maxima(daily) |>
  filter(year >= YEARS[1], year <= YEARS[2]) |>
  mutate(season = qc_season(month), period = qc_period(year, PERIOD_CUT)) |>
  filter(!is.na(season)) |>
  left_join(anom, by = c("year", "month")) |>
  mutate(station = factor(station_name, levels = STATIONS$key),
         date = as.Date(sprintf("%d-%02d-01", year, month))) |>
  arrange(station, date)
stopifnot(!anyNA(mm$anomaly))
readr::write_csv(mm, file.path(PROC, "monthly_maxima.csv"))
log_step("  monthly maxima: ", nrow(mm), " station-months")

# ---- Table 3.3 / 3.4: counts, means and standard errors ----------------------
desc <- mm |>
  group_by(period, season) |>
  summarise(n = n(), mean_mm = mean(max_precip),
            se_mm = sd(max_precip) / sqrt(n()), .groups = "drop")
save_tab(desc, "tab_03_seasonal_summary")

welch <- mm |>
  group_by(period) |>
  summarise(F = oneway.test(max_precip ~ season)$statistic,
            p_value = oneway.test(max_precip ~ season)$p.value)
shift <- mm |>
  group_by(season) |>
  summarise(diff_post_minus_pre = mean(max_precip[period == "post"]) -
              mean(max_precip[period == "pre"]),
            p_value = t.test(max_precip ~ period)$p.value)
save_tab(welch, "tab_03_welch_anova_by_period")
save_tab(shift, "tab_03_period_shift_by_season")

p_means <- ggplot(desc, aes(season, mean_mm, colour = period, group = period)) +
  geom_line(linewidth = 0.8) +
  geom_pointrange(aes(ymin = mean_mm - 1.96 * se_mm,
                      ymax = mean_mm + 1.96 * se_mm),
                  position = position_dodge(0.15)) +
  scale_colour_manual(values = PAL_PERIOD,
                      labels = c(pre = "1930-1975", post = "1976-2017")) +
  labs(title = "Mean monthly maximum precipitation (±95% CI)",
       x = NULL, y = "mm", colour = NULL)
save_fig(p_means, "fig_03_seasonal_means", 6, 4)

# ---- temperature anomaly motivating the 1975 split (Fig. 3.2) ---------------
anom_y <- anom |> group_by(year) |> summarise(anomaly = mean(anomaly))
p_anom <- ggplot(anom_y, aes(year, anomaly)) +
  geom_col(aes(fill = anomaly > 0), width = 0.9, show.legend = FALSE) +
  geom_smooth(method = "loess", span = 0.4, se = FALSE, colour = "black",
              linewidth = 0.7, formula = y ~ x) +
  geom_vline(xintercept = PERIOD_CUT + 0.5, linetype = 2) +
  annotate("text", x = PERIOD_CUT + 2, y = max(anom_y$anomaly),
           label = "1975 split", hjust = 0, size = 3.2) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC")) +
  labs(title = "North-America land temperature anomaly (Apr-Nov mean)",
       subtitle = "NOAA NCEI Climate at a Glance, relative to 1910-2000",
       x = NULL, y = "°C")
save_fig(p_anom, "fig_03_temperature_anomaly", 8, 3.6)

# ---- station map (Fig. 3.1) ---------------------------------------------------
if (requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  prov <- rnaturalearth::ne_states(country = c("canada", "united states of america"),
                                   returnclass = "sf")
  p_map <- ggplot() +
    geom_sf(data = prov, fill = "grey96", colour = "grey55", linewidth = 0.3) +
    geom_point(data = STATIONS, aes(lon, lat), size = 3, colour = "#0072B2") +
    geom_text(data = STATIONS, aes(lon, lat, label = label), size = 3.4,
              nudge_y = 0.09, fontface = "bold") +
    coord_sf(xlim = c(-74.8, -71.0), ylim = c(45.0, 46.8), expand = FALSE) +
    labs(title = "Six ECCC stations in southern Québec", x = NULL, y = NULL)
  save_fig(p_map, "fig_03_station_map", 7.5, 4.6)
}
readr::write_csv(STATIONS, file.path(PROC, "stations.csv"))
log_step("01 done")
