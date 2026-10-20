# Step 2 - non-stationary GEV margins (Chapter 3.2) and pseudo-observations
#   * backward elimination per station, in both orders (robustness check)
#   * common station model  mu ~ period + season, sigma ~ season, xi ~ 1
#   * pooled model           mu ~ period * station + season, sigma ~ season
#   * Gumbel standardisation with the pooled fit -> rank pseudo-observations
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("02 marginal GEV")

mm <- readr::read_csv(file.path(PROC, "monthly_maxima.csv"),
                      show_col_types = FALSE) |>
  mutate(season = factor(season, levels = SEASONS),
         period = factor(period, levels = c("pre", "post")),
         station = factor(station_name, levels = STATIONS$key))
by_station <- split(as.data.frame(mm), mm$station)

# ---- 2a. backward elimination, both orders ------------------------------------
fmt_terms <- function(x) if (length(x)) paste(x, collapse = " + ") else "1"
sel <- lapply(STATIONS$key, function(k) {
  d <- by_station[[k]]
  a <- gev_backward(d, order = "scale_first")
  b <- gev_backward(d, order = "location_first")
  path <- rbind(cbind(order = "scale_first", a$path),
                cbind(order = "location_first", b$path))
  list(
    summary = data.frame(
      station = STATIONS$label[STATIONS$key == k],
      location_scale_first = fmt_terms(a$location),
      scale_scale_first    = fmt_terms(a$scale),
      BIC_scale_first      = gev_ic(a$fit)$BIC,
      location_loc_first   = fmt_terms(b$location),
      scale_loc_first      = fmt_terms(b$scale),
      BIC_loc_first        = gev_ic(b$fit)$BIC,
      same_model = setequal(a$location, b$location) && setequal(a$scale, b$scale)),
    path = cbind(station = STATIONS$label[STATIONS$key == k], path))
})
save_tab(do.call(rbind, lapply(sel, `[[`, "summary")), "tab_04_backward_elimination")
save_tab(do.call(rbind, lapply(sel, `[[`, "path")), "tab_04_backward_elimination_path")
log_step("  backward elimination done")

# ---- 2b. common station-specific model ------------------------------------------
st_fits <- lapply(by_station, fit_gev, location = ~ period + season,
                  scale = ~ season)
st_coef <- do.call(rbind, lapply(names(st_fits), function(k) {
  cbind(station = STATIONS$label[STATIONS$key == k],
        gev_coef_table(st_fits[[k]]), gev_ic(st_fits[[k]])[, c("loglik", "BIC")])
}))
save_tab(st_coef, "tab_05_gev_station_coefficients")

# ---- 2c. pooled model --------------------------------------------------------------
pooled <- fit_gev(as.data.frame(mm), location = ~ period * station + season,
                  scale = ~ season)
# The pooled likelihood treats the six stations as independent (composite
# likelihood). Godambe sandwich SEs, clustered by calendar month (spatial
# dependence) and by year (spatial + within-season serial dependence), show
# how optimistic the naive SEs are.
tab_pool <- gev_coef_table(pooled, cluster = mm$date)
tab_year <- gev_coef_table(pooled, cluster = mm$year)
tab_pool$se_sandwich_year <- tab_year$se_sandwich
tab_pool$se_ratio_month <- tab_pool$se_sandwich / tab_pool$se
save_tab(cbind(tab_pool, gev_ic(pooled)[rep(1, nrow(tab_pool)), c("loglik", "BIC")]),
         "tab_06_gev_pooled_coefficients")
log_step(sprintf("  pooled GEV: sandwich/naive SE ratio median %.2f (range %.2f-%.2f)",
                 median(tab_pool$se_ratio_month), min(tab_pool$se_ratio_month),
                 max(tab_pool$se_ratio_month)))

pars <- gev_params(pooled)
mm <- mm |>
  mutate(mu = pars$mu, sigma = pars$sigma, xi = pars$xi,
         z = gev_to_gumbel(max_precip, mu, sigma, xi))

# QQ plot on the common Gumbel scale (Fig. 3.12)
qq <- mm |> arrange(z) |>
  mutate(theoretical = -log(-log(ppoints(n()))))
p_qq <- ggplot(qq, aes(theoretical, z)) +
  geom_abline(colour = "grey50", linetype = 2) +
  geom_point(size = 0.8, alpha = 0.5, colour = "#0072B2") +
  labs(title = "Pooled non-stationary GEV: Gumbel QQ plot",
       subtitle = "mu ~ period x station + season, sigma ~ season, xi constant",
       x = "Theoretical Gumbel quantile", y = "Standardised observation")
save_fig(p_qq, "fig_04_pooled_gev_qq", 5.5, 5)

# Fitted GEV mean by season / period: stations vs pooled (Fig. 3.13)
grid <- expand.grid(season = factor(SEASONS, levels = SEASONS),
                    period = factor(c("pre", "post"), levels = c("pre", "post")))
gev_grid_mean <- function(fit, newdata) {
  p <- fit$results$par
  Xl <- model.matrix(fit$par.models$location, newdata)
  Xs <- model.matrix(fit$par.models$scale, newdata)
  mu <- Xl %*% p[grep("^mu", names(p))]
  sg <- Xs %*% p[grep("^sigma", names(p))]
  gev_mean(as.numeric(mu), as.numeric(sg), p[["shape"]])
}
means <- do.call(rbind, lapply(STATIONS$key, function(k) {
  nd <- transform(grid, station = factor(k, levels = STATIONS$key))
  rbind(cbind(nd, model = "Station-specific", mean = gev_grid_mean(st_fits[[k]], nd)),
        cbind(nd, model = "Pooled", mean = gev_grid_mean(pooled, nd)))
})) |>
  mutate(station = factor(STATIONS$label[match(station, STATIONS$key)],
                          levels = STATIONS$label))
save_tab(means, "tab_06_gev_fitted_means")
p_mean <- ggplot(means, aes(season, mean, colour = period, linetype = model,
                            group = interaction(period, model))) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  facet_wrap(~ station, nrow = 2) +
  scale_colour_manual(values = PAL_PERIOD,
                      labels = c(pre = "1930-1975", post = "1976-2017")) +
  labs(title = "Fitted GEV mean of monthly maxima",
       subtitle = "Season drives the level everywhere; the post-1975 shift is station-specific",
       x = NULL, y = "mm", colour = NULL, linetype = NULL)
save_fig(p_mean, "fig_04_gev_mean_station_vs_pooled", 9, 5.2)

# ---- 2d. pseudo-observations ---------------------------------------------------
# Rank pseudo-observations within each station, then a wide month x station
# table with season / period indicators (Chapter 4.1). Two versions:
#   * "raw"    : ranks of the monthly maxima - the input behind every number
#                reported in Chapters 4-5 of the thesis (default)
#   * "gumbel" : ranks of the pooled-GEV Gumbel residuals z (removes season /
#                period effects from the margins) - used as a sensitivity check
make_pobs <- function(value) {
  mm |>
    mutate(v = value) |>
    group_by(station) |>
    mutate(u = rank(v, ties.method = "average") / (n() + 1)) |>
    ungroup() |>
    mutate(col = U_COLS[as.integer(station)]) |>
    select(date, year, month, season, period, col, u) |>
    pivot_wider(names_from = col, values_from = u) |>
    select(date, year, month, season, period, all_of(U_COLS)) |>
    arrange(date)
}
pobs <- make_pobs(mm$max_precip)
readr::write_csv(pobs, file.path(PROC, "pseudo_obs.csv"))
readr::write_csv(make_pobs(mm$z), file.path(PROC, "pseudo_obs_gumbel.csv"))
readr::write_csv(select(mm, -station), file.path(PROC, "monthly_maxima_gumbel.csv"))
log_step("  pseudo-observations: ", nrow(pobs), " months, ",
         sum(complete.cases(pobs[U_COLS])), " complete")
log_step("02 done")
