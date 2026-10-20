# Step 5 - Clayton-Archimax inference across the 15 station pairs (Ch. 5.3.1)
#   * moment estimator of theta per pair (90% Wald CI) and global theta
#   * season x period linear model for A(omega | season, period)
#       -> distribution of the period coefficient beta_4(omega)
#       -> model-implied lambda_L / lambda_U per season-period cell
#   * season-only model: beta(0.5) table and A(omega | season) curves
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("05 pairwise Clayton-Archimax")

POBS <- Sys.getenv("POBS", "raw")     # "raw" (thesis) or "gumbel" (sensitivity)
sfx  <- if (POBS == "raw") "" else paste0("_", POBS)
r <- load_pobs(POBS)
pairs <- pair_index()

# ---- theta ------------------------------------------------------------------------
U <- as.matrix(r[U_COLS])
colnames(U) <- STATIONS$label
th <- estimate_theta_pairs(U, alpha = 0.10)
THETA <- th$global
save_tab(rbind(th$pairs,
               data.frame(station_1 = "Global (mean of pairs)", station_2 = "",
                          n = NA, theta = THETA, se = NA, ci_lower = NA,
                          ci_upper = NA)), paste0("tab_11_theta_pairwise", sfx))
log_step(sprintf("  global theta = %.3f; %d/15 CIs contain it", THETA,
                 sum(th$pairs$ci_lower <= THETA & th$pairs$ci_upper >= THETA)))

p_th <- ggplot(transform(th$pairs, pair = factor(pairs$pair, rev(pairs$pair))),
               aes(theta, pair)) +
  geom_vline(xintercept = THETA, linetype = 2, colour = "#D55E00") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25,
                 colour = "grey40") +
  geom_point(size = 2.2, colour = "#0072B2") +
  annotate("text", x = THETA, y = 15.6, hjust = -0.1, size = 3.2,
           label = sprintf("global theta = %.2f", THETA), colour = "#D55E00") +
  labs(title = "Clayton generator parameter by station pair (90% CI)",
       subtitle = sprintf("%d of 15 intervals cover the global estimate; all overlap -> a common theta is tenable",
                          sum(th$pairs$ci_lower <= THETA & th$pairs$ci_upper >= THETA)),
       x = expression(hat(theta)), y = NULL)
save_fig(p_th, paste0("fig_07_theta_pairwise", sfx), 7, 5)

# Conditional copula given the covariates: pseudo-observations are re-ranked
# inside each covariate cell, otherwise the seasonal shift in the *margins*
# (e.g. all April-May maxima being low) masquerades as dependence.
pair_data <- function(i, j, by = c("season", "period")) {
  d <- r[c(U_COLS[c(i, j)], "season", "period")]
  d <- d[stats::complete.cases(d), ]
  g <- interaction(d[by], drop = TRUE)
  d[1:2] <- pseudo_obs_by_group(d[1:2], g)
  d
}

# ---- season x period model --------------------------------------------------------
cells <- expand.grid(season = factor(SEASONS, SEASONS),
                     period = factor(c("pre", "post"), c("pre", "post")))
design_sp <- function(d) cbind(intercept = 1, JJ = d$season == "JJ",
                               AS = d$season == "AS", ON = d$season == "ON",
                               post = d$period == "post") * 1
X_cells <- design_sp(cells)

fit_sp <- lapply(seq_len(nrow(pairs)), function(k) {
  d <- pair_data(pairs$i[k], pairs$j[k])
  pickands_cov_mle(as.matrix(d[1:2]), design_sp(d), THETA, W_GRID)
})
log_step("  season x period model fitted")

beta4 <- data.frame(pair = rep(pairs$pair, each = length(W_GRID)),
                    w = rep(W_GRID, nrow(pairs)),
                    beta4 = unlist(lapply(fit_sp, function(f) f$beta[, "post"])))
save_tab(beta4, paste0("tab_12_beta4_period_effect", sfx))
p_b4 <- ggplot(beta4, aes(beta4)) +
  geom_histogram(binwidth = 0.005, boundary = 0, fill = "#0072B2", colour = "white") +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = expression("Period effect " * hat(beta)[4](omega) *
                            " across 15 pairs and the full " * omega * "-grid"),
       subtitle = sprintf("Median %.3f; %.0f%% of estimates within +/-0.02 (season effects reach 0.1-0.15)",
                          median(beta4$beta4), 100 * mean(abs(beta4$beta4) < 0.02)),
       x = expression(hat(beta)[4](omega)), y = "Count")
save_fig(p_b4, paste0("fig_08_beta4_histogram", sfx), 7, 4)

k05 <- which.min(abs(W_GRID - 0.5))
lam_sp <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
  A_half <- as.numeric(X_cells %*% fit_sp[[k]]$beta[k05, ])
  cbind(pair = pairs$pair[k], i = pairs$i[k], j = pairs$j[k], cells,
        A_half = A_half, tail_dependence_archimax(A_half, THETA))
}))
save_tab(lam_sp, paste0("tab_13_lambda_season_period", sfx))

p_lam <- lam_sp |>
  pivot_longer(c(lambda_L, lambda_U), names_to = "coef") |>
  mutate(coef = factor(coef, c("lambda_U", "lambda_L"),
                       c("Upper tail  lambda_U", "Lower tail  lambda_L")),
         pair = factor(pair, rev(pairs$pair)),
         cell = factor(paste(season, period), paste(rep(SEASONS, each = 2),
                                                    c("pre", "post")))) |>
  ggplot(aes(cell, pair, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.3) +
  facet_wrap(~ coef) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = expression(hat(lambda))) +
  labs(title = "Model-implied tail dependence by season and period",
       subtitle = "Read across: JJ is near tail-independent and ON/AM strongest for most pairs; pre/post-1975 shifts are smaller and pair-specific",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
save_fig(p_lam, paste0("fig_09_lambda_season_period", sfx), 11, 6)

# ---- season-only model -------------------------------------------------------------
design_s <- function(d) cbind(intercept = 1, JJ = d$season == "JJ",
                              AS = d$season == "AS", ON = d$season == "ON") * 1
fit_s <- lapply(seq_len(nrow(pairs)), function(k) {
  d <- pair_data(pairs$i[k], pairs$j[k], by = "season")
  pickands_cov_mle(as.matrix(d[1:2]), design_s(d), THETA, W_GRID)
})
log_step("  season-only model fitted")

# Optimiser diagnostics for every pair x grid point (never silently dropped)
diag_all <- rbind(
  do.call(rbind, lapply(seq_len(nrow(pairs)), function(k)
    cbind(model = "season x period", pair = pairs$pair[k], fit_sp[[k]]$diagnostics))),
  do.call(rbind, lapply(seq_len(nrow(pairs)), function(k)
    cbind(model = "season only", pair = pairs$pair[k], fit_s[[k]]$diagnostics))))
save_tab(diag_all, paste0("tab_19_optimizer_diagnostics", sfx))
log_step(sprintf("  optimiser: %d fits, %d failed, %d hit maxit, %d with a profile on the A-bound",
                 nrow(diag_all), sum(diag_all$status == "failed"),
                 sum(diag_all$status == "maxit"), sum(diag_all$at_bound)))

beta05 <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
  b <- fit_s[[k]]$beta[k05, ]
  data.frame(station_1 = STATIONS$short[pairs$i[k]],
             station_2 = STATIONS$short[pairs$j[k]],
             beta0_AM = b[1], beta1_JJ = b[2], beta2_AS = b[3], beta3_ON = b[4])
}))
save_tab(beta05, paste0("tab_14_season_beta_w05", sfx), digits = 3)

X_season <- design_s(data.frame(season = factor(SEASONS, SEASONS)))
curves <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
  A <- fit_s[[k]]$beta %*% t(X_season)
  data.frame(pair = pairs$pair[k], w = rep(W_GRID, 4),
             season = factor(rep(SEASONS, each = length(W_GRID)), SEASONS),
             A = as.numeric(A))
}))
save_tab(curves, paste0("tab_15_season_pickands_curves", sfx))
p_A <- ggplot(transform(curves, pair = factor(pair, pairs$pair)),
              aes(w, A, colour = season)) +
  geom_line(aes(y = pmax(w, 1 - w)), colour = "grey70", linetype = 3) +
  geom_hline(yintercept = 1, colour = "grey70", linetype = 3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ pair, ncol = 5) +
  scale_colour_manual(values = PAL_SEASON) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75)) +
  labs(title = "Season-specific Pickands dependence functions, Clayton-Archimax MLE",
       subtitle = "Lower curve = stronger extremal dependence. JJ (orange) is closest to independence, ON (purple) the most coherent",
       x = expression(omega), y = expression(hat(A)(omega)), colour = NULL)
save_fig(p_A, paste0("fig_10_pickands_by_season", sfx), 11, 7)
log_step("05 done")
