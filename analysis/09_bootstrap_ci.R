# Step 9 - year-block bootstrap confidence intervals (90 %, percentile)
#   theta (global), beta(0.5) of the season-only model, the season contrast
#   A_JJ(0.5) - A_ON(0.5), the period coefficient beta_4(omega) of the
#   season x period model, season-specific A(omega) bands, and the 6-d
#   extremal coefficients eta (incl. eta_JJ - eta_ON).
#
# Design: resample the 88 calendar years with replacement (blocks keep the
# six stations and the eight months of a year together); every replicate
# recomputes the pseudo-observations, theta, and all fits, so the intervals
# include the uncertainty of theta and of the ranks. Fits inside replicates
# use the same estimator with fewer optimiser starts (n_grid_starts = 1);
# the effect of that shortcut on the point estimate is reported below.
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("09 bootstrap (B = ", BOOT_B, ", cores = ", N_CORES, ")")

r <- load_pobs()
pairs <- pair_index()
W_BOOT <- if (QUICK) c(0.3, 0.5, 0.7) else seq(0.1, 0.9, by = 0.1)
LEVEL <- 0.90

design_s <- function(d) cbind(1, JJ = d$season == "JJ", AS = d$season == "AS",
                              ON = d$season == "ON") * 1
design_sp <- function(d) cbind(design_s(d), post = (d$period == "post") * 1)

count_warnings <- function(expr) {
  n <- 0L
  val <- withCallingHandlers(expr, warning = function(w) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  attr(val, "n_warnings") <- n
  val
}

# One evaluation of every bootstrapped quantity on a (possibly resampled)
# month x station table. Shared by the point estimate and the replicates.
statistic <- function(dat, n_grid_starts = 1) {
  dat[U_COLS] <- rank_columns(dat[U_COLS])
  theta <- estimate_theta_pairs(as.matrix(dat[U_COLS]), se = FALSE)$global
  per_pair <- lapply(seq_len(nrow(pairs)), function(k) {
    d <- dat[c(U_COLS[c(pairs$i[k], pairs$j[k])], "season", "period")]
    d <- d[stats::complete.cases(d), ]
    Us <- pseudo_obs_by_group(d[1:2], d$season)
    Usp <- pseudo_obs_by_group(d[1:2], interaction(d$season, d$period, drop = TRUE))
    fs <- count_warnings(pickands_cov_mle(Us, design_s(d), theta, W_BOOT,
                                          n_grid_starts = n_grid_starts))
    fsp <- count_warnings(pickands_cov_mle(Usp, design_sp(d), theta, W_BOOT,
                                           n_grid_starts = n_grid_starts))
    list(beta_s = fs$beta, beta4 = fsp$beta[, "post"],
         n_failed = sum(fs$diagnostics$status == "failed") +
           sum(fsp$diagnostics$status == "failed"))
  })
  cc <- dat[stats::complete.cases(dat[U_COLS]), ]
  eta <- vapply(SEASONS, function(s)
    extremal_coefficient(pseudo_obs(cc[cc$season == s, U_COLS]), theta), numeric(1))
  list(theta = theta, eta = eta,
       beta_s = simplify2array(lapply(per_pair, `[[`, "beta_s")),   # w x 4 x pair
       beta4 = sapply(per_pair, `[[`, "beta4"),                      # w x pair
       n_failed = sum(vapply(per_pair, `[[`, numeric(1), "n_failed")))
}

key <- list(digest_key(r), BOOT_B, SEED, W_BOOT)
point <- cached("bootstrap_point", key, list(
  full = statistic(r, n_grid_starts = 5),
  fast = statistic(r, n_grid_starts = 1)))
boot <- cached("bootstrap_replicates", key,
               boot_years(r, statistic, B = BOOT_B, seed = SEED, cores = N_CORES))
reps <- Filter(Negate(is.null), boot$replicates)
log_step(sprintf("  %d/%d replicates ok; %d failed grid fits inside replicates; %.0f s",
                 length(reps), boot$B, sum(vapply(reps, `[[`, numeric(1), "n_failed")),
                 attr(boot, "compute_seconds")))
shortcut_gap <- max(abs(point$full$beta_s - point$fast$beta_s), na.rm = TRUE)
log_step(sprintf("  max |beta(full starts) - beta(1 start)| on the original data: %.2g",
                 shortcut_gap))

ci_row <- function(est, draws) {
  ci <- boot_ci(draws, LEVEL)
  data.frame(estimate = est, ci_lower = ci[["lower"]], ci_upper = ci[["upper"]],
             boot_se = stats::sd(draws, na.rm = TRUE), n_boot = ci[["n_ok"]])
}
P <- point$full
k05 <- which.min(abs(W_BOOT - 0.5))

# ---- theta and eta ---------------------------------------------------------------
th_draw <- vapply(reps, `[[`, numeric(1), "theta")
eta_draw <- t(vapply(reps, `[[`, numeric(4), "eta"))
tab_eta <- rbind(
  cbind(quantity = "theta (global)", ci_row(P$theta, th_draw)),
  do.call(rbind, lapply(SEASONS, function(s)
    cbind(quantity = paste("eta", s), ci_row(P$eta[[s]], eta_draw[, s])))),
  cbind(quantity = "eta JJ - eta ON", ci_row(P$eta[["JJ"]] - P$eta[["ON"]],
                                             eta_draw[, "JJ"] - eta_draw[, "ON"])),
  cbind(quantity = "eta JJ - eta AM", ci_row(P$eta[["JJ"]] - P$eta[["AM"]],
                                             eta_draw[, "JJ"] - eta_draw[, "AM"])))
save_tab(tab_eta, "tab_20_bootstrap_theta_eta")

# ---- season-only beta(0.5) and the JJ - ON contrast ---------------------------
bs_draw <- simplify2array(lapply(reps, `[[`, "beta_s"))   # w x 4 x pair x B
coef_nm <- c("beta0 (AM level)", "beta1 (JJ - AM)", "beta2 (AS - AM)", "beta3 (ON - AM)")
tab_beta <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(p) {
  rows <- lapply(1:4, function(c) cbind(pair = pairs$pair[p], coefficient = coef_nm[c],
                                        ci_row(P$beta_s[k05, c, p], bs_draw[k05, c, p, ])))
  contrast <- cbind(pair = pairs$pair[p], coefficient = "A_JJ - A_ON (= beta1 - beta3)",
                    ci_row(P$beta_s[k05, 2, p] - P$beta_s[k05, 4, p],
                           bs_draw[k05, 2, p, ] - bs_draw[k05, 4, p, ]))
  do.call(rbind, c(rows, list(contrast)))
}))
tab_beta$excludes_zero <- tab_beta$ci_lower > 0 | tab_beta$ci_upper < 0
save_tab(tab_beta, "tab_21_bootstrap_season_beta_w05")

# ---- period coefficient beta_4(omega) ------------------------------------------
b4_draw <- simplify2array(lapply(reps, `[[`, "beta4"))    # w x pair x B
tab_b4 <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(p)
  do.call(rbind, lapply(seq_along(W_BOOT), function(k)
    cbind(pair = pairs$pair[p], w = W_BOOT[k],
          ci_row(P$beta4[k, p], b4_draw[k, p, ]))))))
tab_b4$excludes_zero <- tab_b4$ci_lower > 0 | tab_b4$ci_upper < 0
save_tab(tab_b4, "tab_22_bootstrap_beta4")

summ <- data.frame(
  B = boot$B, replicates_ok = length(reps), level = LEVEL,
  grid = paste(W_BOOT, collapse = " "),
  compute_seconds = round(attr(boot, "compute_seconds")),
  max_gap_full_vs_fast_starts = shortcut_gap,
  pairs_JJ_minus_ON_excludes_0 = sum(tab_beta$excludes_zero &
                                       grepl("A_JJ - A_ON", tab_beta$coefficient)),
  beta4_cells_excluding_0 = sum(tab_b4$excludes_zero),
  beta4_cells = nrow(tab_b4))
save_tab(summ, "tab_23_bootstrap_summary")
log_step(sprintf("  JJ-ON contrast excludes 0 for %d/15 pairs; beta4 CI excludes 0 in %d/%d cells",
                 summ$pairs_JJ_minus_ON_excludes_0, summ$beta4_cells_excluding_0,
                 summ$beta4_cells))

# ---- figures ------------------------------------------------------------------------
X_season <- rbind(AM = c(1, 0, 0, 0), JJ = c(1, 1, 0, 0),   # season profiles
                  AS = c(1, 0, 1, 0), ON = c(1, 0, 0, 1))
bands <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(p) {
  do.call(rbind, lapply(1:4, function(s) {
    est <- P$beta_s[, , p] %*% X_season[s, ]
    draws <- apply(bs_draw[, , p, , drop = FALSE], c(1, 4), function(b) sum(b * X_season[s, ]))
    ci <- t(apply(draws, 1, boot_ci, level = LEVEL))
    data.frame(pair = pairs$pair[p], season = factor(SEASONS[s], SEASONS),
               w = W_BOOT, A = as.numeric(est), lower = ci[, "lower"], upper = ci[, "upper"])
  }))
}))
save_tab(bands, "tab_24_bootstrap_A_bands")
p_bands <- ggplot(transform(bands, pair = factor(pair, pairs$pair)),
                  aes(w, A, colour = season, fill = season)) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ pair, ncol = 5) +
  scale_colour_manual(values = PAL_SEASON, aesthetics = c("colour", "fill")) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75)) +
  labs(title = sprintf("Season-specific Pickands functions with %d%% year-block bootstrap bands",
                       round(100 * LEVEL)),
       subtitle = sprintf("B = %d replicates; theta, ranks and fits re-estimated in every replicate",
                          length(reps)),
       x = expression(omega), y = expression(hat(A)(omega)), colour = NULL, fill = NULL)
save_fig(p_bands, "fig_14_bootstrap_A_bands", 11, 7)

p_b4 <- ggplot(transform(tab_b4, pair = factor(pair, pairs$pair)), aes(w, estimate)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#B2182B", alpha = 0.18) +
  geom_line(colour = "#B2182B") + geom_point(colour = "#B2182B", size = 0.9) +
  facet_wrap(~ pair, ncol = 5) +
  labs(title = expression("Period coefficient " * hat(beta)[4](omega) *
                            " (post-1975 shift in A) with 90% bootstrap bands"),
       subtitle = sprintf("%d of %d (pair, omega) intervals exclude 0",
                          summ$beta4_cells_excluding_0, summ$beta4_cells),
       x = expression(omega), y = expression(hat(beta)[4](omega)))
save_fig(p_b4, "fig_15_bootstrap_beta4", 11, 6)
finish_step("09")
