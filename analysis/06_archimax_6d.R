# Step 6 - six-dimensional extension (Chapter 5.3.2)
#   * extremal coefficient eta = 6 A(omega_0) at the simplex centroid, by season
#   * A along paths from bivariate faces to the centroid
# Pseudo-observations are re-ranked within each season so that every
# subgroup is a proper rank sample.
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("06 six-dimensional analysis")

POBS <- Sys.getenv("POBS", "raw")
sfx  <- if (POBS == "raw") "" else paste0("_", POBS)
r <- load_pobs(POBS)
THETA <- estimate_theta_pairs(as.matrix(r[U_COLS]))$global
cc <- r[stats::complete.cases(r[U_COLS]), ]
d <- length(U_COLS)
log_step(sprintf("  complete months: %d, theta = %.3f", nrow(cc), THETA))

U_all <- pseudo_obs(cc[U_COLS])
U_season <- lapply(setNames(SEASONS, SEASONS), function(s)
  pseudo_obs(cc[cc$season == s, U_COLS]))

eta <- rbind(
  data.frame(group = SEASONS, n = vapply(U_season, nrow, 1L),
             A_centroid = vapply(U_season, function(U)
               pickands_at(U, rep(1 / d, d), THETA), numeric(1))),
  data.frame(group = "All months (no covariates)", n = nrow(U_all),
             A_centroid = pickands_at(U_all, rep(1 / d, d), THETA))
)
eta$eta <- d * eta$A_centroid
save_tab(eta, paste0("tab_16_extremal_coefficient_6d", sfx))
log_step("  eta: ", paste(sprintf("%s=%.2f", eta$group, eta$eta), collapse = ", "))

# ---- simplex paths ------------------------------------------------------------------
t_grid <- seq(0, 1, by = if (QUICK) 0.1 else 0.02)
faces <- list(c(1, 2), c(1, 6), c(2, 4), c(4, 5), c(3, 6))
path_a <- do.call(rbind, lapply(faces, function(f) {
  data.frame(face = paste(STATIONS$short[f], collapse = " & "), t = t_grid,
             A = vapply(t_grid, function(t)
               pickands_at(U_all, simplex_path(f[1], f[2], t, d), THETA),
               numeric(1)))
}))
path_b <- do.call(rbind, lapply(SEASONS, function(s) {
  data.frame(season = factor(s, SEASONS), t = t_grid,
             A = vapply(t_grid, function(t)
               pickands_at(U_season[[s]], simplex_path(3, 6, t, d), THETA),
               numeric(1)))
}))
save_tab(path_a, paste0("tab_17_simplex_paths_pooled", sfx))
save_tab(path_b, paste0("tab_17_simplex_paths_season", sfx))

pa <- ggplot(path_a, aes(t, A, colour = face)) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  scale_colour_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0.6, 1.02)) +
  guides(colour = guide_legend(ncol = 2)) +
  labs(title = "(a) All months, five starting faces",
       x = "t  (0 = bivariate face, 1 = 6-d centroid)",
       y = expression(hat(A)(omega(t))), colour = NULL)
pb <- ggplot(path_b, aes(t, A, colour = season)) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = PAL_SEASON) +
  coord_cartesian(ylim = c(0.6, 1.02)) +
  labs(title = "(b) By season, from the Les Cèdres-Nicolet face",
       x = "t  (0 = bivariate face, 1 = 6-d centroid)", y = NULL, colour = NULL)
save_fig(pa + pb, paste0("fig_11_simplex_paths_6d", sfx), 11, 4.6)

pe <- ggplot(eta, aes(reorder(group, eta), eta,
                      fill = ifelse(group %in% SEASONS, group, "all"))) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", eta)), hjust = -0.15, size = 3.4) +
  geom_hline(yintercept = c(1, d), linetype = 2, colour = "grey50") +
  annotate("text", x = 5.45, y = c(1, d), label = c("complete dependence", "independence"),
           hjust = c(-0.05, 1.05), size = 3, colour = "grey40") +
  scale_fill_manual(values = c(PAL_SEASON, all = "grey60")) +
  coord_flip(ylim = c(0, d + 0.3), xlim = c(0.6, 5.6)) +
  labs(title = expression("Network-wide extremal coefficient " * eta == 6 * hat(A)(omega[0])),
       subtitle = "Six stations jointly: June-July extremes behave almost independently",
       x = NULL, y = expression(hat(eta)))
save_fig(pe, paste0("fig_12_extremal_coefficient", sfx), 7, 3.6)
log_step("06 done")
