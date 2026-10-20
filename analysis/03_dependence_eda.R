# Step 3 - exploratory analysis of the dependence structure (Chapter 4)
#   * Kendall's tau overall / by season / by period
#   * empirical upper & lower tail dependence across thresholds
#   * test of extreme-value dependence (Ben Ghorbal et al., 2009; evTestK)
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("03 dependence EDA")

r <- load_pobs()
pairs <- pair_index()
lab_pair <- function(i, j) paste(STATIONS$abbr[i], STATIONS$abbr[j], sep = "-")

# ---- Kendall's tau ---------------------------------------------------------------
tau_long <- function(dat, group) {
  U <- as.matrix(dat[U_COLS])
  do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
    i <- pairs$i[k]; j <- pairs$j[k]
    data.frame(group = group, i = i, j = j,
               tau = cor(U[, i], U[, j], method = "kendall",
                         use = "pairwise.complete.obs"))
  }))
}
tau <- rbind(
  tau_long(r, "All months"),
  do.call(rbind, lapply(SEASONS, function(s) tau_long(filter(r, season == s), s))),
  tau_long(filter(r, period == "pre"), "1930-1975"),
  tau_long(filter(r, period == "post"), "1976-2017")
) |>
  mutate(group = factor(group, levels = c("All months", SEASONS,
                                          "1930-1975", "1976-2017")),
         x = factor(STATIONS$abbr[j], levels = STATIONS$abbr),
         y = factor(STATIONS$abbr[i], levels = rev(STATIONS$abbr)))
save_tab(transmute(tau, group, station_1 = STATIONS$label[i],
                   station_2 = STATIONS$label[j], tau), "tab_07_kendall_tau")

p_tau <- ggplot(tau, aes(x, y, fill = tau)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", tau)), size = 2.6) +
  facet_wrap(~ group, nrow = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-0.5, 0.5), name = expression(tau)) +
  coord_fixed() +
  labs(title = "Kendall's tau between stations",
       subtitle = "Autumn (ON) is the most spatially coherent season; early summer (JJ) the least",
       x = NULL, y = NULL,
       caption = "D Drummondville, A L'Assomption, C Les Cèdres, S St-Charles, L Lennoxville, N Nicolet") +
  theme(panel.grid = element_blank(), panel.border = element_blank())
save_fig(p_tau, "fig_05_kendall_tau", 10, 5.6)

# ---- empirical tail dependence across thresholds (Fig. 4.4) ----------------------
q_grid <- seq(0.50, 0.95, by = 0.01)
tdc <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
  i <- pairs$i[k]; j <- pairs$j[k]
  cbind(pair = pairs$pair[k],
        tail_dependence_empirical(r[[U_COLS[i]]], r[[U_COLS[j]]], q_grid))
})) |>
  pivot_longer(c(lambda_U, lambda_L), names_to = "tail") |>
  mutate(pair = factor(pair, levels = pairs$pair),
         tail = factor(tail, c("lambda_U", "lambda_L"), c("upper", "lower")))
p_tdc <- ggplot(tdc, aes(q, value, colour = tail)) +
  geom_vline(xintercept = 0.8, colour = "grey70", linetype = 3) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ pair, ncol = 5) +
  scale_colour_manual(values = c(upper = "#B2182B", lower = "#2166AC")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Empirical tail-dependence coefficients vs threshold",
       subtitle = "Non-negligible lower-tail dependence rules out a pure extreme-value copula",
       x = "Threshold q", y = expression(hat(lambda)(q)), colour = NULL)
save_fig(p_tdc, "fig_05_tail_dependence_curves", 10, 6.2)

# ---- lambda by season (q = .80) and period (q = .85), Fig. 4.5 ------------------
lam_group <- function(dat, group, q) {
  do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
    i <- pairs$i[k]; j <- pairs$j[k]
    cbind(group = group, pair = lab_pair(i, j), i = i, j = j,
          tail_dependence_empirical(dat[[U_COLS[i]]], dat[[U_COLS[j]]], q))
  }))
}
lam <- rbind(
  do.call(rbind, lapply(SEASONS, function(s) lam_group(filter(r, season == s), s, 0.80))),
  lam_group(filter(r, period == "pre"), "1930-1975", 0.85),
  lam_group(filter(r, period == "post"), "1976-2017", 0.85)
)
save_tab(lam, "tab_07_tail_dependence_empirical")

# ---- extreme-value dependence test (Table 4.1) -----------------------------------
set.seed(SEED)
ev <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
  i <- pairs$i[k]; j <- pairs$j[k]
  x <- na.omit(as.matrix(r[U_COLS[c(i, j)]]))
  t <- copula::evTestK(copula::pobs(x))
  data.frame(station_1 = STATIONS$label[i], station_2 = STATIONS$label[j],
             n = nrow(x), statistic = unname(t$statistic), p_value = t$p.value)
}))
save_tab(ev, "tab_08_ev_dependence_test")
log_step("  evTestK: ", sum(ev$p_value < 0.05), "/15 pairs reject EV dependence at 5%")
finish_step("03")
