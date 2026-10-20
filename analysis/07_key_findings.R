# Step 7 - one-panel summary of the thesis' main finding (README hero figure)
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("07 key-findings figure")

desc <- readr::read_csv(file.path(TAB, "tab_03_seasonal_summary.csv"),
                        show_col_types = FALSE) |>
  mutate(season = factor(season, SEASONS),
         period = factor(period, c("pre", "post"), c("1930-1975", "1976-2017")))
# point estimates and 90% year-block bootstrap intervals (step 09)
boot <- readr::read_csv(file.path(TAB, "tab_20_bootstrap_theta_eta.csv"),
                        show_col_types = FALSE)
eta <- boot |>
  filter(quantity %in% paste("eta", SEASONS)) |>
  mutate(season = factor(sub("eta ", "", quantity), SEASONS), eta = estimate)
d_jo <- boot[boot$quantity == "eta JJ - eta ON", ]
shift <- readr::read_csv(file.path(TAB, "tab_03_period_shift_by_season.csv"),
                         show_col_types = FALSE)

pa <- ggplot(desc, aes(season, mean_mm, fill = period)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = mean_mm - 1.96 * se_mm, ymax = mean_mm + 1.96 * se_mm),
                position = position_dodge(0.75), width = 0.2, colour = "grey30") +
  scale_fill_manual(values = unname(PAL_PERIOD)) +
  coord_cartesian(ylim = c(15, 34)) +
  labs(title = "MAGNITUDE: how big are extremes?",
       subtitle = sprintf(paste0("Mean monthly maximum, 95%% CI. Post-1975 minus pre-1975:\n",
                                 "AS %+.1f mm (p = %.3f), ON %+.1f mm (p < 0.001)"),
                          shift$diff_post_minus_pre[shift$season == "AS"],
                          shift$p_value[shift$season == "AS"],
                          shift$diff_post_minus_pre[shift$season == "ON"]),
       x = NULL, y = "mm / day", fill = NULL)

pb <- ggplot(eta, aes(season, eta, fill = season)) +
  geom_hline(yintercept = 6, linetype = 2, colour = "grey50") +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, colour = "grey25") +
  geom_text(aes(y = ci_upper, label = sprintf("%.2f", eta)), vjust = -0.5, size = 3.6) +
  annotate("text", x = 4.5, y = 6, label = "independence", hjust = 1,
           vjust = -0.5, size = 3, colour = "grey40") +
  scale_fill_manual(values = PAL_SEASON) +
  coord_cartesian(ylim = c(1, 6.8)) +
  labs(title = "CO-OCCURRENCE: do they hit together?",
       subtitle = sprintf(paste0("6-station extremal coefficient (1 = together, 6 = independent),\n",
                                 "90%% block-bootstrap CI. JJ - ON = %.2f (%.2f, %.2f)"),
                          d_jo$estimate, d_jo$ci_lower, d_jo$ci_upper),
       x = NULL, y = expression(hat(eta) == 6 * hat(A)(omega[0])))

hero <- (pa | pb) + plot_annotation(
  title = "Southern Québec precipitation extremes, 1930-2017",
  subtitle = "Magnitude and spatial co-occurrence of extremes respond to different covariates",
  theme = theme(plot.title = element_text(face = "bold", size = 15)))
save_fig(hero, "fig_00_key_findings", 10.5, 4.6)
finish_step("07")
