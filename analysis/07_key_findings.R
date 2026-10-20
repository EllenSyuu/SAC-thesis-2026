# Step 7 - one-panel summary of the thesis' main finding (README hero figure)
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("07 key-findings figure")

desc <- readr::read_csv(file.path(TAB, "tab_03_seasonal_summary.csv"),
                        show_col_types = FALSE) |>
  mutate(season = factor(season, SEASONS),
         period = factor(period, c("pre", "post"), c("1930-1975", "1976-2017")))
eta <- readr::read_csv(file.path(TAB, "tab_16_extremal_coefficient_6d.csv"),
                       show_col_types = FALSE) |>
  filter(group %in% SEASONS) |>
  mutate(season = factor(group, SEASONS))

pa <- ggplot(desc, aes(season, mean_mm, fill = period)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = mean_mm - 1.96 * se_mm, ymax = mean_mm + 1.96 * se_mm),
                position = position_dodge(0.75), width = 0.2, colour = "grey30") +
  scale_fill_manual(values = unname(PAL_PERIOD)) +
  coord_cartesian(ylim = c(15, 34)) +
  labs(title = "MAGNITUDE: how big are extremes?",
       subtitle = "Mean monthly maximum at a station. Late-summer and\nautumn extremes grew after 1975",
       x = NULL, y = "mm / day", fill = NULL)

pb <- ggplot(eta, aes(season, eta, fill = season)) +
  geom_hline(yintercept = 6, linetype = 2, colour = "grey50") +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", eta)), vjust = -0.4, size = 3.6) +
  annotate("text", x = 4.5, y = 6, label = "independence", hjust = 1,
           vjust = -0.5, size = 3, colour = "grey40") +
  scale_fill_manual(values = PAL_SEASON) +
  coord_cartesian(ylim = c(1, 6.6)) +
  labs(title = "CO-OCCURRENCE: do they hit together?",
       subtitle = "6-station extremal coefficient (1 = all together,\n6 = independent). Driven by season, not era",
       x = NULL, y = expression(hat(eta) == 6 * hat(A)(omega[0])))

hero <- (pa | pb) + plot_annotation(
  title = "Southern Québec precipitation extremes, 1930-2017",
  subtitle = "Non-stationarity in the margins and in the dependence structure are different stories",
  theme = theme(plot.title = element_text(face = "bold", size = 15)))
save_fig(hero, "fig_00_key_findings", 10.5, 4.6)
log_step("07 done")
