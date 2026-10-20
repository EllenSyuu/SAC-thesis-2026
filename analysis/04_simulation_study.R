# Step 4 - Monte-Carlo study of the estimators of A (Chapter 5.2)
#   Study 1: homogeneous Clayton-Archimax, theta = 0.4, logistic rho = 4.5,
#            n = 1000, B = 100  -> Pickands vs CFG vs MLE
#   Study 2: binary covariate, rho0 = 1.5 / rho1 = 4.5, n = 2000, B = 50
#            -> group-wise Pickands / CFG vs pooled covariate MLE
#               (identity link, constrained) vs MLE (logit link)
#   Study 2b: rank-scheme check. Same copulas, but group z = 1 has its
#            margins shifted by a monotone transform (copula unchanged), as
#            seasons shift the level of precipitation. Pseudo-observations
#            are formed either by pooled ranks or by within-group ranks (the
#            scheme used for the real data in steps 05, 06 and 09).
#
# Every replicate draws from its own pre-generated seed, so results do not
# depend on the number of cores; results are cached (see _setup.R).
source(file.path(rprojroot::find_root("DESCRIPTION"), "analysis", "_setup.R"))
log_step("04 simulation study (QUICK = ", QUICK, ", cores = ", N_CORES, ")")

THETA <- 0.4
w <- seq(0, 1, by = 0.01)
w_in <- w[w > 0 & w < 1]
fill_ends <- function(A) c(1, A, 1)          # A(0) = A(1) = 1 by definition
B1 <- if (QUICK) 6 else 100
B2 <- if (QUICK) 4 else 50
N1 <- if (QUICK) 400 else 1000
N2 <- if (QUICK) 600 else 2000

# Replicate b always uses seeds[b]: identical results for any core count.
par_lapply <- function(B, FUN, seed) {
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  run <- function(b) {
    set.seed(seeds[b])
    FUN(b)
  }
  if (.Platform$OS.type == "unix" && N_CORES > 1) {
    parallel::mclapply(seq_len(B), run, mc.cores = N_CORES, mc.preschedule = FALSE)
  } else lapply(seq_len(B), run)
}

# ---- Study 1 ---------------------------------------------------------------------
A_true1 <- pickands_logistic(w, 4.5)
study1 <- cached("sim_study1", list(B1, N1, THETA, w, SEED), par_lapply(B1, seed = SEED, function(b) {
  U <- pseudo_obs(rclayton_archimax(N1, THETA, 4.5))
  th <- estimate_theta(U[, 1], U[, 2], se = FALSE)$theta
  curves <- lapply(c(MLE = "mle", CFG = "cfg", Pickands = "pickands"),
                   function(m) fill_ends(pickands_curve(U, th, w_in, m)))
  list(theta = th,
       err = do.call(rbind, lapply(names(curves), function(m)
         data.frame(rep = b, method = m, ISE = ise(w, curves[[m]], A_true1),
                    IRAE = irae(w, curves[[m]], A_true1)))),
       curves = if (b == 1) curves else NULL)
}))
err1 <- do.call(rbind, lapply(study1, `[[`, "err"))
th1  <- vapply(study1, `[[`, numeric(1), "theta")
tab1 <- err1 |> group_by(method) |>
  summarise(mean_ISE = mean(ISE), mean_IRAE = mean(IRAE)) |> arrange(mean_ISE)
tab1 <- rbind(tab1, data.frame(method = "theta_hat (mean, sd)",
                               mean_ISE = mean(th1), mean_IRAE = sd(th1)))
save_tab(tab1, "tab_09_sim_homogeneous")
log_step("  study 1 done: mean theta_hat = ", round(mean(th1), 3))

c1 <- study1[[1]]$curves
df1 <- rbind(data.frame(w = w, A = A_true1, method = "True"),
             do.call(rbind, lapply(names(c1), function(m)
               data.frame(w = w, A = c1[[m]], method = m))))
p1 <- ggplot(df1, aes(w, A, colour = method, linewidth = method == "True")) +
  geom_line(aes(y = pmax(w, 1 - w)), colour = "grey75", linetype = 3,
            linewidth = 0.4) +
  geom_line() +
  scale_linewidth_manual(values = c(`TRUE` = 1.1, `FALSE` = 0.6), guide = "none") +
  scale_colour_manual(values = PAL_METHOD) +
  labs(title = "One replicate: estimators vs true logistic A",
       x = expression(omega), y = expression(A(omega)), colour = NULL)

# ---- Study 2 / 2b ------------------------------------------------------------------
A0 <- pickands_logistic(w, 1.5)
A1 <- pickands_logistic(w, 4.5)
w_coarse <- seq(0.02, 0.98, by = 0.02)        # 2b uses a coarser grid (cost)
shift_margins <- function(x) x^(1 / 3)          # monotone: copula unchanged

study2 <- cached("sim_study2", list(B2, N2, THETA, w, w_coarse, SEED + 1),
                 par_lapply(B2, seed = SEED + 1, function(b) {
  z <- stats::rbinom(N2, 1, 0.5)
  raw <- rbind(rclayton_archimax(sum(z == 0), THETA, 1.5),
               rclayton_archimax(sum(z == 1), THETA, 4.5))
  z <- sort(z)
  X <- cbind(intercept = 1, z = z)
  err_rows <- function(est, grid, setting) {
    do.call(rbind, lapply(names(est), function(nm) {
      g <- as.integer(sub(".*\\.", "", nm))
      At <- pickands_logistic(grid, if (g == 0) 1.5 else 4.5)
      data.frame(rep = b, setting = setting, group = g,
                 method = sub("\\.[01]$", "", nm),
                 ISE = ise(grid, est[[nm]], At), IRAE = irae(grid, est[[nm]], At))
    }))
  }
  fit_all <- function(U, th, grid, with_link) {
    gi <- grid[grid > 0 & grid < 1]
    fe <- function(A) c(1, A, 1)
    est <- list()
    for (g in 0:1) {
      Ug <- U[z == g, ]
      est[[paste0("Pickands.", g)]] <- fe(pickands_curve(Ug, th, gi, "pickands"))
      est[[paste0("CFG.", g)]] <- fe(pickands_curve(Ug, th, gi, "cfg"))
    }
    f_id <- pickands_cov_mle(U, X, th, gi, link = "identity")
    for (g in 0:1) est[[paste0("MLE.", g)]] <- fe(f_id$A[, g + 1])
    if (with_link) {
      f_lg <- pickands_cov_mle(U, X, th, gi, link = "logit")
      for (g in 0:1) est[[paste0("MLE (link).", g)]] <- fe(f_lg$A[, g + 1])
    }
    est
  }
  # Study 2 (thesis design): identical margins, pooled ranks
  U <- pseudo_obs(raw)
  th <- estimate_theta(U[, 1], U[, 2], se = FALSE)$theta
  main <- err_rows(fit_all(U, th, w, TRUE), w, "identical margins / pooled ranks")
  # Study 2b: shifted margins in group 1; theta from pooled ranks as in the data
  raw_s <- raw
  raw_s[z == 1, ] <- shift_margins(raw_s[z == 1, ])
  Up <- pseudo_obs(raw_s)
  th_s <- estimate_theta(Up[, 1], Up[, 2], se = FALSE)$theta
  Ug <- pseudo_obs_by_group(raw_s, z)
  wc <- c(0, w_coarse, 1)
  extra <- rbind(
    err_rows(fit_all(pseudo_obs_by_group(raw, z), th, wc, FALSE), wc,
             "identical margins / within-group ranks"),
    err_rows(fit_all(Up, th_s, wc, FALSE), wc, "shifted margins / pooled ranks"),
    err_rows(fit_all(Ug, th_s, wc, FALSE), wc, "shifted margins / within-group ranks"))
  list(err = rbind(main, extra))
}))
err2_all <- do.call(rbind, lapply(study2, `[[`, "err"))
err2 <- subset(err2_all, setting == "identical margins / pooled ranks")
tab2 <- err2 |> group_by(group, method) |>
  summarise(mean_ISE = mean(ISE), mean_IRAE = mean(IRAE), .groups = "drop") |>
  arrange(group, mean_ISE)
save_tab(tab2, "tab_10_sim_binary_covariate", digits = 5)
tab2b <- err2_all |> filter(method %in% c("CFG", "MLE")) |>
  group_by(setting, group, method) |>
  summarise(mean_ISE = mean(ISE), mean_IRAE = mean(IRAE), .groups = "drop")
save_tab(tab2b, "tab_10b_sim_rank_scheme", digits = 5)
save_tab(err2_all, "tab_10_sim_binary_covariate_replicates", digits = 5)
log_step("  study 2 / 2b done")

err2$method <- factor(err2$method, levels = c("Pickands", "CFG", "MLE", "MLE (link)"))
p2 <- ggplot(err2, aes(method, IRAE, fill = method)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  facet_wrap(~ group, labeller = labeller(group = c(`0` = "z = 0 (rho = 1.5)",
                                                    `1` = "z = 1 (rho = 4.5)"))) +
  scale_fill_manual(values = PAL_METHOD) +
  labs(title = sprintf("Binary covariate: IRAE over B = %d replicates", B2),
       x = NULL, y = "IRAE")
save_fig(p1 + p2 + plot_layout(widths = c(1, 1.3)), "fig_06_simulation_study", 11, 4.4)

p2b <- ggplot(transform(subset(err2_all, method == "MLE"),
                        setting = factor(setting, unique(err2_all$setting))),
              aes(setting, IRAE, fill = setting)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  facet_wrap(~ group, scales = "free_y",
             labeller = labeller(group = c(`0` = "z = 0 (rho = 1.5)",
                                           `1` = "z = 1 (rho = 4.5, margins shifted in 2b)"))) +
  scale_fill_manual(values = c("grey70", "grey85", "#D55E00", "#0072B2")) +
  coord_flip() +
  labs(title = "Why the real-data analysis re-ranks within covariate cells",
       subtitle = "Covariate MLE; a monotone shift of one group's margins leaves its copula unchanged",
       x = NULL, y = "IRAE")
save_fig(p2b, "fig_06b_sim_rank_scheme", 10, 4)
finish_step("04")
