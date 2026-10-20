# Step 4 - Monte-Carlo study of the estimators of A (Chapter 5.2)
#   Study 1: homogeneous Clayton-Archimax, theta = 0.4, logistic rho = 4.5,
#            n = 1000, B = 100  -> Pickands vs CFG vs MLE
#   Study 2: binary covariate, rho0 = 1.5 / rho1 = 4.5, n = 2000, B = 50
#            -> group-wise Pickands / CFG vs pooled covariate MLE
#               (identity link, constrained) vs MLE (logit link)
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

par_lapply <- function(X, FUN) {
  if (.Platform$OS.type == "unix" && N_CORES > 1) {
    parallel::mclapply(X, FUN, mc.cores = N_CORES, mc.set.seed = TRUE)
  } else lapply(X, FUN)
}
RNGkind("L'Ecuyer-CMRG")

# ---- Study 1 ---------------------------------------------------------------------
set.seed(SEED)
A_true1 <- pickands_logistic(w, 4.5)
study1 <- par_lapply(seq_len(B1), function(b) {
  U <- pseudo_obs(rclayton_archimax(N1, THETA, 4.5))
  th <- estimate_theta(U[, 1], U[, 2], se = FALSE)$theta
  curves <- lapply(c(MLE = "mle", CFG = "cfg", Pickands = "pickands"),
                   function(m) fill_ends(pickands_curve(U, th, w_in, m)))
  list(theta = th,
       err = do.call(rbind, lapply(names(curves), function(m)
         data.frame(rep = b, method = m, ISE = ise(w, curves[[m]], A_true1),
                    IRAE = irae(w, curves[[m]], A_true1)))),
       curves = if (b == 1) curves else NULL)
})
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

# ---- Study 2 ---------------------------------------------------------------------
A0 <- pickands_logistic(w, 1.5)
A1 <- pickands_logistic(w, 4.5)
study2 <- par_lapply(seq_len(B2), function(b) {
  z <- stats::rbinom(N2, 1, 0.5)
  raw <- rbind(rclayton_archimax(sum(z == 0), THETA, 1.5),
               rclayton_archimax(sum(z == 1), THETA, 4.5))
  z <- sort(z)
  U <- pseudo_obs(raw)
  th <- estimate_theta(U[, 1], U[, 2], se = FALSE)$theta
  X <- cbind(intercept = 1, z = z)
  est <- list()
  for (g in 0:1) {
    Ug <- U[z == g, ]
    est[[paste0("Pickands.", g)]] <- fill_ends(pickands_curve(Ug, th, w_in, "pickands"))
    est[[paste0("CFG.", g)]] <- fill_ends(pickands_curve(Ug, th, w_in, "cfg"))
  }
  f_id <- pickands_cov_mle(U, X, th, w_in, link = "identity")
  f_lg <- pickands_cov_mle(U, X, th, w_in, link = "logit")
  for (g in 0:1) {
    est[[paste0("MLE.", g)]] <- fill_ends(f_id$A[, g + 1])
    est[[paste0("MLE (link).", g)]] <- fill_ends(f_lg$A[, g + 1])
  }
  err <- do.call(rbind, lapply(names(est), function(nm) {
    g <- as.integer(sub(".*\\.", "", nm))
    At <- if (g == 0) A0 else A1
    data.frame(rep = b, group = g, method = sub("\\.[01]$", "", nm),
               ISE = ise(w, est[[nm]], At), IRAE = irae(w, est[[nm]], At))
  }))
  list(err = err, curves = if (b == 1) est else NULL)
})
err2 <- do.call(rbind, lapply(study2, `[[`, "err"))
tab2 <- err2 |> group_by(group, method) |>
  summarise(mean_ISE = mean(ISE), mean_IRAE = mean(IRAE), .groups = "drop") |>
  arrange(group, mean_ISE)
save_tab(tab2, "tab_10_sim_binary_covariate", digits = 5)
save_tab(err2, "tab_10_sim_binary_covariate_replicates", digits = 5)
log_step("  study 2 done")

err2$method <- factor(err2$method, levels = c("Pickands", "CFG", "MLE", "MLE (link)"))
p2 <- ggplot(err2, aes(method, IRAE, fill = method)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  facet_wrap(~ group, labeller = labeller(group = c(`0` = "z = 0 (rho = 1.5)",
                                                    `1` = "z = 1 (rho = 4.5)"))) +
  scale_fill_manual(values = PAL_METHOD) +
  labs(title = sprintf("Binary covariate: IRAE over B = %d replicates", B2),
       x = NULL, y = "IRAE")
save_fig(p1 + p2 + plot_layout(widths = c(1, 1.3)), "fig_06_simulation_study", 11, 4.4)
finish_step("04")
