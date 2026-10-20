# Reproducibility notes: thesis vs. this repository

The original analysis lived in roughly 140 exploratory R scripts with hard-coded local paths. This repository re-implements it as a tested package (`R/`) and a seven-step pipeline (`analysis/`). Every headline number was then re-derived from the raw ECCC download. This page records what matches, what differs, and why.

## Reproduced exactly

| Thesis item | Repository output | Status |
|---|---|---|
| Table 3.3 – sample counts per season × period | `tab_03_seasonal_summary.csv` | ✅ identical |
| Table 3.4 – means and SEs; Welch F = 51.12 / 46.78 | `tab_03_seasonal_summary.csv`, `tab_03_welch_anova_by_period.csv` | ✅ identical |
| Tables 3.5–3.17 – backward elimination, selected station models | `tab_04_backward_elimination*.csv` | ✅ same selected model at all six stations |
| Table 3.6 – Drummondville coefficients | `tab_05_gev_station_coefficients.csv` | ✅ identical to 3 d.p. |
| Table 3.19 – pooled GEV coefficients | `tab_06_gev_pooled_coefficients.csv` | ✅ identical to 3 d.p. |
| Table 4.1 – EV-dependence test rejects for all 15 pairs | `tab_08_ev_dependence_test.csv` | ✅ 15/15 rejected (see note 3) |
| Table 5.3 – pairwise θ, SE, 90 % CI; global θ = 0.35 | `tab_11_theta_pairwise.csv` | ✅ identical |
| Table 5.4 – season-only β(0.5) | `tab_14_season_beta_w05.csv` | ✅ identical to 3 d.p. |
| Table 5.5 – 6-d extremal coefficients | `tab_16_extremal_coefficient_6d.csv` | ✅ identical (η = 4.03 / 5.73 / 4.12 / 3.90 / 4.12) |
| Table 5.2 – ranking MLE (link) < MLE < CFG < Pickands | `tab_10_sim_binary_covariate.csv` | ✅ same ranking, same order of magnitude |

## Clarifications surfaced by the rebuild

1. **Pseudo-observations.** Section 4.1 describes ranking the Gumbel-standardised values from the pooled GEV fit. The numbers in Chapters 4–5 were, however, computed from ranks of the *raw* monthly maxima within each station. The two coincide only if the margins are stationary. Both versions are produced here: `pseudo_obs.csv` (raw, the default, which reproduces the thesis) and `pseudo_obs_gumbel.csv`. Steps 05–06 rerun on the latter (`make sensitivity`, suffix `_gumbel`), and the qualitative conclusions are unchanged (see the README).

2. **Re-ranking within covariate cells.** The season-specific estimates (Table 5.4, Figure 5.6, Table 5.5) re-rank the pseudo-observations inside each season. This matters: with ranks pooled over the whole year, every April–May maximum is low at *every* station. The group-wise MLE then puts AM at the comonotone bound and JJ / AS at independence, which is an artefact of the margins, not of dependence. It is now built into the package as `pseudo_obs_by_group()`.

3. **Standard errors of θ.** Table 5.3 cites Eq. (7.5) of Chatelain et al. (2020), but its SEs come from a delete-one jackknife. The plug-in version of Eq. (7.5) gives negative variance estimates for 6 of 15 pairs at this sample size. The package implements the jackknife in O(n²) with closed-form leave-one-out U-statistic updates.

4. **Common θ.** The thesis states that the global θ̂ = 0.35 lies inside all 15 pairwise 90 % intervals. With unrounded values it lies in 14: St-Charles–Lennoxville has the interval (0.087, 0.336) against θ̂ = 0.348. All 15 intervals still overlap each other, so the common-θ assumption remains reasonable.

5. **EV test pairing.** The original script aligned the two series of a pair by row index, so station series of different lengths were truncated rather than matched by date. Here pairs are matched by month. Test statistics move slightly; all 15 rejections stand.

6. **Homogeneous simulation (Table 5.1).** The script behind Table 5.1 is not in the project archive. The rebuilt study (θ = 0.4, ρ = 4.5, n = 1000, B = 100) reproduces Pickands as clearly the worst estimator, but **does not** reproduce a clear advantage of the MLE over CFG: their IRAE values are 0.0030 and 0.0032, and CFG has the lower ISE. The thesis's absolute ISE values are also larger than here. Study 2b adds that with within-group ranks, CFG is slightly more accurate than the MLE; the MLE's edge appears under pooled ranks.

7. **Period effect, now with intervals.** The histogram of β̂₄(ω) is centred near zero (median ≈ −0.006), but the year-block bootstrap (step 09) shows the period effect is not negligible everywhere. 29 of 135 (pair, ω) intervals exclude 0, against about 13 expected by chance, and 21 of the 29 come from the three pairs St-Charles–L'Assomption, St-Charles–Nicolet and St-Charles–Drummondville, all with stronger post-1975 dependence. For the other 12 pairs the thesis statement that period has a much weaker effect than season holds. The README states it in this qualified form.

8. **Season effect, now with intervals.** At the network level the June–July vs October–November difference in the 6-d extremal coefficient is 1.83, 90 % CI (1.24, 2.28). Pairwise, A_JJ(½) − A_ON(½) has an interval excluding 0 for 8 of 15 pairs. The other seven are the five Lennoxville pairs and Drummondville with L'Assomption and with St-Charles.

9. **Incomplete months (not in the thesis).** 34 of 4,174 station-months have fewer than 20 valid days. The thesis keeps them, and so does the default here (`MIN_DAYS = 0`). Step 08 refits the headline quantities for thresholds from 0 to 25 days: θ̂ and η̂ move by at most 0.01 and 0.11, while the Drummondville post-1975 location shift falls from 1.35 to 1.06 mm.

10. **Composite-likelihood SEs (not in the thesis).** The thesis notes that the pooled GEV SEs are optimistic but does not correct them. Godambe sandwich SEs are now reported. The season, scale and shape effects remain clearly non-zero; the post-1975 shift at Drummondville is borderline (95 % sandwich CI 0.02 to 2.67).

## Engineering conventions

- `analysis/.lintr` relaxes three style rules for plotting scripts: long label strings, one-line anonymous functions and compact pipes. Package code is linted with the stricter root `.lintr`. Both are enforced in CI.
- The cache in `results/cache/` is content-addressed: the key hashes the input data, the parameters and every file in `R/`. A stale cache is therefore recomputed automatically, never silently reused.
