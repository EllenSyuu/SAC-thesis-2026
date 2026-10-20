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

6. **Homogeneous simulation (Table 5.1).** The script behind Table 5.1 is not in the project archive. The rebuilt study (θ = 0.4, ρ = 4.5, n = 1000, B = 100) reproduces the ranking MLE < CFG < Pickands, but its absolute ISE values are smaller than those printed in the thesis.

7. **Period effect.** The histogram of β̂₄(ω) is centred near zero (median ≈ −0.006). Its tail is longer than the thesis wording suggests, though: about 60 % of the 15 × 99 estimates fall within ±0.02, and a few pairs reach −0.14. The conclusion that season dominates period still holds; it is just less one-sided than stated.
