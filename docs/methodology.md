# Methodology

A compact version of Chapters 2 and 5 of the thesis. Each section points to the function that implements it.

## 1. Margins: non-stationary GEV

For station $j$ and month $i$, the monthly maximum is modelled as

$$
M_{ij} \sim \mathrm{GEV}(\mu_{ij}, \sigma_i, \xi), \qquad
\mu_{ij} = \beta_0 + \beta_{\text{post}} \mathbb 1_{\text{post-1975}} + \sum_{s\in\{JJ,AS,ON\}} \beta_s \mathbb 1_{s},\qquad
\sigma_i = \gamma_0 + \sum_s \gamma_s \mathbb 1_s .
$$

The structure is selected per station by backward elimination with likelihood-ratio tests, respecting marginality (`gev_backward()`), starting from all main effects and two-way interactions of `anomaly`, `period` and `season`. Both elimination orders (scale first, location first) are run as a robustness check. In the pooled model, $\mu \sim \text{period}\times\text{station} + \text{season}$. Because the six stations are dependent, its likelihood is a composite likelihood, so its standard errors are optimistic. Implemented by `fit_gev()`, `gev_coef_table()`, `gev_ic()` and `gev_mean()`.

## 2. Dependence: the Clayton–Archimax copula

$$
C_{\theta,A}(u_1,\dots,u_d) = \psi_\theta\!\Big( \textstyle\sum_j \phi_\theta(u_j)\; A\big(\tfrac{\phi_\theta(u_1)}{S},\dots,\tfrac{\phi_\theta(u_d)}{S}\big)\Big),
\quad \psi_\theta(t) = (1+\theta t)^{-1/\theta},\; S = \textstyle\sum_j\phi_\theta(u_j).
$$

The Pickands function $A$ governs extremal dependence. It satisfies $\max_j\omega_j \le A(\omega)\le 1$, with the upper bound meaning independence. The Clayton generator adds lower-tail dependence below the asymptotic regime. In the bivariate case:

$$
\lambda_U = 2-2A(\tfrac12), \qquad \lambda_L = \{2A(\tfrac12)\}^{-1/\theta}
\qquad(\texttt{tail\_dependence\_archimax()}).
$$

## 3. The key identity

For a fixed $\omega$ on the unit simplex, let $\xi(\omega) = \min_j \phi_\theta(U_j)/\omega_j$ (`xi_transform()`). Then

$$
\Pr\{\xi(\omega) > x\} = \psi_\theta\{x\,A(\omega)\}
\;\Longrightarrow\; \xi(\omega) \overset{d}{=} Z / A(\omega),\qquad \Pr(Z>x)=\psi_\theta(x).
$$

So $A(\omega)$ is simply a scale parameter of $\xi(\omega)$. With $\mathbb E Z = 1/(1-\theta)$ and $\mathbb E\log Z = \Psi(1)-\Psi(1/\theta)-\log\theta$, this gives three estimators, implemented in `pickands_at()` and `pickands_curve()`:

| Estimator | Formula | Notes |
|---|---|---|
| Pickands-type | $\hat A = \mathbb E Z / \overline{\hat\xi}$ | requires $\theta<1$; poor when $1/\theta \le 2$ |
| CFG-type | $\hat A = \exp\{\mathbb E\log Z - \overline{\log\hat\xi}\}$ | |
| **MLE (new)** | $\hat A=\arg\max_{a\in[\max\omega,1]} n\log a - (1/\theta+1)\sum_i\log(1+\theta a\hat\xi_i)$ | always inside the admissible band |

## 4. Covariates in $A$ (the thesis's main methodological contribution)

At each fixed $\omega$, the model is $A(\omega\mid z_i) = z_i^\top\beta(\omega)$. The pseudo-log-likelihood is maximised subject to $\max\omega \le z^\top\beta(\omega) \le 1$ for every observed covariate profile:

$$
\ell\{\beta(\omega)\} = \sum_i \log a_i - (1/\theta+1)\log\{1+\theta a_i\hat\xi_i(\omega)\},\qquad a_i = z_i^\top\beta(\omega).
$$

There is also an unconstrained alternative, the logit link: $A(\omega\mid z) = m + (1-m)\,\mathrm{logit}^{-1}(z^\top\beta)$ with $m=\max\omega$. Both are implemented in `pickands_cov_mle(link = "identity" | "logit")`. For a discrete covariate, the pseudo-observations are re-ranked within each covariate cell (`pseudo_obs_by_group()`). This targets the conditional copula; otherwise, seasonal shifts in the *margins* would show up as spurious dependence.

## 5. The generator parameter $\theta$

The moment estimator of Chatelain et al. (2020, Eq. 7.4) is $\hat\theta = S_n/R_n$, where $S_n = 8m_{n,1}-9m_{n,2}-1$, $R_n = 1-4m_{n,1}+3m_{n,2}$, and $m_{n,k}$ are U-statistic estimates of the moments of the Kendall distribution. Standard errors come from a delete-one jackknife computed in $O(n^2)$ via closed-form updates (`estimate_theta()`). The global $\theta$ is the average of the 15 pairwise estimates (`estimate_theta_pairs()`).

## 6. Beyond pairs: $d = 6$

`pickands_at()` works on any face of the simplex. Two summaries are used:
- the **extremal coefficient** $\eta = d\,\hat A(\omega_0)$ at the centroid, which runs from 1 (complete dependence) to $d$ (independence) (`extremal_coefficient()`);
- **simplex paths** $\omega(t) = (1-t)\,\omega_{jk} + t\,\omega_0$, from a bivariate face to the centroid (`simplex_path()`).

## References

- Chatelain, S., Fougères, A.-L., Nešlehová, J. G. (2020). Inference for Archimax copulas. *Annals of Statistics* 48(2), 1025–1051.
- Mhalla, L., Chavez-Demoulin, V., Naveau, P. (2017). Non-linear models for extremal dependence. *Journal of Multivariate Analysis* 159, 49–66.
- Ben Ghorbal, N., Genest, C., Nešlehová, J. (2009). On the Ghoudi, Khoudraji, and Rivest test for extreme-value dependence. *Canadian Journal of Statistics* 37(4), 534–552.
- Charpentier, A., Fougères, A.-L., Genest, C., Nešlehová, J. G. (2014). Multivariate Archimax copulas. *Journal of Multivariate Analysis* 126, 118–136.
- Gilleland, E., Katz, R. W. (2016). extRemes 2.0. *Journal of Statistical Software* 72(8).
