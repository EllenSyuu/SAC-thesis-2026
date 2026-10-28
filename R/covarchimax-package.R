#' covarchimax: Covariate-Dependent Clayton-Archimax Copulas for Spatial Extremes
#'
#' Semiparametric inference for Clayton-Archimax copulas
#' \deqn{C(u_1,\dots,u_d) = \psi_\theta\Big(\ell\big(\phi_\theta(u_1),\dots,
#'   \phi_\theta(u_d)\big)\Big),\qquad \psi_\theta(t) = (1+\theta t)^{-1/\theta},}
#' where \eqn{\ell} is a stable tail dependence function with Pickands
#' dependence function \eqn{A} on the unit simplex, and where \eqn{A} may
#' depend on covariates.
#'
#' All estimators of \eqn{A} rest on one identity (Chatelain et al., 2020):
#' for fixed \eqn{\omega} in the simplex, \eqn{\xi(\omega) =
#' \min_j \phi_\theta(U_j)/\omega_j} has survival function
#' \eqn{x \mapsto \psi_\theta\{x A(\omega)\}}, so \eqn{A(\omega)} is a scale
#' parameter of \eqn{\xi(\omega)} ([xi_transform()]).
#'
#' @section Main functions:
#' \describe{
#'   \item{Generator parameter}{[estimate_theta()], [estimate_theta_pairs()]}
#'   \item{Pickands dependence function}{[pickands_at()], [pickands_curve()],
#'     [extremal_coefficient()], [simplex_path()]}
#'   \item{Covariate model}{[pickands_cov_mle()], [pseudo_obs_by_group()]}
#'   \item{Tail dependence}{[tail_dependence_archimax()],
#'     [tail_dependence_empirical()]}
#'   \item{Simulation and accuracy}{[rclayton_archimax()],
#'     [pickands_logistic()], [ise()], [irae()]}
#'   \item{Uncertainty}{[boot_years()], [boot_ci()]}
#'   \item{Margins}{[fit_gev()], [gev_backward()], [gev_coef_table()],
#'     [gev_sandwich_vcov()], [gev_to_gumbel()]}
#'   \item{Data checks}{[validate_daily()], [validate_monthly()],
#'     [monthly_maxima()]}
#' }
#'
#' @references
#' Capéraà, P., Fougères, A.-L. and Genest, C. (1997). A nonparametric
#' estimation procedure for bivariate extreme value copulas.
#' \emph{Biometrika}, 84(3), 567--577. \doi{10.1093/biomet/84.3.567}
#'
#' Charpentier, A., Fougères, A.-L., Genest, C. and Nešlehová, J. G. (2014).
#' Multivariate Archimax copulas. \emph{Journal of Multivariate Analysis},
#' 126, 118--136. \doi{10.1016/j.jmva.2013.12.013}
#'
#' Chatelain, S., Fougères, A.-L. and Nešlehová, J. G. (2020). Inference for
#' Archimax copulas. \emph{The Annals of Statistics}, 48(2), 1025--1051.
#' \doi{10.1214/19-AOS1836}
#'
#' Mhalla, L., Chavez-Demoulin, V. and Naveau, P. (2017). Non-linear models
#' for extremal dependence. \emph{Journal of Multivariate Analysis}, 159,
#' 49--66. \doi{10.1016/j.jmva.2017.04.006}
#'
#' Zhou, Z. (2026). \emph{Semiparametric Analysis of Precipitation Extremes in
#' Southern Québec using Archimax Copulas}. M.Sc. thesis, McGill University.
#'
#' @keywords internal
"_PACKAGE"
