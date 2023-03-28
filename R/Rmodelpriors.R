#' Prior specification for model hyperparameters
#'
#' @description
#' Specification of the hyperparameters for the Gamma prior on the roughness
#' penalty parameter associated to the P-spline model and the Gamma prior on the
#' overdispersion parameter of the negative binomial model underlying the
#' incidence data.
#'
#' @param listcontrol A list specifying the hyperparameters in the Gamma priors
#' for the roughness penalty parameter of the P-spline model (named
#' \eqn{a_{\delta}}, \eqn{b_{\delta}} and \eqn{\phi}) and the overdispersion
#' parameter of the negative binomial model for the incidence data
#' (named \eqn{a_{\rho}} and \eqn{b_{\rho}}).
#'
#' @return A list with the specified hyperparameter components. By default,
#' \eqn{a_{\delta} = b_{\delta} = 10}, \eqn{\phi = 2} and
#'  \eqn{a_{\rho}=b_{\rho}=10^{-4}}.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Gressani, O., Wallinga, J., Althaus, C. L., Hens, N. and Faes, C.
#' (2022). EpiLPS: A fast and flexible Bayesian tool for estimation of the
#' time-varying reproduction number. \emph{Plos Computational Biology},
#' \strong{18}(10): e1010618.
#'
#' @export

Rmodelpriors <- function(listcontrol = list(a_delta = 10, b_delta = 10, phi = 2,
                                        a_rho = 1e-04, b_rho = 1e-04)) {
  return(listcontrol)
}
