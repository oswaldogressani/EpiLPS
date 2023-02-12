#' Posterior of B-spline parameter vector and related routines
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerPtheta <- function(Dobs, BB, Pen) {
  # Kernel routine
  # Author: Oswaldo Gressani (oswaldo_gressani@hotmail.fr)
  logptheta <- function(theta, rho, lambda) {# logpost theta
    equal <- KerLikelihood(Dobs, BB)$loglik(theta, rho) - 0.5 *
      lambda * sum((theta * Pen) %*% theta)
    return(equal)
  }
  Dlogptheta <- function(theta, rho, lambda) {# Gradient
    equal <- KerLikelihood(Dobs, BB)$Dloglik(theta, rho) -
      lambda * as.numeric(Pen %*% theta)
    return(equal)
  }
  D2logptheta <- function(theta, rho, lambda) {# Hessian
    equal <- KerLikelihood(Dobs, BB)$D2loglik(theta, rho) - lambda * Pen
    return(equal)
  }

  outlist <- list(logptheta = logptheta, Dlogptheta = Dlogptheta,
                  D2logptheta = D2logptheta)
  return(outlist)
}
