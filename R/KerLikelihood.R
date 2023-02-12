#' Likelihood and related routines
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerLikelihood <- function(Dobs, BB) {
  # Kernel routine
  # Author: Oswaldo Gressani (oswaldo_gressani@hotmail.fr)
  n <- length(Dobs)
  loglik <- function(theta, rho) {# Log-likelihood
    Btheta <- as.numeric(BB %*% theta)
    equal <-  sum(lgamma(rho + Dobs) - lgamma(rho)) +
      sum(Dobs * Btheta) + n * (rho * log(rho)) -
      sum((Dobs + rho) * log(rho + exp(Btheta)))
    return(equal)
  }
  Dloglik <- function(theta, rho) {# Gradient
    Btheta <- as.numeric(BB %*% theta)
    res <-
      colSums((Dobs - (Dobs + rho) / (1 + rho * exp(-Btheta))) * BB)
    return(res)
  }
  D2loglik <- function(theta, rho) {# Hessian
    Btheta <- as.numeric(BB %*% theta)
    midvec <-
      rho * (Dobs + rho) * (exp(Btheta) / (exp(Btheta) + rho) ^ 2)
    Hess <- (-1) * (t(BB) %*% (midvec * BB))
    return(Hess)
  }

  outlist <- list(loglik = loglik, Dloglik = Dloglik, D2loglik = D2loglik)
  return(outlist)
}
