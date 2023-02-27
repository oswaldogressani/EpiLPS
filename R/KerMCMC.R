#' Kernel for MCMC sampling
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerMCMC <- function(Dobs, BB, Pen, Covar, thetaoptim, penoptim, overdispoptim,
                    progress) {
  # Kernel routine
  # Author: Oswaldo Gressani (oswaldo_gressani@hotmail.fr)
  K <- ncol(BB)

  logtar <- function(zeta, lambda) {# Log of target
    theta <- zeta[1:K]
    w <- zeta[(K + 1)]
    rho <- exp(w)
    equal <- KerLikelihood(Dobs, BB)$loglik(theta, rho) - 0.5 * lambda *
      sum((theta * Pen) %*% theta) - 1e-04 * exp(w) + 1e-04 * w
    return(equal)
  }

  Dlogtar <- function(zeta, lambda) {# Gradient of log target
    theta <- zeta[1:K]
    w <- zeta[(K + 1)]
    rho <- exp(w)
    Btheta <- as.numeric(BB %*% theta)
    grad_theta <- KerLikelihood(Dobs, BB)$Dloglik(theta = theta, rho = rho) -
      lambda * as.numeric(Pen %*% theta)
    Dloglik_w <- sum(exp(w) * (digamma(Dobs + exp(w)) - digamma(exp(w)) +
                                 (1 + w) - (log(exp(Btheta) + exp(w)) +
                                              (1 / (1 + exp(Btheta - w))))) -
                       Dobs * (1 / (1 + exp(Btheta - w))))
    deriv_w <- Dloglik_w - 1e-04 * exp(w) + 1e-04
    res <- c(grad_theta, deriv_w)
    return(res)

  }

  MALA <- function(M) {# Metropolis-adjusted Langevin algorithm
    SigLH <- matrix(0, nrow = (K + 1), ncol = (K + 1))
    SigLH[1:K, 1:K] <- Covar
    SigLH[(K + 1), (K + 1)] <- 1
    counter <- 0
    zetamat <- matrix(0, nrow = M, ncol = (K + 1))
    lambvec <- c()
    deltavec <- c()

    # Initial values
    lambda_cur <- penoptim
    w_cur <- log(overdispoptim)
    theta_cur <- thetaoptim
    zeta_cur <- c(theta_cur, w_cur)
    tun <- 0.15

    for (m in 1:M) {
      # New proposal
      meanLH <- zeta_cur + 0.5 * tun * as.numeric(SigLH %*%
                                                    Dlogtar(zeta_cur, lambda_cur))
      zeta_prop <- as.numeric(Rcpp_KerMVN(mu = meanLH, Sigma = (tun * SigLH)))

      # Accept/Reject decision
      G_cur  <- Dlogtar(zeta_cur, lambda_cur)
      G_prop <- Dlogtar(zeta_prop, lambda_cur)
      ldiffq <- as.numeric((-0.5) * t(G_prop + G_cur) %*%
                             ((zeta_prop - zeta_cur) +
                                (tun / 4) * SigLH %*% (G_prop - G_cur)))
      ldiffp <- logtar(zeta_prop, lambda_cur) - logtar(zeta_cur, lambda_cur)
      logr <- ldiffp + ldiffq

      if (logr >= 0) {
        zetamat[m, ] <- zeta_prop
        counter <- counter + 1
        zeta_cur <- zeta_prop
      } else if (logr < 0) {
        u <- stats::runif(1)
        if (u <= exp(logr)) {
          zetamat[m, ] <- zeta_prop
          counter <- counter + 1
          zeta_cur <- zeta_prop
        } else{
          zetamat[m,] <- zeta_cur
        }
      }

      # Gibbs step for delta
      gdelta_shape <- 0.5 * 2 + 10
      gdelta_rate <- 0.5 * 2 * lambda_cur + 10
      deltavec[m] <- stats::rgamma(n = 1,
                                   shape = gdelta_shape,
                                   rate = gdelta_rate)

      # Gibbs step for lambda
      glambda_shape <- 0.5 * (K + 2)
      glambda_rate <- 0.5 * (as.numeric(t(zetamat[m, ][1:K]) %*% Pen
                                      %*% zetamat[m, ][1:K]) + deltavec[m] * 2)
      lambvec[m] <- stats::rgamma(n = 1, shape = glambda_shape,
                                  rate = glambda_rate)
      lambda_cur <- lambvec[m]

      # Automatic tuning of Langevin algorithm
      accept_prob <- min(c(1, exp(logr)))
      heval <- sqrt(tun) + (1 / m) * (accept_prob - 0.57)

      hfun <- function(x) {
        epsil <- 1e-04
        Apar <- 10 ^ 4
        if (x < epsil) {
          val <- epsil
        } else if (x >= epsil && x <= Apar) {
          val <- x
        } else{
          val <- Apar
        }
        return(val)
      }

      tun <- (hfun(heval)) ^ 2
      utils::setTxtProgressBar(progress, m)
    }
    close(progress)

    accept_rate <- round(counter / M * 100, 2)

    outlist <- list(
      theta_mcmc   = zetamat[, (1:K)],
      rho_mcmc     = exp(zetamat[, (K + 1)]),
      lambda_mcmc  = lambvec,
      delta_mcmc = deltavec,
      accept_rate  = accept_rate,
      niter = M
    )

    return(outlist)
  }

  outlist <- list(logtar = logtar, Dlogtar = Dlogtar, MALA = MALA)
  return(outlist)
}
