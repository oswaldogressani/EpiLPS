#' Computing R summary statistics for epilps
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerRpostmap <- function(BB, theta, Covar, sinter, MVvec, Tdom) {
  # Kernel routine
  # Author: Oswaldo Gressani (oswaldo_gressani@hotmail.fr)
  K <- ncol(BB)
  muhat <- as.numeric(exp(BB %*% theta))
  n <- length(muhat)
  simax <- length(sinter)
  tt <- seq_len(n)

  # Point estimation of reproduction number at time point t
  Restim <- function(t) {
    if (t == 1) {
      val <- muhat[t]
    } else if (t >= 2 && t <= simax) {
      val <- muhat[t] * ((sum(rev(muhat[1:(simax - 1)][1:(t - 1)]) *
                                sinter[1:(simax - 1)][1:(t - 1)])) ^ (-1))
    } else if (t > simax && t <= n) {
      val <- muhat[t] * (sum(rev(muhat[(t - simax):(t - 1)]) *
                               sinter) ^ (-1)) * (t > simax && t <= n)
    }
    return(val)
  }
  R <- sapply(tt, Restim)

  # Posterior standard deviation of reproduction number at time point t
  Rsdpost <- function(t) {
    dhstar_t <- c()
    for (k in 1:K) {
      if (t == 1) {
        add <- 0
      } else if (t >= 2 && t <= simax) {
        add <- (-1) * ((sum(rev(muhat[1:(simax - 1)][1:(t - 1)]) *
                              sinter[1:(simax - 1)][1:(t - 1)])) ^ (-1)) *
          (sum(rev(muhat[1:(simax - 1)][1:(t - 1)]) *
                 sinter[1:(simax - 1)][1:(t - 1)] *
                 rev(BB[, k][1:(simax - 1)][1:(t - 1)])))
      } else if (t > simax && t <= n) {
        add <- (-1) * (sum(rev(muhat[(t - simax):(t - 1)]) * sinter) ^ (-1)) *
          sum(rev(muhat[(t - simax):(t - 1)]) *
                sinter * rev(BB[, k][(t - simax):(t - 1)]))
      }
      dhstar_t[k] <- BB[t, k] + add
    }
    val <- as.numeric(sqrt(MVvec[t] * t(dhstar_t) %*% Covar %*% dhstar_t))
    return(val)
  }
  meanlog <- log(R)
  sdlog <- sapply(tt, Rsdpost)
  Rsd <- sqrt(exp(sdlog ^ 2 - 1) * exp(2 * meanlog + sdlog ^ 2))
  Rq0.025 <- stats::qlnorm(0.025, meanlog = meanlog, sdlog = sdlog)
  Rq0.05  <- stats::qlnorm(0.05,  meanlog = meanlog, sdlog = sdlog)
  Rq0.25  <- stats::qlnorm(0.25,  meanlog = meanlog, sdlog = sdlog)
  Rq0.50  <- stats::qlnorm(0.50,  meanlog = meanlog, sdlog = sdlog)
  Rq0.75  <- stats::qlnorm(0.75,  meanlog = meanlog, sdlog = sdlog)
  Rq0.95  <- stats::qlnorm(0.95,  meanlog = meanlog, sdlog = sdlog)
  Rq0.975 <- stats::qlnorm(0.975, meanlog = meanlog, sdlog = sdlog)

  out <- data.frame(Time = Tdom, R = R, Rsd = Rsd, Rq0.025 = Rq0.025,
                    Rq0.05 = Rq0.05, Rq0.25 = Rq0.25, Rq0.50 = Rq0.50,
                    Rq0.75 = Rq0.75, Rq0.95 = Rq0.95, Rq0.975 = Rq0.975)
  return(out)
}
