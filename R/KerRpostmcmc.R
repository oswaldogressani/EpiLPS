#' Kernel for R output based on MCMC sample
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerRpostmcmc <- function(BB, sinter, thetasample, Tdom){
  # Kernel routine
  # Author: Oswaldo Gressani (oswaldo_gressani@hotmail.fr)
  K <- ncol(BB)
  n <- length(Tdom)
  simax <- length(sinter)
  tt <- seq_len(n)
  nitereff <- nrow(thetasample)

  Rsampling <- function(t, theta) {
    musample <- as.numeric(exp(BB %*% theta))
    if (t == 1) {
      val <- musample[t]
    } else if (t >= 2 && t <= simax) {
      val <- musample[t] * ((sum(rev(musample[1:(simax - 1)][1:(t - 1)]) *
                                   sinter[1:(simax - 1)][1:(t - 1)])) ^ (-1))
    } else if (t > simax && t <= n) {
      val <- musample[t] * (sum(rev(musample[(t - simax):(t - 1)]) *
                                  sinter) ^ (-1)) * (t > simax && t <= n)
    }
    return(val)
  }

  # Obtaining the samples for R at each time point
  Rsamplemat <- matrix(0, nrow = nitereff, ncol = n)
  for (j in 1:nitereff) {
    Rsamplemat[j, ] <- sapply(tt, Rsampling, theta = thetasample[j,])
  }

  # Create the output matrix
  Time <- Tdom
  R <- colMeans(Rsamplemat)
  Rsd <- apply(Rsamplemat, MARGIN = 2, FUN = stats::sd)
  Rq0.025 <- apply(Rsamplemat, MARGIN = 2, FUN = stats::quantile, probs = 0.025)
  Rq0.05  <-  apply(Rsamplemat, MARGIN = 2, FUN = stats::quantile, probs = 0.05)
  Rq0.25  <- apply(Rsamplemat, MARGIN = 2, FUN = stats::quantile, probs = 0.25)
  Rq0.50  <- apply(Rsamplemat, MARGIN = 2, FUN = stats::quantile, probs = 0.50)
  Rq0.75  <- apply(Rsamplemat, MARGIN = 2, FUN = stats::quantile, probs = 0.75)
  Rq0.95  <- apply(Rsamplemat, MARGIN = 2, FUN = stats::quantile, probs = 0.95)
  Rq0.975 <- apply(Rsamplemat, MARGIN = 2, FUN = stats::quantile, probs = 0.975)

  out <- data.frame(Time = Time, R = R, Rsd = Rsd, Rq0.025 = Rq0.025,
                    Rq0.05 = Rq0.05, Rq0.25 = Rq0.25, Rq0.50 = Rq0.50,
                    Rq0.75 = Rq0.75, Rq0.95 = Rq0.95, Rq0.975 = Rq0.975)

  return(out)
}
