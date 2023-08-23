#' Estimation of the incubation density based on coarse data
#'
#' @description
#' This function computes an estimate of the incubation density based on
#' coarse data constructed from symptom onset times and exposure windows. It
#' uses the Laplacian-P-splines methodology with a Langevinized Gibbs algorithm
#' to sample from the posterior distribution.
#'
#' @usage estimIncub(x, K = 10, niter = 1000, tmax = max(x), tgridlen = 500, verbose = FALSE)
#'
#' @param x A data frame with the lower and upper bound of incubation interval.
#' @param K Number of B-splines in the basis.
#' @param niter The number of MCMC samples.
#' @param tmax The upper bound for the B-spline basis. Default is the largest
#'  data point in \code{x}.
#' @param tgridlen The number of grid points on which to evaluate the density.
#' @param verbose Should a message be printed? Default is FALSE.
#'
#' @return A list of class \code{incubestim} containing summary values for
#'  the estimated incubation density.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @examples
#' set.seed(123)
#' simdat <- incubsim(n = 30, tmax = 20) # Simulate incubation data
#' data <- simdat$Dobsincub              # Incubation bounds
#' incubfit <- estimIncub(x = data, niter = 500, tmax = 20, verbose = TRUE)
#'
#' @export

estimIncub <- function(x, K = 10, niter = 1000, tmax = max(x), tgridlen = 500, verbose = FALSE){

  tic <- proc.time()
  nobs <- nrow(x)
  Bl <- 0
  Br <- tmax
  tg <- seq(Bl, Br, length = tgridlen)
  dtg <- tg[2] - tg[1]

  #---- Partitioning the time domain
  nbins <- 300
  tie <- TRUE
  while(tie == TRUE){
    bins <- seq(Bl, Br, length = nbins + 1)
    binwidth <- bins[2] - bins[1]
    binmid <- bins[1:nbins] + binwidth * 0.5
    uptoL <- as.integer(x$tL / binwidth) + 1
    uptoL[which(uptoL == nbins + 1)] <- nbins
    uptoR <- as.integer(x$tR / binwidth) + 1
    uptoR[which(uptoR == nbins + 1)] <- nbins
    tie <- any(uptoL==uptoR)
    if(tie == TRUE){
      nbins <- nbins + 100
    }
  }

  #---- Penalty and B-spline basis construction
  D <- diag(K)
  penorder <- 3
  for(k in 1:penorder){D <- diff(D)}
  P <- t(D) %*% D           # Penalty matrix of dimension c(K,K)
  P <- P + diag(1e-06, K)   # Perturbation to ensure P is full rank
  Bmid <- Rcpp_KercubicBspline(binmid, lower = Bl, upper = Br, K = K)
  Btg <- Rcpp_KercubicBspline(tg, lower = Bl, upper = Br, K = K)

  #---- Priors on penalty parameter lambda
  a_pen <- 1e-04
  b_pen <- 1e-04

  #---- Log-likelihood, gradient and Hessian
  loglik <- function(theta){
    hmidbw <- as.numeric(exp(Bmid %*% theta)) * binwidth
    cmsum <- cumsum(hmidbw)
    SL <- exp(-cmsum[uptoL])
    SR <- exp(-cmsum[uptoR])
    val <- sum(log(SL-SR))
    return(val)
  }

  Dloglik <- function(theta) {
    hmidbw <- as.numeric(exp(Bmid %*% theta)) * binwidth
    cmsum <- cumsum(hmidbw)
    hBb <- hmidbw * Bmid
    hBbcmsum <- apply(hBb, 2, cumsum)
    SL <- exp(-cmsum[uptoL])
    SR <- exp(-cmsum[uptoR])
    DSLSR <- SL - SR
    psiR <- hBbcmsum[uptoR, ]
    psiL <- hBbcmsum[uptoL, ]
    grad <- colSums(((SR * psiR) - (SL * psiL)) * (1 / DSLSR))
    return(grad)
  }

  D2loglik <- function(theta){
    hessval <- matrix(0, nrow = K, ncol = K)
    hmidbw <- as.numeric(exp(Bmid %*% theta)) * binwidth
    cmsum <- cumsum(hmidbw)
    SL <- exp(-cmsum[uptoL])
    SR <- exp(-cmsum[uptoR])
    hBb <- hmidbw * Bmid
    hBbcmsum <- apply(hBb, 2, cumsum)
    psiR <- hBbcmsum[uptoR, ]
    psiL <- hBbcmsum[uptoL, ]
    zn <- (SR * psiR) - (SL * psiL)

    for(k in 1:K){
      psi_kL <- psiL[,k]
      psi_kR <- psiR[,k]
      znk <- zn[,k]
      crosscm <- apply(hBb * Bmid[, k], 2, cumsum)
      psicrossL <- crosscm[uptoL,]
      psicrossR <- crosscm[uptoR,]
      hessval[,k] <- colSums(((SL-SR)^(-1)) * (SR *   (psicrossR-(psiR*psi_kR))-
                        SL *  (psicrossL - (psiL*psi_kL)))-
                        znk * zn * ((SL-SR)^(-2)))
    }

    return(hessval)
  }

  #---- Gradient/Hessian of conditional log-posterior of B-spline coeffs.
  Dlogptheta <- function(theta, lambda){
    val <- as.numeric(Dloglik(theta) - lambda * (P%*%theta))
    return(val)
  }

  D2logptheta <- function(theta, lambda){
    val <- D2loglik(theta) - lambda * P
    return(val)
  }

  thetainit <- rep(0.5,K) # Initial value for theta vector

  logplamb <- function(lambda){
    LAcall <- Rcpp_KerLaplaceIncub(theta0 = thetainit, lambda, K,
                               Dlogptheta, D2logptheta)
    thetastar <- as.numeric(LAcall$Lapmode)

    val <- (0.5*K+a_pen - 1) * log(lambda) + 0.5 * Re(LAcall$logdetSigma) +
      loglik(thetastar) -
      lambda * (0.5 * sum((thetastar*P)%*%thetastar) + b_pen)
    return(val)

  }

  if(K < 15) {
    lamblow <- 10
  } else{
    lamblow <- 1.5
  }

  #----- Exploration of penalty over a grid
  loglambdas <- seq(log(lamblow), 4, length = 15)
  lamb_grid <- 10 ^ loglambdas
  logplamb_val <- sapply(lamb_grid, logplamb)
  if (anyNA(logplamb_val)) {
    nanpos <- which(is.na(logplamb_val))
    lamb_grid <- lamb_grid[-nanpos]
    logplamb_val <- logplamb_val[-nanpos]
  }
  lambmax <- lamb_grid[which(logplamb_val == max(logplamb_val))]

  if(lambmax == lamb_grid[1]) {
    penpos <- "boundary"
  } else if (lambmax == utils::tail(lamb_grid, 1)) {
    penpos <- "boundary"
  } else{
    penpos <- "interior"
  }

  maxplamb <- list(penpos = penpos, lambmax = lambmax)
  Lap_approx <- Rcpp_KerLaplaceIncub(theta0 = thetainit, lambmax, K,
                                 Dlogptheta, D2logptheta)
  thetaoptim <- as.numeric(Lap_approx$Lapmode)
  Sighat <- round(Lap_approx$Lapvar, 6)
  Sighat <- KerAdjustPD(Sighat)$PD

  #----------------------------------------------------------------------------#
  #                      MCMC-Metropolis-adjusted Langevin                     #
  #----------------------------------------------------------------------------#

  logtar <- function(theta, lambda){
    val <- loglik(theta) - 0.5 * lambda * sum((theta * P) %*% theta)
    return(val)
  }

  Dlogtar<- function(theta, lambda){
    val <- as.numeric(Dloglik(theta) - lambda * (P%*%theta))
    return(val)
  }

  MALA <- function(M) {# Metropolis-adjusted Langevin algorithm
    SigLH <- Sighat
    counter <- 0
    thetamat <- matrix(0, nrow = M, ncol = K)
    lambvec <- c()

    # Initial values
    lambda_cur <- lambmax
    theta_cur <- thetaoptim
    tun <- 0.15

    for (m in 1:M) {
      # New proposal
      meanLH <- theta_cur + 0.5 * tun * as.numeric(SigLH %*%
                                            Dlogtar(theta_cur, lambda_cur))
      theta_prop <- as.numeric(Rcpp_KerMVN(mu = meanLH, Sigma = (tun * SigLH)))

      # Accept/Reject decision
      G_cur  <- Dlogtar(theta_cur, lambda_cur)
      G_prop <- Dlogtar(theta_prop, lambda_cur)
      ldiffq <- as.numeric((-0.5) * t(G_prop + G_cur) %*%
                             ((theta_prop - theta_cur) +
                                (tun / 4) * SigLH %*% (G_prop - G_cur)))
      ldiffp <- logtar(theta_prop, lambda_cur) - logtar(theta_cur, lambda_cur)
      logr <- ldiffp + ldiffq

      if (logr >= 0) {
        thetamat[m, ] <- theta_prop
        counter <- counter + 1
        theta_cur <- theta_prop
      } else if (logr < 0) {
        u <- stats::runif(1)
        if (u <= exp(logr)) {
          thetamat[m, ] <- theta_prop
          counter <- counter + 1
          theta_cur <- theta_prop
        } else{
          thetamat[m,] <- theta_cur
        }
      }

      # Gibbs step for lambda
      glambda_shape <- 0.5 * K + a_pen
      glambda_rate <- 0.5 * as.numeric(t(thetamat[m, ]) %*% P
                                       %*% thetamat[m, ]) + b_pen
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
    }

    accept_rate <- round(counter / M * 100, 2)

    outlist <- list(
      theta_mcmc   = thetamat,
      lambda_mcmc  = lambvec,
      accept_rate  = accept_rate,
      niter = M,
      tunS = tun * SigLH
    )

    return(outlist)
  }

  MCMC <- MALA(niter)

  thetaMCMC <- MCMC$theta_mcmc
  penMCMC <- MCMC$lambda_mcmc
  thetahatMCMC <- apply(thetaMCMC, 2, "median")

  Incubfit <- function(theta) {
    htg <- as.numeric(exp(Btg %*% theta))
    Stg <- exp(-cumsum(htg * dtg))
    Itg <- htg * Stg
    Itg <- Itg/sum(Itg * dtg)
    return(Itg)
  }

  mixelem <- 2
  fmix <- matrix(0, nrow = length(tg), ncol = mixelem)
  xmid <- stats::qunif(p = 0.50, min = x$tL, max = x$tR)
  xfmid <- histosmooth(x = xmid, xl = 0, xr = tmax + 2)
  fmid <- sapply(tg, xfmid$fdens)
  fmid <- fmid / sum(fmid * dtg)
  fmix[, 1] <- Incubfit(thetahatMCMC)
  fmix[, 2] <- fmid

  # Computation of BIC
  ILL <- (-1) * D2loglik(thetahatMCMC)
  EDLPS <- sum(diag(solve(ILL + stats::median(penMCMC) * P) %*% ILL))
  BIC_LPS <- (-2) * loglik(thetahatMCMC) + EDLPS * log(nobs)

  #----------------------------------------------------------------------------#
  #                             Statistics for LPS                             #
  #----------------------------------------------------------------------------#
  qqgrid <- seq(0.05, 0.95, by = 0.05)
  qqlen <- length(qqgrid)

  LPSincub <- matrix(0, nrow = length(qqgrid) + 2, ncol = 5)
  colnames(LPSincub) <- c("PE","CI90L","CI90R","CI95L","CI95R")
  rownames(LPSincub) <- c("Mean","SD", paste0("q", qqgrid))

  ftgMCMC_LPS <- matrix(0, nrow = niter, ncol = length(tg))
  MeanLPS <- c()
  SDLPS <- c()
  qmatLPS <- matrix(0, nrow = niter, ncol = length(qqgrid))
  colnames(qmatLPS) <- paste0("q", qqgrid)

  for(j in 1:niter){
    fmixx <- fmix
    fmixx[,1] <- Incubfit(thetaMCMC[j,])
    ftgLPS <- rowSums(fmixx * (1/mixelem))
    ftgMCMC_LPS[j, ] <- ftgLPS
    # Statistics
    MeanLPS[j] <- sum(tg  * ftgLPS) * dtg
    SDLPS[j] <- sqrt(sum(((tg - MeanLPS[j])^2) * ftgLPS * dtg))
    cmsum <- cumsum(ftgLPS * dtg)
    for(l in 1:length(qqgrid)){
      qmatLPS[j,l] <- tg[utils::tail(which(cmsum<qqgrid[l]),1)]
    }
  }

  ftgLPS <- rowSums(fmix * (1/mixelem))
  EftgLPS <- sum(tg * ftgLPS) * dtg
  SDftgLPS <- sqrt(sum(((tg - EftgLPS)^2) * ftgLPS * dtg))
  cmsum <- cumsum(ftgLPS * dtg)
  qftgLPS <- c()
  for(l in 1:qqlen){
    qftgLPS[l] <- tg[utils::tail(which(cmsum < qqgrid[l]), 1)]
  }
  LPSincub[1,] <- c(EftgLPS,
                    stats::quantile(MeanLPS, probs = c(0.05,0.95,0.025,0.975)))
  LPSincub[2,] <- c(SDftgLPS,
                    stats::quantile(SDLPS, probs = c(0.05,0.95,0.025,0.975)))

  for(j in 3:nrow(LPSincub)){
    LPSincub[j,] <- c(qftgLPS[j-2],
                    stats::quantile(qmatLPS[,j-2], probs = c(0.05,0.95,0.025,
                                                             0.975)))
  }

  qMeanLPS <- stats::quantile(MeanLPS, probs = c(0.05, 0.95, 0.025, 0.975))
  qSDLPS <- stats::quantile(SDLPS, probs = c(0.05, 0.95, 0.025, 0.975))


  #---------------------------------------------------------------------------#
  #             Moment matching to Log-Normal, Gamma and Weibull              #
  #---------------------------------------------------------------------------#


  #---------------------------------- Log-Normal fit ---------------------------

  LN_mean_sd <- matrix(0, nrow = niter, ncol = 2)
  colnames(LN_mean_sd) <- c("meanlog","sdlog")
  qmatLN <- matrix(0, nrow = niter, ncol = length(qqgrid))
  colnames(qmatLN) <- paste0("q", qqgrid)
  ftgMCMC_LN <- matrix(0, nrow = niter, ncol = length(tg))

  for(j in 1:niter){
    # Moment matching
    dLN <- EpiLPS::Idist(mean = MeanLPS[j], sd = SDLPS[j], dist = "lognorm")
    LN_mean_sd[j,1] <- dLN$location
    LN_mean_sd[j,2] <- dLN$scale
    ftgMCMC_LN[j, ] <- stats::dlnorm(tg, meanlog = dLN$location,
                                     sdlog = dLN$scale)
    for(l in 1:qqlen){
      qmatLN[j,l] <- stats::qlnorm(p = qqgrid[l],
                            meanlog = dLN$location, sdlog = dLN$scale)
    }
  }

  meanlogMoM <- stats::median(LN_mean_sd[,1])
  sdlogMoM <- stats::median(LN_mean_sd[,2])
  ftgLN <- stats::dlnorm(tg, meanlog = meanlogMoM, sdlog = sdlogMoM)

  EftgLN <- exp(meanlogMoM + 0.5 * (sdlogMoM ^ 2))
  SDftgLN <- sqrt(exp(2 * meanlogMoM + sdlogMoM ^ 2) *
                    (exp(sdlogMoM ^ 2) - 1))
  qftgLN <- stats::qlnorm(p = qqgrid, meanlog = meanlogMoM, sdlog = sdlogMoM)

  LNincub <- matrix(0, nrow = length(qqgrid) + 2, ncol = 5)
  colnames(LNincub) <- c("PE","CI90L","CI90R","CI95L","CI95R")
  rownames(LNincub) <- c("Mean","SD", paste0("q", qqgrid))

  LNincub[1, ] <- c(EftgLN, qMeanLPS)
  LNincub[2,] <- c(SDftgLN, qSDLPS)

  for(j in 3:nrow(LNincub)){
    LNincub[j,] <- c(qftgLN[j-2],
                    stats::quantile(qmatLN[,j-2], probs = c(0.05,0.95,0.025,
                                                            0.975)))
  }

  LNloglik <- c()
  for(j in 1:nobs){
    LNloglik[j] <- log(stats::integrate(stats::dlnorm,
                                        meanlog = meanlogMoM, sdlog = sdlogMoM,
                                 lower = x$tL[j], upper = x$tR[j])$value)
  }

  BIC_LN <- (-2) * sum(LNloglik) + 2 * log(nobs)


  #---------------------------------- Gamma fit -------------------------------

  G_shape_rate <- matrix(0, nrow = niter, ncol = 2)
  colnames(G_shape_rate) <- c("shape","rate")
  qmatG <- matrix(0, nrow = niter, ncol = length(qqgrid))
  colnames(qmatG) <- paste0("q", qqgrid)
  ftgMCMC_G <- matrix(0, nrow = niter, ncol = length(tg))

  for(j in 1:niter){
    # Moment matching
    dG <- EpiLPS::Idist(mean = MeanLPS[j], sd = SDLPS[j], dist = "gamma")
    G_shape_rate[j,1] <- dG$shape
    G_shape_rate[j,2] <- dG$rate
    ftgMCMC_G[j, ] <- stats::dgamma(tg, shape = dG$shape, rate = dG$rate)
    for(l in 1:qqlen){
      qmatG[j,l] <- stats::qgamma(p = qqgrid[l], shape = dG$shape,
                                  rate = dG$rate)
    }
  }
  shapeGMoM <- stats::median(G_shape_rate[,1])
  rateMoM <- stats::median(G_shape_rate[,2])
  ftgG <- stats::dgamma(tg, shape = shapeGMoM, rate = rateMoM)
  EftgG <- shapeGMoM / rateMoM
  SDftgG <- sqrt(shapeGMoM / (rateMoM ^ 2))
  qftgG <- stats::qgamma(p = qqgrid, shape = shapeGMoM, rate = rateMoM)

  Gincub <- matrix(0, nrow = length(qqgrid) + 2, ncol = 5)
  colnames(Gincub) <- c("PE","CI90L","CI90R","CI95L","CI95R")
  rownames(Gincub) <- c("Mean","SD", paste0("q", qqgrid))
  Gincub[1,] <- c(EftgG, qMeanLPS)
  Gincub[2,] <- c(SDftgG, qSDLPS)

  for(j in 3:nrow(Gincub)){
    Gincub[j,] <- c(qftgG[j-2],
              stats::quantile(qmatG[,j-2], probs = c(0.05,0.95,0.025, 0.975)))
  }

  Gloglik <- c()
  for(j in 1:nobs){
    Gloglik[j] <- log(stats::integrate(stats::dgamma, shape = shapeGMoM,
                                rate = rateMoM, lower = x$tL[j],
                                upper = x$tR[j])$value)
  }
  BIC_G <- (-2) * sum(Gloglik) + 2 * log(nobs)

  #---------------------------------- Weibull fit -----------------------------

  W_shape_scale <- matrix(0, nrow = niter, ncol = 2)
  colnames(W_shape_scale) <- c("shape","scale")
  qmatW <- matrix(0, nrow = niter, ncol = length(qqgrid))
  colnames(qmatW) <- paste0("q", qqgrid)
  ftgMCMC_W <- matrix(0, nrow = niter, ncol = length(tg))

  for(j in 1:niter){
    # Moment matching
    dW <- EpiLPS::Idist(mean = MeanLPS[j], sd = SDLPS[j], dist = "weibull")
    W_shape_scale[j,1] <- dW$shape
    W_shape_scale[j,2] <- dW$scale
    ftgMCMC_W[j, ] <- stats::dweibull(tg, shape = dW$shape, scale = dW$scale)
    for(l in 1:qqlen){
      qmatW[j,l] <- stats::qweibull(p = qqgrid[l], shape = dW$shape,
                                    scale = dW$scale)
    }
  }

  shapeWMoM <- stats::median(W_shape_scale[,1])
  scaleMoM <- stats::median(W_shape_scale[,2])
  ftgW <- stats::dweibull(tg, shape = shapeWMoM, scale = scaleMoM)
  EftgW <- scaleMoM * gamma(1 + 1 / shapeWMoM)
  SDftgW <- sqrt((scaleMoM ^ 2) * ((gamma(1 + 2 / shapeWMoM) - (gamma(1 + 1 / shapeWMoM) ^ 2))))
  qftgW <- stats::qweibull(p = qqgrid, shape = shapeWMoM, scale = scaleMoM)

  Wincub <- matrix(0, nrow = length(qqgrid) + 2, ncol = 5)
  colnames(Wincub) <- c("PE","CI90L","CI90R","CI95L","CI95R")
  rownames(Wincub) <- c("Mean","SD", paste0("q", qqgrid))
  Wincub[1,] <- c(EftgW, qMeanLPS)
  Wincub[2,] <- c(SDftgW, qSDLPS)

  for(j in 3:nrow(Wincub)){
    Wincub[j,] <- c(qftgW[j-2],
                  stats::quantile(qmatW[,j-2], probs = c(0.05,0.95,0.025,
                                                         0.975)))
  }

  Wloglik <- c()
  for(j in 1:nobs){
    Wloglik[j] <- log(stats::integrate(stats::dweibull,
                                       shape = shapeWMoM, scale = scaleMoM,
                                lower = x$tL[j],
                                upper = x$tR[j])$value)
  }
  BIC_W <- (-2) * sum(Wloglik) + 2 * log(nobs)

  #----------------------------------------------------------------------------#
  #                          Model selection                                   #
  #----------------------------------------------------------------------------#

  BICs <- c(BIC_LPS, BIC_LN, BIC_G, BIC_W)
  modselect <- which(BICs == min(BICs))

  fdens <- function(t){
    if(t < Bl){
      val <- 0
    } else {
      if(modselect == 2){ # Log-Normal
        val <- stats::dlnorm(t, meanlog = meanlogMoM, sdlog = sdlogMoM)
      } else if (modselect == 3){ # Gamma
        val <- stats::dgamma(t, shape = shapeGMoM, rate = rateMoM)
      } else if (modselect == 4){ # Weibull
        val <- stats::dweibull(t, shape = shapeWMoM, scale = scaleMoM)
      } else if (modselect == 1){ # Flexible non-parametric
        if(t > Br){
          stop("Chosen time value not in range of B-spline domain.")
        } else{
          tgidx_up <- which(t <= tg)[1]
          tgidx_low <- tgidx_up - 1
          val <- 0.5 * (ftgLPS[tgidx_low] + ftgLPS[tgidx_up])
        }
      }
    }
    return(val)
  }

  #------------ Return output

  toc <- proc.time()-tic

  if(isTRUE(verbose)){
    if(modselect == 2){# Log-Normal
      cat("----------------------------------------------------------------\n")
      cat(paste0("Time elapsed: ", round(toc[3],3)," seconds.\n"))
      cat(paste0("Fitted density is Log-Normal with meanlog=",
                 round(meanlogMoM,3), " and sdlog=", round(sdlogMoM,3),".\n"))
      cat(paste0("Mean incubation period (days): ", round(LNincub[1,1],3),
          " with 95% CI: ", round(LNincub[1,4],3), "-",
          round(LNincub[1,5],3),".\n"))
      cat(paste0("95th percentile (days): ",
                 round(LNincub[21,1],3),
                 " with 95% CI: ", round(LNincub[21,4],3), "-",
                 round(LNincub[21,5],3),".\n"))
      cat("----------------------------------------------------------------\n")
    } else if (modselect == 3){# Gamma
      cat("----------------------------------------------------------------\n")
      cat(paste0("Time elapsed: ", round(toc[3],3)," seconds.\n"))
      cat(paste0("Fitted density is Gamma with shape=",
                 round(shapeGMoM,3), " and rate=", round(rateMoM,3),".\n"))
      cat(paste0("Mean incubation period (days): ", round(Gincub[1,1],3),
                 " with 95% CI: ", round(Gincub[1,4],3), "-",
                 round(Gincub[1,5],3),".\n"))
      cat(paste0("95th percentile (days): ",
                 round(Gincub[21,1],3),
                 " with 95% CI: ", round(Gincub[21,4],3), "-",
                 round(Gincub[21,5],3),".\n"))
      cat("----------------------------------------------------------------\n")
    } else if (modselect == 4){# Weibull
      cat("----------------------------------------------------------------\n")
      cat(paste0("Time elapsed: ", round(toc[3],3)," seconds.\n"))
      cat(paste0("Fitted density is Weibull with shape=",
                 round(shapeWMoM,3), " and scale=", round(scaleMoM,3),".\n"))
      cat(paste0("Mean incubation period (days): ", round(Wincub[1,1],3),
                 " with 95% CI: ", round(Wincub[1,4],3), "-",
                 round(Wincub[1,5],3),".\n"))
      cat(paste0("95th percentile (days): ",
                 round(Wincub[21,1],3),
                 " with 95% CI: ", round(Wincub[21,4],3), "-",
                 round(Wincub[21,5],3),".\n"))
      cat("----------------------------------------------------------------\n")
    }else if (modselect == 1){# Flexible parametric P-spline model
      cat("----------------------------------------------------------------\n")
      cat(paste0("Time elapsed: ", round(toc[3],3)," seconds.\n"))
      cat(paste0("Fitted density with flexible parametric approach.\n"))
      cat(paste0("Mean incubation period (days): ", round(LPSincub[1,1],3),
                 " with 95% CI: ", round(LPSincub[1,4],3), "-",
                 round(LPSincub[1,5],3),".\n"))
      cat(paste0("95th percentile (days): ",
                 round(LPSincub[21,1],3),
                 " with 95% CI: ", round(LPSincub[21,4],3), "-",
                 round(LPSincub[21,5],3),".\n"))
      cat("----------------------------------------------------------------\n")
    }
  }

  if(modselect == 1){ # LPS model
    outlist <- list(tg = tg, ftg = ftgLPS, stats = LPSincub,
                    mcmcrate = MCMC$accept_rate, timing = toc,
                    modselect = modselect,
                    pengpos = maxplamb$penpos,
                    penval = lambmax,
                    fdens = fdens,
                    ftgMCMC = ftgMCMC_LPS,
                    incubdat = x)
  } else if(modselect == 2){ # Log-Normal
    outlist <- list(tg = tg, ftg = ftgLN, stats = LNincub,
                    mcmcrate = MCMC$accept_rate, timing = toc,
                    modselect = modselect,
                    pengpos = maxplamb$penpos,
                    penval = lambmax,
                    fdens = fdens,
                    meanlog = meanlogMoM,
                    sdlog = sdlogMoM,
                    ftgMCMC = ftgMCMC_LN,
                    incubdat = x)
  } else if(modselect == 3){ # Gamma
    outlist <- list(tg = tg, ftg = ftgG, stats = Gincub,
                    mcmcrate = MCMC$accept_rate, timing = toc,
                    modselect = modselect,
                    pengpos = maxplamb$penpos,
                    penval = lambmax,
                    fdens = fdens,
                    shape = shapeGMoM,
                    rate = rateMoM,
                    ftgMCMC = ftgMCMC_G,
                    incubdat = x)
  } else if(modselect == 4){# Weibull
    outlist <- list(tg = tg, ftg = ftgW, stats = Wincub,
                    mcmcrate = MCMC$accept_rate, timing = toc,
                    modselect = modselect,
                    pengpos = maxplamb$penpos,
                    penval = lambmax,
                    fdens = fdens,
                    shape = shapeWMoM,
                    scale = scaleMoM,
                    ftgMCMC = ftgMCMC_W,
                    incubdat = x)
  }

  attr(outlist, "class") <- "incubestim"
  outlist
}
