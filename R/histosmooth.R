#' Histogram smoothing with Laplacian-P-splines
#'
#' @description
#' This function provides a smooth density estimate to a histogram using
#' Laplacian-P-splines. The B-spline basis is computed on the midpoints of
#' the histogram bins. The default number of (cubic) B-splines is 30 and
#' a third-order penalty is specified. The negative binomial distribution is
#' used to model the number of observations falling in each bin.
#'
#' @usage histosmooth(x, xl = min(x), xr = max(x), K = 30)
#'
#' @param x A vector of real numbers from which the histogram will be constructed.
#' @param xl The left bound for the domain of \code{x} which also coincides
#'  with the left bound of the domain of the B-spline basis. The default is
#'  taken to be the minimum of \code{x}.
#' @param xr The right bound for the domain of \code{x} which also coincides
#'  with the right bound of the domain of the B-spline basis. The default is
#'  taken to be the maximum of \code{x}.
#' @param K Number of B-splines in the basis.
#'
#' @return A list containing the left (\code{xl})
#' and right (\code{xr}) bounds of the domain of the estimated density, the
#' binwidth and a function to be evaluated between \code{xl} and \code{xr}.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Gressani, O. and Lambert, P. (2018). Fast Bayesian inference
#'  using Laplace approximations in a flexible promotion time cure model based
#'  on P-splines. \emph{Computational Statistical & Data Analysis} \strong{124}:
#'  151-167.
#' @references Gressani, O., Wallinga, J., Althaus, C. L., Hens, N. and Faes, C.
#'  (2022). EpiLPS: A fast and flexible Bayesian tool for estimation of the
#'  time-varying reproduction number. \emph{Plos Computational Biology},
#'  \strong{18}(10): e1010618.
#'
#' @examples
#' # Old Faithful geyser application
#' data(eruptions)
#' x <- eruptions
#' ffit <- histosmooth(x, xl = 1, xr = 6)
#' tt <- seq(ffit$xl, ffit$xr, length = 500)
#' dtt <- tt[2] - tt[1]
#' graphics::hist(x, breaks = seq(ffit$xl, ffit$xr, by = ffit$binwidth),
#'      freq = FALSE, ylim = c(0, 0.8), main = "Old Faithful Geyser",
#'     xlab = "Eruption time (minutes)")
#' densfit <- sapply(tt, ffit$fdens)
#' densfit <- densfit / (sum(densfit * dtt))
#' graphics::lines(tt, densfit, col = "red", lwd = 2)
#'
#' @export

histosmooth <- function(x, xl = min(x), xr = max(x), K = 30){

  # Compute histogram and define B-spline basis
  kerndens <- stats::density(x)
  binwidth <- kerndens$bw
  n <- length(x)
  histo <- graphics::hist(x, breaks = seq(xl, xr, by = binwidth), plot = FALSE)
  xmid  <- histo$mids
  y <- histo$counts
  B <- Rcpp_KercubicBspline(xmid, lower = xl, upper = xr, K = K)

  # Penalty matrix specification (order 3)
  D <- diag(K)
  penorder <- 3
  for(k in 1:penorder){D <- diff(D)}
  P <- t(D) %*% D           # Penalty matrix of dimension c(K,K)
  P <- P + diag(1e-06, K)   # Perturbation to ensure P is full rank

  # Prior specification
  priorspec <- Rmodelpriors()
  a_delta <- priorspec$a_delta
  b_delta <- priorspec$b_delta
  phi <- priorspec$phi
  a_rho <- priorspec$a_rho
  b_rho <- priorspec$b_rho

  logphyper <- function(x) {
    w <- x[1] # w = log(rho)
    v <- x[2] # v = log(lambda)

    # Laplace approximation
    LL <- Rcpp_KerLaplace(theta0 = rep(1.5,K), exp(w), exp(v), K,
                          KerPtheta(Dobs = y, BB = B, Pen = P)$Dlogptheta,
                          KerPtheta(Dobs = y, BB = B, Pen = P)$D2logptheta)
    thetastar <- as.numeric(LL$Lapmode)
    logdetSigstar <- Re(LL$logdetSigma)

    equal <- (-1) * (0.5 * logdetSigstar + 0.5 * (K + phi) * v + a_rho * w -
                       (0.5 * phi + a_delta) * log(0.5 * phi * exp(v) + b_delta) +
                       KerLikelihood(Dobs = y, BB = B)$loglik(thetastar, exp(w)) -
                       0.5 * exp(v) * sum((thetastar * P) %*% thetastar) -
                       b_rho * exp(w))

    return(equal)
  }

  Dlogphyper <- function(x, thetavec){
    w <- x[1] # w = log(rho)
    v <- x[2] # v = log(lambda)

    # Laplace approximation
    LL <- Rcpp_KerLaplace(theta0 = thetavec, exp(w), exp(v), K,
                          KerPtheta(Dobs = y, BB = B, Pen = P)$Dlogptheta,
                          KerPtheta(Dobs = y, BB = B, Pen = P)$D2logptheta)
    thetastar <- as.numeric(LL$Lapmode)
    Btheta <- as.numeric(B%*%thetastar)

    Dloglik_w <- sum(exp(w) * (digamma(y + exp(w)) - digamma(exp(w)) +
                                 (1 + w) - (log(exp(Btheta) + exp(w)) +
                                              (1 / (1 + exp(Btheta - w))))) -
                       y * (1 / (1 + exp(Btheta - w))))

    grad_w <- 0.5 * Re(LL$dSigw) + a_rho - b_rho * exp(w) + Dloglik_w

    grad_v <- 0.5 * Re(LL$dSigv) + 0.5 * (K + phi) - ((0.5 * phi + a_delta) /
                                                        (0.5 * phi * exp(v) + b_delta)) *
      (0.5 * phi * exp(v)) - 0.5 * exp(v) * sum((thetastar * P) %*% thetastar)

    return(list(grad=c(grad_w, grad_v), thetastar = thetastar))
  }

  optimhyper <- stats::optim(c(1, 5), fn = logphyper)
  hypermap <- optimhyper$par
  optimconverged <- (optimhyper$convergence == 0)
  disphat <- exp(hypermap[1])
  lambhat <- exp(hypermap[2])
  Lap_approx <- Rcpp_KerLaplace(theta0 = rep(1.5,K), disphat, lambhat, K,
                              KerPtheta(Dobs = y, BB = B, Pen = P)$Dlogptheta,
                              KerPtheta(Dobs = y, BB = B, Pen = P)$D2logptheta)
  thetahat <- as.numeric(Lap_approx$Lapmode)

  fdens <- function(t){
    if(t < xl | t > xr){
      fhat <- 0
    } else{
      Bt <- Rcpp_KercubicBspline(t, lower = xl, upper = xr, K = K)
      fhat <- as.numeric((1 / (n * binwidth)) * exp(Bt %*% thetahat))
    }
    return(fhat)
  }

  # Computation of BIC via effective dimension
  I <- (-1) * KerLikelihood(Dobs = y, BB = B)$D2loglik(thetahat, rho = disphat)
  ED <- sum(diag(solve(I + lambhat * P) %*% I))
  BIC <- (-2) * KerLikelihood(Dobs = y, BB = B)$loglik(thetahat,
                  rho = disphat) + ED * log(n)

  outlist <- list(xl = xl, xr = xr, binwidth = binwidth, fdens = fdens,
                  BIC = BIC)
  return(outlist)
}
