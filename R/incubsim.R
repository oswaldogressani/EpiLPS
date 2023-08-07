#' Simulation of incubation times
#'
#' @description
#' This routine simulates incubation times for a given generation times density
#' and incubation times density.
#'
#' @usage incubsim(incubdist = c("LogNormal","Weibull","MixWeibull", "Gamma"),
#' coarseness = 1, n = 100, tmax  = 20, tgridlen = 500, plotsim = FALSE)
#'
#' @param incubdist The distribution of the incubation period.
#' @param coarseness The average coarseness of the data.
#' @param n The number of observations.
#' @param tmax The upper bound on which to evaluate the \code{incubdist} density.
#' @param plotsim Create a plot of the simulated data.
#' @param tgridlen The number of grid points on which to evaluate the density.
#'
#' @return An object of class \code{incubsim} consisting of....
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Ferretti, L., Wymant, C., Kendall, et al. (2020). Quantifying
#' SARS-CoV-2 transmission suggests epidemic control with digital contact
#' tracing. \emph{Science}, \strong{368}(6491), eabb6936.
#' @references Backer, J. A., Klinkenberg, D., & Wallinga, J. (2020). Incubation
#' period of 2019 novel coronavirus (2019-nCoV) infections among
#' travellers from Wuhan, China, 20–28 January 2020. \emph{Eurosurveillance},
#' \strong{25}(5), 2000062.
#' @references Donnelly, C. A., Ghani, A. C., Leung, G. M., et al. (2003).
#' Epidemiological determinants of spread of causal agent of severe acute
#' respiratory syndrome in Hong Kong. \emph{The Lancet}, \strong{361}(9371),
#' 1761-1766.
#'
#' @examples
#' incubsim()
#'
#' @export


incubsim <- function(incubdist = c("LogNormal","Weibull","MixWeibull", "Gamma"),
                     coarseness = 1, n = 100, tmax = 20, tgridlen = 500, plotsim = FALSE){

  mixWeibull <- function(n, alpha, Ishape1 = 2.356, Iscale1 = 3.564,
                         Ishape2 = 9.351, Iscale2 = 12.563) {
    #Probability weights in the mixture
    w1 <- alpha
    w2 <- 1 - alpha

    #Sample scheme with replacement and weights on the elements being sampled
    w <- sample(c(1, 2), n, replace = TRUE, prob = c(w1, w2))

    #Sample of size n from the mixture of 2 Weibulls
    sample <- (w == 1) * stats::rweibull(n = n, shape = Ishape1, scale = Iscale1) +
      (w == 2) * stats::rweibull(n = n, shape = Ishape2, scale = Iscale2)

    # set.seed(123)
    # mixsample <- mixWeibull(n = 50000000, alpha = 0.5)
    # summary(mixsample$sample)
    # sd(mixsample$sample)
    # quantile(mixsample$sample, probs = c(0.05,0.25,0.5,0.75,0.95))

    outlist <- list(sample = sample, Ishape1 = Ishape1, Iscale1 = Iscale1,
                    Ishape2 = Ishape2, Iscale2 = Iscale2)
    return(outlist)
  }

  # Domain of Generation time/Incubation time density
  tt <- seq(0, tmax, length = tgridlen)

  # Generation interval Weibull parameters #Ferretti et al. 2020
  GIshape <- 2.826
  # GIshape <- 3.4395
  GIscale <- 5.665

  if (match.arg(incubdist) == "LogNormal") { # Ferretti et al. 2020
    Incubmean <- 1.644
    Incubsd <- 0.363
    if (coarseness == 1) {
      dirforce <- 0.185
    } else if (coarseness == 2) {
      dirforce <- 0.375
    }
    Ipdf <- stats::dlnorm(tt, meanlog = Incubmean, sdlog = Incubsd)
    meanI <- exp(Incubmean + 0.5 * (Incubsd ^ 2))
    medianI <- exp(Incubmean)
    sdI <- sqrt(exp(2 * Incubmean + Incubsd ^ 2) * (exp(Incubsd ^ 2) - 1))
    Iparams <- c(Incubmean, Incubsd)
  } else if (match.arg(incubdist) == "Weibull") { # Backer et al. 2020
    Ishape <- 3.00060
    Iscale <- 7.170345
    if (coarseness == 1) {
      dirforce <- 0.160
    } else if (coarseness == 2) {
      dirforce <- 0.320
    }
    Ipdf <- stats::dweibull(tt, shape = Ishape, scale = Iscale)
    meanI <- Iscale * gamma(1 + 1 / Ishape)
    medianI <- stats::qweibull(p = 0.5, shape = Ishape, scale = Iscale)
    sdI <- sqrt((Iscale ^ 2) * ((gamma(1 + 2 / Ishape) - (gamma(1 + 1 / Ishape) ^ 2))))
    Iparams <- c(Ishape, Iscale)
  } else if (match.arg(incubdist) == "MixWeibull"){ # Artificial
    if (coarseness == 1) {
      dirforce <- 0.132
    } else if (coarseness == 2) {
      dirforce <- 0.264
    }
    mmix <- mixWeibull(n=1, alpha = 0.5)
    Ishape1 <- mmix$Ishape1
    Iscale1 <- mmix$Iscale1
    Ishape2 <- mmix$Ishape2
    Iscale2 <- mmix$Iscale2

    Ipdf <- 0.5 * stats::dweibull(tt, shape = Ishape1, scale = Iscale1) +
      0.5 * stats::dweibull(tt, shape = Ishape2, scale = Iscale2)
    meanmix <- 0.5 * (Iscale1 * gamma(1 + 1 / Ishape1) +
                      Iscale2 * gamma(1 + 1 / Ishape2))
    sdmix <- 4.622
    qmix <- c(1.371, 1.885,2.301,2.680,3.050,3.434,3.856,4.361,5.075,7.191,
              9.876,10.701,11.251,11.692,12.080,12.446,12.814,13.218,13.734)
  } else if(match.arg(incubdist) == "Gamma") { # Donnelly 2003
    Ishape <- 1.739646
    Irate <- 0.4566
    if (coarseness == 1) {
      dirforce <- 0.262
    } else if (coarseness == 2) {
      dirforce <- 0.558
    }
    Ipdf <- stats::dgamma(tt, shape = Ishape, rate = Irate)
    meanI <- Ishape / Irate
    medianI <- stats::qgamma(p = 0.5, shape = Ishape, rate = Irate)
    sdI <- sqrt(Ishape / (Irate ^ 2))
    Iparams <- c(Ishape, Irate)
  }

  #-- Generation of symptom onset data and infection exposure intervals
  GI  <-  c()  # Vector to host generation times
  SI  <-  c()  # Vector to host serial interval times
  TI1 <-  c()  # Vector to host infection times of source/infector
  TI2 <-  c()  # Vector to host infection times of recipient/infectee
  I1  <-  c()  # Vector to host incubation period for source/infector
  I2  <-  c()  # Vector to host incubation period for recipient/infectee
  SO1 <-  c()  # Vector to host symptom onset times of source/infector
  SO2 <-  c()  # Vector to host symptom onset times of recipient/infectee
  EL1 <-  c()  # Vector to host lower bound of exposure time of source/infector
  ER1 <-  c()  # Vector to host upper bound of exposure time of source/infector
  EL2 <-  c()  # Vector to host lower bound of exposure time of recipient/infectee
  ER2 <-  c()  # Vector to host upper bound of exposure time of recipient/infectee

  i <- 1

  while(i <= n){
    # Step 1: Draw generation time from its pdf (Weibull)
    GI[i] <- stats::rweibull(n = 1, shape = GIshape, scale = GIscale)

    # Step 2: Set time of infection of infector/source (0,365)
    TI1[i] <- stats::runif(1, min = 0, max = 365)

    # Step 3: Compute the time of infection of infectee/recipient
    TI2[i] <- TI1[i] + GI[i]

    # Step 4: Generate a pair of incubation periods from pdf (LogNormal)
    if (match.arg(incubdist) == "LogNormal") {
      I1[i] <- stats::rlnorm(1, meanlog = Incubmean, sdlog = Incubsd)
      I2[i] <- stats::rlnorm(1, meanlog = Incubmean, sdlog = Incubsd)
    } else if (match.arg(incubdist) == "Weibull") {
      I1[i] <- stats::rweibull(1, shape = Ishape, scale = Iscale)
      I2[i] <- stats::rweibull(1, shape = Ishape, scale = Iscale)
    } else if (match.arg(incubdist) == "MixWeibull") {
      I1[i] <- mixWeibull(n = 1, alpha = 0.5)$sample
      I2[i] <- mixWeibull(n = 1, alpha = 0.5)$sample
    } else if (match.arg(incubdist) == "Gamma") {
      I1[i] <- stats::rgamma(1, shape = Ishape, rate = Irate)
      I2[i] <- stats::rgamma(1, shape = Ishape, rate = Irate)
    }

    # Step 5: Compute times of symptom onset
    SO1[i] <- TI1[i] + I1[i]
    SO2[i] <- TI2[i] + I2[i]
    SI[i] <- SO2[i] - SO1[i]

    # Step 6: Compute lower/upper bounds of exposure times
    minexpo1 <- TI1[i] - dirforce * I1[i]
    minexpo2 <- TI2[i] - dirforce * I2[i]

    exposure_draw <- function(tlow,tup){
      norm_mean <- 0.5 * (tlow + tup)
      sdnew <- 0
      probtar <- 0.999
      ptol <- stats::pnorm(q = tup, mean = norm_mean, sd = sdnew)
      while(ptol >= probtar){
        sdnew <- sdnew + 0.001
        ptol <- stats::pnorm(q = tup, mean = norm_mean, sd = sdnew)
      }
      tdraw <- 1e+10
      while(tdraw >= tup || tdraw <= tlow) {
        tdraw <- stats::rnorm(n = 1, mean = norm_mean, sd = sdnew)
      }
      return(tdraw)
    }

    # Draw left bound of exposure for infector/infectee
    EL1_temp <- exposure_draw(tlow = minexpo1, tup = TI1[i])
    EL2_temp <- exposure_draw(tlow = minexpo2, tup = TI2[i])

    while(EL2_temp < EL1_temp){# Ferretti et al. 2020 constraints
      EL1_temp <- exposure_draw(tlow = minexpo1, tup = TI1[i])
      EL2_temp <- exposure_draw(tlow = minexpo2, tup = TI2[i])
    }
    EL1[i] <- EL1_temp
    EL2[i] <- EL2_temp

    maxexpo1 <- TI1[i] + dirforce * I1[i]
    maxexpo2 <- TI2[i] + dirforce * I2[i]

    # Draw right bound of exposure for infector/infectee
    ER1_temp <- exposure_draw(tlow = TI1[i], tup = maxexpo1)
    ER2_temp <- exposure_draw(tlow = TI2[i], tup = maxexpo2)

    while(ER1_temp > min(SO1[i],ER2_temp)){# Ferretti et al. 2020 constraints
      ER1_temp <- exposure_draw(tlow = TI1[i], tup = maxexpo1)
      ER2_temp <- exposure_draw(tlow = TI2[i], tup = maxexpo2)
    }

    ER1[i] <- ER1_temp
    ER2[i] <- ER2_temp

    ediff1 <- ER1_temp - EL1_temp
    ediff2 <- ER2_temp - EL2_temp
    if(ediff1 > 7 || ediff2 > 7){
      next
    }
    i <- i + 1
  }

  # Gather data in matrix (observable)
  Dobs <- matrix(0, nrow = n, ncol = 9)
  colnames(Dobs) <- c("PairIndex", "SympOnset1", "SympOnset2", "SI",
                      "ExpoL1", "ExpoR1", "ExpoL2", "ExpoR2", "GI")
  Dobs[,1] <- seq_len(n)
  Dobs[,2] <- SO1
  Dobs[,3] <- SO2
  Dobs[,4] <- SO2-SO1
  Dobs[,5] <- EL1
  Dobs[,6] <- ER1
  Dobs[,7] <- EL2
  Dobs[,8] <- ER2
  Dobs[,9] <- GI
  Dobs <- data.frame(Dobs)

  # Width of infection exposure windows
  expoL <- c(Dobs$ExpoL1, Dobs$ExpoL2)
  expoR <- c(Dobs$ExpoR1, Dobs$ExpoR2)
  expo_diff <- expoR-expoL
  exposure_mean <- mean(expo_diff)

  Dobs_incub <- matrix(0, nrow = (2*n), ncol = 2)
  colnames(Dobs_incub) <- c("tL","tR")
  Dobs_incub[,1] <- c(Dobs$SympOnset1-Dobs$ExpoR1,
                      Dobs$SympOnset2-Dobs$ExpoR2)
  Dobs_incub[,2] <- c(Dobs$SympOnset1-Dobs$ExpoL1,
                      Dobs$SympOnset2-Dobs$ExpoL2)
  Dobs_incub <- Dobs_incub[sample(seq_len(2*n), replace = FALSE)[1:n],]
  Dobs_incub <- as.data.frame(Dobs_incub)

  # Check constraints on exposure intervals of Ferretti et al. (2020)
  C1 <- all(ER1 <= apply(cbind(SO1,ER2), 1, min))
  C2 <- all(EL2 >= EL1)
  C3 <- all(SO2 >= ER2)
  Ferretti_constr <- all(c(C1,C2,C3))

  #-- Plot of generation interval and incubation time pdf

  # Generation interval distribution
  GTpdf <- stats::dweibull(tt, shape = GIshape, scale = GIscale)
  meanGT <- GIscale * gamma(1 + 1 / GIshape)
  sdGT <- sqrt((GIscale ^ 2) * (gamma(1 + 2 / GIshape) -
                                  (gamma(1 + 1 / GIshape) ^ 2)))

  # Plot
  if(plotsim == TRUE){
    graphics::par(mfrow = c(2,2))
    plot(tt, GTpdf, type = "l", col = "black", lwd = 2,
         xlab = "Time (days)", ylab = "Density", ylim = c(0,0.25))
    graphics::rug(GI)
    graphics::lines(tt, Ipdf, type = "l", col = "red", lwd = 2)
    graphics::rug(I1, col = "red", side = 3)
    graphics::legend("topright", c("Generation interval", "Incubation period"),
           lwd = c(2,2), lty = c(1,1), col = c("black","red"), bty = "n")

    graphics::hist(Dobs$SI, breaks = sqrt(n) + 5, freq = FALSE,
                   xlab = "Time (days)",
         ylab = "Serial interval distribution", main = "")
    graphics::lines(stats::density(Dobs$SI), type = "l",
                    col = "darkgreen", lwd = 2)

    graphics::hist(Dobs$GI, breaks = sqrt(n) + 5, freq = FALSE,
                   xlab = "Time (days)",
         ylab = "Generation interval distribution", main = "")
    graphics::lines(stats::density(Dobs$GI), type = "l",
                    col = "purple", lwd = 2)

    graphics::hist(expo_diff, xlab = "Time (days)", ylab = "",
         main ="Duration of exposure", freq = TRUE, col = "cornflowerblue")
  }
  graphics::par(mfrow = c(1,1))

  if (match.arg(incubdist) == "MixWeibull"){
    outlist <- list(Dobs = Dobs,
                    Ferrconst = Ferretti_constr,
                    Constr = paste0("Ferretti constraints passed: ",
                                    Ferretti_constr),
                    GIdist = paste0("Generation interval with mean: ",
                                    round(meanGT,2), " days and standard deviation of: ",
                                    round(sdGT,2), " days."),
                    Incubdist = paste0("Incubation distribution is a Weibull mixture"),
                    Expoduration = paste0("Mean duration of exposure window: ",
                                          round(exposure_mean,2), " days."),
                    Expomean = exposure_mean,
                    Dobsincub = Dobs_incub,
                    coarseness = coarseness,
                    meanmix = meanmix,
                    sdmix = sdmix,
                    qmix = qmix,
                    Ipdf = Ipdf,
                    tdom = tt)
  } else {
    outlist <- list(Dobs = Dobs,
                  Ferrconst = Ferretti_constr,
                  Constr = paste0("Ferretti constraints passed: ",
                                  Ferretti_constr),
                  GIdist = paste0("Generation interval with mean: ",
                                  round(meanGT,2), " days and standard deviation of: ",
                                  round(sdGT,2), " days."),
                  Incubdist = paste0("Incubation distribution is ",
                                     match.arg(incubdist),"(",Iparams[1],Iparams[2],")",
                                     " with mean ",
                                     round(meanI,2), " days and standard deviation of: ",
                                     round(sdI,2), " days."),
                  Iparams = Iparams,
                  Expoduration = paste0("Mean duration of exposure window: ",
                                        round(exposure_mean,2), " days."),
                  Expomean = exposure_mean,
                  Dobsincub = Dobs_incub,
                  coarseness = coarseness,
                  tdom = tt)
  }

  return(outlist)
}
