#' Simulation of incidence count data
#'
#' @description
#' Based on a serial interval and a functional input for the reproduction number
#' over T days, the routine generates a set of incidence counts following a
#' Poisson or negative binomial model. The link between the reproduction number
#' and the generated incidence data is governed by the renewal equation. The
#' baseline (mean) number of cases at day 1 is fixed at 10. The mean number of
#' cases for the remaining days of the epidemic are generated following
#' equation (2) of Azmon et al. (2013).
#'
#' @usage episim(serial_interval, endepi = 50, Rpattern = 1, Rconst = 2.5,
#'        dist = c("poiss", "negbin"), overdisp = 1, verbose = FALSE, plotsim = FALSE)
#'
#' @param serial_interval A vector of values for the discrete serial interval
#'  (must sum to 1).
#' @param endepi The total number of days of the epidemic.
#' @param Rpattern Different scenarios for the true underlying curve of
#'  Rt. Six scenarios are possible with 1,2,3,4,5,6.
#' @param Rconst The constant value of R (if scenario 1 is selected), default is 2.5.
#' @param dist The distribution from which to sample the incidence of cases.
#'  Either Poisson (default) or negative binomial.
#' @param overdisp Overdispersion parameter for the negative binomial setting.
#' @param verbose Should metadata on simulated epidemic be printed?
#' @param plotsim Create a plot of the incidence time series, the true
#'  reproduction number curve and the serial interval.
#'
#' @return An object of class \code{episim} consisting of a list with the
#'  generated time series of cases, the mean vector of the Poisson/negative binomial
#'  distribution, the true underlying R function for the data generating process and the
#'  chosen serial interval distribution.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Azmon, A., Faes, C., Hens, N. (2014). On the estimation of the
#' reproduction number based on misreported epidemic data. \emph{Statistics in
#' medicine}, \strong{33}(7):1176-1192.
#'
#' @examples
#' si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
#' epidemic <- episim(serial_interval = si, Rpattern = 1)
#'
#' @export

episim <- function(serial_interval, endepi = 50, Rpattern = 1, Rconst = 2.5,
                   dist = c("poiss", "negbin"), overdisp = 1,
                   verbose = FALSE, plotsim = FALSE) {

  #-- Scenarios for the true R(t)
  if (Rpattern == 1) Rtrue <- function(t) {Rconst} else if
     (Rpattern == 2) Rtrue <- function(t){if (t < 20) {res <- 2}
                                                else if (t >= 20) {res <- 0.9}
                                                return(res)}     else if
     (Rpattern == 3) Rtrue <- function(t){0.25 + exp(cos(t / 7))} else if
     (Rpattern == 4) Rtrue <- function(t){exp(cos(t / 15))} else if
     (Rpattern == 5) Rtrue <- function(t){0.5 * (exp(sin(pi * t / 9)) +
                                              1.5 * exp(cos(t / 4)))} else if
     (Rpattern == 6) Rtrue <- function(t){1.5 + cos(0.8 * pi * t / 15) +
                                          0.5 * t ^ 2 / 400
     }

  smax <- length(serial_interval)      # Length of serial_interval
  mu_y <- c()                          # Mean of y
  y <- c()                             # Incidence count

  sampling_dist <- match.arg(dist)     # Chosen distribution to sample from

  for (t in 1:endepi) {
    if (t == 1) {
      mu_y[t] <- 10 # Mean of incidence count at day = 1.
      y[t] <- 10
    } else if (t >= 2 && t <= smax) {
      mu_y[t] <- Rtrue(t) *
        sum(rev(y[1:(smax - 1)][1:(t - 1)]) *
              serial_interval[1:(smax - 1)][1:(t - 1)])
      if (sampling_dist == "poiss") {
        y[t] <- stats::rpois(n = 1, lambda = mu_y[t])
      } else if (sampling_dist == "negbin") {
        y[t] <- stats::rnbinom(n = 1, mu = mu_y[t], size = overdisp)
      }
    } else if (t > smax && t <= endepi) {
      mu_y[t] <-
        Rtrue(t) * sum(rev(y[(t - smax):(t - 1)]) * serial_interval)
      if (sampling_dist == "poiss") {
        y[t] <- stats::rpois(n = 1, lambda = mu_y[t])
      } else if (sampling_dist == "negbin") {
        y[t] <- stats::rnbinom(n = 1, mu = mu_y[t], size = overdisp)
      }
    }
  }

  #-- Print information on generated data
  if (verbose == TRUE) {
    if (Rpattern == 1) {
      pattern = "Constant Rt"
    } else if (Rpattern == 2) {
      pattern = "Step function Rt (intervention mimick)"
    } else if (Rpattern == 3) {
      pattern = "Wiggly Rt"
    } else if (Rpattern == 4) {
      pattern = "Decaying Rt"
    } else if (Rpattern == 5) {
      pattern = "Wiggly then stable Rt"
    } else if (Rpattern == 6) {
      pattern = "Increasing Rt"
    }
    cat("Chosen scenario: ", Rpattern, " '", pattern, "'",  ".\n", sep = "")
    if (sampling_dist == "poiss") {
      cat("Incidence of cases generated from a Poisson distribution. \n")
    } else if (sampling_dist == "negbin"){
   cat("Incidence of cases generated from a negative binomial distribution. \n")
   cat("Overdispersion parameter value: ", overdisp, ".\n", sep = "")
    }
    cat("Total number of days of epidemic: ", endepi, ".\n", sep = "")
  }

  #-- Plot simulation-based results
  if (plotsim == TRUE) {
    # Plot 1 (incidence data)
    daysvec <- seq_len(endepi)
    incidata <- data.frame(daysvec = daysvec, y = y)

    plotinci <- ggplot2::ggplot(data = incidata,
                                ggplot2::aes(x = daysvec, y = y)) +
      ggplot2::geom_bar(stat = "identity", width = 0.35, color = "steelblue",
                          fill = "steelblue") +
      ggplot2::xlab("Days") +
      ggplot2::ylab("Incidence") +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 14),
        axis.title.y = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14)
      )

    # Plot 2 (serial interval)
    silen <- seq_len(smax)
    sispec <- data.frame(silen = silen, serial_interval = serial_interval)

    plotsint <- ggplot2::ggplot(data = sispec,
                                ggplot2::aes(x = silen, y = serial_interval)) +
                ggplot2::scale_x_discrete(name = "",
                                          limits = as.character(silen)) +
                ggplot2::geom_bar(stat = "identity", width = 0.35,
                                  color = "forestgreen", fill = "forestgreen") +
      ggplot2::xlab("Serial interval index") +
      ggplot2::ylab("Serial interval distribution") +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 14),
        axis.title.y = ggplot2::element_text(size = 14),
        axis.text.x =  ggplot2::element_text(size = 14),
        axis.text.y =  ggplot2::element_text(size = 14)
      )

    # Plot 3 (true R function)
    tdom <- seq(1, endepi, length = 500)
    Rtdom <- sapply(tdom, Rtrue)
    Rfunc <- data.frame(tdom = tdom, Rtdom = Rtdom)

    plotR <- ggplot2::ggplot(data = Rfunc, ggplot2::aes(x = tdom, y = Rtdom)) +
              ggplot2::geom_line(size = 1.1) +
              ggplot2::xlab("Days") +
              ggplot2::ylab("R") +
              ggplot2::theme(
                axis.title.x = ggplot2::element_text(size = 14),
                axis.title.y = ggplot2::element_text(size = 14),
                axis.text.x = ggplot2::element_text(size = 14),
                axis.text.y = ggplot2::element_text(size = 14)
              )
    plot_simsummary <- gridExtra::grid.arrange(plotinci, plotsint, plotR,
                                               nrow = 1)
  }


  outlist <- list(y = y, mu_y = mu_y, Rtrue = Rtrue,
                  serial_interval = serial_interval)
  attr(outlist, "class") <- "episim"
  outlist
}




