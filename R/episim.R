#' Simulation of epidemic data based on Poisson counts
#'
#' @description
#' Based on a serial interval and a functional input for the reproduction number
#' over T days, the routine generates a set of incidence counts following a
#' Poisson model. The link between the reproduction number and the generated
#' incidence data is governed by the renewal equation. The baseline mean number
#' of cases at day 1 is fixed at 15. The mean number of cases for the remaining
#' days of the epidemic are generated following equation (2) of Azmon et al.
#' (2013).
#'
#' @usage episim(serial_interval, endepi = 50, Rpattern = 1, verbose = FALSE, plotsim = FALSE)
#'
#' @param serial_interval A vector of values for the serial interval.
#' @param endepi The total number of days of the epidemic.
#' @param Rpattern Different scenarios for the true underlying curve of
#'  Rt. Four scenarios are possible with 1,2,3,4.
#' @param verbose Should metadata on simulated edpidemic be printed?
#' @param plotsim Create a plot of the incidence time series, the true
#'  reproduction number curve and the serial interval.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Azmon, A., Faes, C., Hens, N. (2014). On the estimation of the
#' reproduction number based on misreported epidemic data. \emph{Statistics in
#' medicine}, \strong{33}(7):1176-1192. \url{https://doi.org/10.1002/sim.6015}
#'
#' @examples
#' si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
#' epidemic <- episim(serial_interval = si, Rpattern = 1)
#'
#' @export

episim <- function(serial_interval, endepi = 50, Rpattern = 1,
                   verbose = FALSE, plotsim = FALSE) {

  #-- Scenarios for the true R(t)
  if (Rpattern == 1) Rtrue <- function(t) {2.5} else if
     (Rpattern == 2) Rtrue <-function(t){if (t < 20) {res <- 2.5}
                                                else if (t >= 20) {res <- 0.7}
                                                return(res)}     else if
     (Rpattern == 3) Rtrue <- function(t){0.4 + exp(cos(t / 5))} else if
     (Rpattern == 4) Rtrue <- function(t){exp(cos(t / 15))}

  smax <- length(serial_interval)      # Length of serial_interval
  mu_y <- c()                          # Mean of y
  y <- c()                             # Incidence count

  for (t in 1:endepi) {
    if (t == 1) {
      mu_y[t] <- 15 # Mean of incidence count at day = 1.
      y[t] <- 0
      while (y[t] <= 10) { # Force a minimum number of cases at day 1
        y[t] <- stats::rpois(n = 1, lambda = mu_y[t])
      }
    } else if (t >= 2 && t <= smax) {
      mu_y[t] <- Rtrue(t) *
        sum(rev(y[1:(smax - 1)][1:(t - 1)]) *
              serial_interval[1:(smax - 1)][1:(t - 1)])
      y[t] <- stats::rpois(n = 1, lambda = mu_y[t])
    } else if (t > smax && t <= endepi) {
      mu_y[t] <-
        Rtrue(t) * sum(rev(y[(t - smax):(t - 1)]) * serial_interval)
      y[t] <- stats::rpois(n = 1, lambda = mu_y[t])
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
    }
    cat("Chosen scenario: ", Rpattern, " '", pattern, "'",  ".\n", sep = "")
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




