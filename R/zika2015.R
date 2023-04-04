#' Data on the 2015 Zika virus disease in Colombia
#'
#' @docType data
#'
#' @description A list containing incidence data for the 2015 Zika
#' disease in Girardot (Colombia) from October 2015 to January 2016,
#' a vector of dates and a discrete serial interval distribution.
#'
#' @usage data(zika2015)
#'
#' @format A list with three components:
#' \describe{
#'  \item{\code{incidence}}{An incidence time series of length 93.}
#'  \item{\code{dates}}{A vector of dates in "YYYY-MM-DD" format.}
#'  \item{\code{si}}{A vector of probabilities corresponding to the serial
#'    interval distribution with a mean of 7 days and standard deviation
#'    of 1.5 days.}
#' }
#'
#' @source \url{https://cran.r-project.org/package=outbreaks}
#'
#' @references Rojas, D. P. et al. (2016). The epidemiology and
#'  transmissibility of Zika virus in Girardot and San Andres island,
#'  Colombia, September 2015 to January 2016. \emph{Eurosurveillance},
#'  \strong{21}(28), 30283.
#'
#'
"zika2015"
