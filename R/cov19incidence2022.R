#' Incidence data for Belgium in 2022
#'
#' @docType data
#'
#' @description A data frame with six columns including the calendar date of
#' COVID-19 infection, the date on which cases have been reported and the number of cases.
#'
#' @usage data(cov19incidence2022)
#'
#' @format
#' \describe{
#'  \item{\code{t}}{A numeric variable associated to a calendar date.}
#'  \item{\code{d}}{A numeric variable indicating the delay of reporting.}
#'  \item{\code{Date}}{The calendar date for incidence.}
#'  \item{\code{Rep.date}}{The calendar date for the reporting of incidence.}
#'  \item{\code{Cases}}{The number of cases.}
#'  \item{\code{Reported}}{Indicates whether cases are already reported or not.}
#' }
#'
#'@source \url{https://epistat.sciensano.be/covid/}
#'
"cov19incidence2022"
