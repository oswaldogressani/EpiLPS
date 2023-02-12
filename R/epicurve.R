#' Plot the epidemic curve
#'
#' @description
#' This routine gives a graphical representation of the epidemic curve based
#' on case incidence data.
#'
#' @param incidence A vector for the time series of incidence.
#' @param dates A vector of dates in format "YY-MM-DD".
#' @param datelab The spacing for ticks on the x-axis. Default "7d".
#' @param col The color of the epidemic curve.
#' @param barwidth The width of the bars. Default is 1.
#' @param title Title of the plot.
#' @param xtickangle The angle of the x-ticks. Default is 0 (horizontal).
#'
#' @return A plot of the epidemic curve based on the provided time series of
#' counts.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @examples
#' si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
#' epidemic <- episim(serial_interval = si, Rpattern = 4)
#' epicurve(epidemic$y)
#' @export

epicurve <- function(incidence, dates = NULL, datelab = "7d",
                     col = "deepskyblue4", barwidth = 1,
                     title = "Epidemic curve", xtickangle = 0){
  numdays <- length(incidence)
  xlabtype <- ggplot2::geom_blank()
  Time  <- seq_len(numdays)
  Incidence <- NULL
  epiData <- data.frame(Time = Time, Incidence = incidence)
  if (!is.null(dates)){
    epiData$Time <- dates
    xlabtype <- KerDateTicker(x = dates, spacing = datelab)
  }
  plotout <- ggplot2::ggplot(data = epiData, ggplot2::aes(x = Time)) +
    ggplot2::geom_bar(stat = "identity", width = barwidth,
                      ggplot2::aes(y = Incidence), col = col, fill = col) +
    ggplot2::ggtitle(title) + xlabtype +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = xtickangle,
                                          vjust = 0.55))
  return(plotout)
}
