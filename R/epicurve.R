#' Plot the epidemic curve
#'
#' @description
#' This routine gives a graphical representation of the epidemic curve based
#' on incidence data.
#'
#' @usage epicurve(epicounts, dates = NULL,
#' datelab = c("7d", "1m", "3m", "6m"), seriesname = "Series 1",
#' themetype = c("gray","classic","light","dark"), colbar = "deepskyblue4",
#' title = "Epidemic curve", poslegend = "right", dirlegend = "vertical",
#' xlabel = "Time (days)", ylabel = "Incidence")
#'
#' @param epicounts A time series of counts.
#' @param dates A vector of dates in format "YY-MM-DD".
#' @param datelab The spacing for ticks on the x-axis. Either 7 days, 1 month,
#'  3 months or 6 months.
#' @param seriesname The legend name of the time series.
#' @param themetype Type of theme for the plot.
#' @param colbar The color of the bars.
#' @param title Title of the plot.
#' @param poslegend The position of the legend in the plot.
#' @param dirlegend The direction of the legend in the plot.
#' @param xlabel The label on the x-axis.
#' @param ylabel The label on the y-axis.
#'
#' @return A plot of the epidemic curve based on the provided time series of
#' counts.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @examples
#' si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
#' epidemic <- episim(serial_interval = si, Rpattern = 4)
#' epi_inci <- as.data.frame(epidemic$y)
#' epicurve(epi_inci)
#' @export

epicurve <- function(epicounts, dates = NULL,
                     datelab = c("7d", "1m", "3m", "6m"),
                     seriesname = "Series 1",
                     themetype = c("gray","classic","light","dark"),
                     colbar = "deepskyblue4",
                     title = "Epidemic curve",
                     poslegend = "right",
                     dirlegend = "vertical",
                     xlabel = "Time (days)",
                     ylabel = "Incidence"){

  if(!is.data.frame(epicounts)) stop("Count data should be a data frame")
  numSeries <- 1
  totDays <- nrow(epicounts)
  if (!is.null(dates)) {
    epicounts$Time <- dates
  } else{
    epicounts$Time <- seq_len(totDays)
  }
  colnames(epicounts) <- c(seriesname, "Time")
  Time <- NULL
  Incidence <- NULL
  Series <- NULL
  epiData <- data.frame(epicounts$Time, rep(seriesname, totDays),
                        epicounts[, 1])
  colnames(epiData) <- c("Time","Series","Incidence")

  # Adapt x-axis labels if dates are provided
  if(!is.null(dates)) {
    datelab <- match.arg(datelab)
    if (datelab == "7d") {
      xlabtype <- eval(parse(text = "ggplot2::scale_x_date(date_breaks = '7 days')"))
    } else if (datelab == "1m") {
      xlabtype <- eval(parse(text = "ggplot2::scale_x_date(date_breaks = '1 month')"))
    } else if (datelab == "3m") {
      xlabtype <- eval(parse(text = "ggplot2::scale_x_date(date_breaks = '3 months')"))
    } else if (datelab == "6m") {
      xlabtype <- eval(parse(text = "ggplot2::scale_x_date(date_breaks = '6 months')"))
    }
  } else{
    xlabtype <- NULL
  }

  # Theme choice
  themetype <- match.arg(themetype)
  if (themetype == "classic") {
    themeval <- eval(parse(text = "ggplot2::theme_classic()"))
  } else if (themetype == "gray") {
    themeval <- eval(parse(text = "ggplot2::theme_gray()"))
  } else if (themetype == "light") {
    themeval <- eval(parse(text = "ggplot2::theme_light()"))
  } else if (themetype == "dark") {
    themeval <- eval(parse(text = "ggplot2::theme_dark()"))
  }

  colpalette <- colbar
  epicurvePlot <-
    ggplot2::ggplot(data = epiData, ggplot2::aes(x = Time))+
    ggplot2::geom_bar(stat = "identity", width = 1,
                      ggplot2::aes(y = Incidence, fill = Series),
                      position = ggplot2::position_dodge(1)) +
    ggplot2::scale_fill_manual(values=colpalette[1:numSeries])+
    ggplot2::labs(fill="Daily incidence")+
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(title) +
    themeval +
    xlabtype +
    ggplot2::theme(legend.position = poslegend,
                   legend.direction = dirlegend,
                   legend.title = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 13),
                   plot.title = ggplot2::element_text(size = 17),
                   axis.title.x = ggplot2::element_text(size = 14),
                   axis.title.y = ggplot2::element_text(size = 14),
                   axis.text.x = ggplot2::element_text(size = 14),
                   axis.text.y = ggplot2::element_text(size = 14))

  return(epicurvePlot)
}
