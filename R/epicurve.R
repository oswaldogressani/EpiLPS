#' Plot the epidemic curve
#'
#' @description
#' This routine gives a graphical representation of the epidemic curve based
#' on an incidence time series.
#'
#' @usage epicurve(incidence, dates = NULL, datelab = "7d", col = "deepskyblue4", barwidth = 1,
#'            title = "Epidemic curve", xtickangle = 0, smooth = NULL,  smoothcol = "orange")
#'
#' @param incidence A vector containing the incidence time series
#' @param dates A vector of dates in format "YYYY-MM-DD".
#' @param datelab The spacing for ticks on the x-axis. Default "7d".
#' @param col The color of the epidemic curve.
#' @param barwidth The width of the bars. Default is 1.
#' @param title Title of the plot.
#' @param xtickangle The angle of the x-ticks. Default is 0 (horizontal).
#' @param smooth An object of class Rt obtained with the \code{estimR} or
#'  \code{estimRmcmc} routine. It is used to draw a smoothed estimate of the
#'  epidemic curve based on the Laplacian-P-splines model.
#' @param smoothcol The color of the smoothed curve and associated credible
#' interval.
#'
#' @return A plot of the epidemic curve.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @examples
#' si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
#' epidemic <- episim(si = si, Rpattern = 4)
#' epicurve(epidemic$y)
#'
#' @export

epicurve <- function(incidence, dates = NULL, datelab = "7d",
                     col = "deepskyblue4", barwidth = 1,
                     title = "Epidemic curve", xtickangle = 0, smooth = NULL,
                     smoothcol = "orange"){
  numdays <- length(incidence)
  xlabtype <- ggplot2::geom_blank()
  Time  <- seq_len(numdays)
  Incidence <- NULL
  CIlow <- NULL
  CIup <- NULL
  epiData <- data.frame(Time = Time, Incidence = incidence)
  if (!is.null(dates)){
    epiData$Time <- dates
    xlabtype <- KerDateTicker(x = dates, spacing = datelab)
  }

    if(is.null(smooth)){
    mulayer <- ggplot2::geom_blank()
    CIlayer <- ggplot2::geom_blank()
  } else{
    thetahat <- smooth$thetahat
    K <- length(thetahat)
    B <- Rcpp_KercubicBspline(x = seq_len(numdays), lower = 1, upper = numdays,
                              K = K)
    musmooth <- as.numeric(exp(B %*% thetahat))
    Sighat <- smooth$Sighat
    CImu <- function(t, alpha = 0.05) {
      bt <- as.numeric(Rcpp_KercubicBspline(t, lower = 1, upper = numdays, K = K))
      zalpha <- stats::qnorm(1 - 0.5 * alpha, lower.tail = T)
      sdd <- sqrt(sum((bt * Sighat) %*% bt))
      logCIlow <- sum(thetahat * bt) - zalpha * sdd
      logCIup <-  sum(thetahat * bt) + zalpha * sdd
      output <- cbind(exp(logCIlow), exp(logCIup))
      return(output)
    }
    CImusmooth <- t(sapply(seq_len(numdays), CImu))
    smoothData <- data.frame(Time = Time, musmooth, CImusmooth)
    if (!is.null(dates)){
      smoothData$Time <- dates
    }
    colnames(smoothData) <- c("Time", "musmooth", "CIlow", "CIup")
    mulayer <- ggplot2::geom_line(data = smoothData, ggplot2::aes(y = musmooth),
                                  color = smoothcol, size = 0.8)
    CIlayer <-  ggplot2::geom_ribbon(data = smoothData,
                                     ggplot2::aes(ymin = CIlow, ymax = CIup),
                                     alpha = 0.3, fill = smoothcol)
  }


  plotout <- ggplot2::ggplot(data = epiData, ggplot2::aes(x = Time)) +
    ggplot2::geom_bar(stat = "identity", width = barwidth,
                      ggplot2::aes(y = Incidence), col = col, fill = col) +
    ggplot2::ggtitle(title) + xlabtype +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = xtickangle,
                                          vjust = 0.55)) + mulayer + CIlayer
  return(plotout)
}
