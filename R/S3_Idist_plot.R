#' Plot the interval distribution from an Idist object
#'
#' @description This routine plots the interval distribution based on an Idist
#' object.
#'
#' @usage
#' \method{plot}{Idist}(x, barcol = "firebrick", denscol = "pink", denstransparent = 0.5,
#'  barwidth = 0.30, title = NULL, themetype = c("gray","classic","light","dark"),
#'  titlesize = 15, xtitlesize = 13, ytitlesize = 13, ...)
#'
#' @param x An object of class \code{Idist}.
#' @param barcol Color of the discretized interval distribution.
#' @param denscol Color of the probability density function (pdf).
#' @param denstransparent The transparency of the pdf.
#' @param barwidth Width of the bars for the discrete distribution.
#' @param title Title of the plot.
#' @param themetype Theme of the plot.
#' @param titlesize Size of the plot title. Default is 15.
#' @param xtitlesize Size of title and text on x axis. Default is 13.
#' @param ytitlesize Size of title and text on y axis. Default is 13.
#' @param ... Further arguments to be passed to plot.
#'
#' @examples
#' x <- Idist(mean = 3, sd = 1.6, dist = "weibull")
#' plot(x)
#'
#' @return A plot of the interval distribution.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

plot.Idist <- function(x, barcol = "firebrick", denscol = "pink",
                 denstransparent = 0.5, barwidth = 0.30, title = NULL,
                 themetype = c("gray","classic","light","dark"),
                 titlesize = 15, xtitlesize = 13, ytitlesize = 13, ...){

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

  idx <- seq_len(length(x$pvec))
  y <- x$pvec
  data <- data.frame(idx = idx, y = y)

  if(length(x) == 4){ # Plot probability density function
  dom <- seq(0, length(y) + 1, length = 500)
  if(is.null(title)){ # Add default title
    if(x$dist == "gamma"){
      dens <- stats::dgamma(dom, shape = x$shape, rate = x$rate)
      data2 <- data.frame(dom = dom, dens = dens)
      densadd <- ggplot2::geom_area(data = data2, ggplot2::aes(x=dom, y=dens),
                                    fill=denscol, alpha = denstransparent)
      titleplot <- paste0("Interval distribution based on a Gamma with shape=",
                           round(x$shape,3)," and rate=", round(x$rate,3))
    } else if(x$dist == "weibull"){
      dens <- stats::dweibull(dom, shape = x$shape, scale = x$scale)
      data2 <- data.frame(dom = dom, dens = dens)
      densadd <- ggplot2::geom_area(data = data2, ggplot2::aes(x=dom, y=dens),
                                    fill=denscol, alpha = denstransparent)
      titleplot <- paste0("Interval distribution based on a Weibull with shape=",
                          round(x$shape,3)," and scale=", round(x$scale,3))
    } else if(x$dist == "lognorm"){
      dens <- stats::dlnorm(dom, meanlog=x$location, sdlog = x$scale)
      data2 <- data.frame(dom = dom, dens = dens)
      densadd <- ggplot2::geom_area(data = data2, ggplot2::aes(x=dom, y=dens),
                                    fill=denscol, alpha = denstransparent)
      titleplot <- paste0("Interval distribution based on a LogNormal with location=",
                          round(x$location,3),
                          " and scale=", round(x$scale,3))
    }
  } else{
    titleplot <- title
    }
  } else{
    densadd <- ggplot2::geom_blank()
    titleplot <- "Interval distribution"
  }

  Iplot <- ggplot2::ggplot(data = data, ggplot2::aes(x = idx, y = y)) +
    ggplot2::scale_x_discrete(name = "Time", limits = as.character(idx)) +
    ggplot2::geom_bar(stat = "identity", width = barwidth,
                      color = barcol, fill = barcol) +
    ggplot2::ylab("Probability") +
    ggplot2::ggtitle(titleplot) +
    themeval +
    densadd +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = titlesize),
      axis.title.x = ggplot2::element_text(size = xtitlesize),
      axis.title.y = ggplot2::element_text(size = ytitlesize),
      axis.text.x = ggplot2::element_text(size = xtitlesize),
      axis.text.y = ggplot2::element_text(size = ytitlesize)
    )

  return(Iplot)
  }



