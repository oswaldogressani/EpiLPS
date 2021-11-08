#' Plot the EpiLPS fitted epidemic curve and reproduction number
#'
#' @description This routine can be used to plot the estimated epidemic curve
#'  and reproduction number with EpiLPS.
#'
#' @usage
#' \method{plot}{epilps}(x, plotout = c("rt", "epicurve"), dates = NULL,
#'      datelab = c("7d", "1m", "3m", "6m"),
#'      overlayEpiestim = FALSE, Rtitle = "", epititle = "", rtcol = "red", cicol = "gray",
#'      transparency = 0.5, epicol = "red",
#'      incibars = FALSE, barwidth = 0.35,
#'      themetype = c("gray", "classic", "light", "dark"), ...)
#'
#' @param x An object of class \code{epilps}.
#' @param plotout The type of plot, either "rt" for showing the reproduction
#'  number or "epicurve" for showing the epidemic curve.
#' @param dates A vector of dates in format "YY-MM-DD".
#' @param datelab The spacing for ticks on the x-axis. Either 7 days, 1 month,
#'  3 months or 6 months.
#' @param overlayEpiestim Should the EpiEstim fit be overlayed?
#' @param Rtitle The title for the plot of R.
#' @param epititle The title for the plot of the epidemic curve.
#' @param rtcol Color for the reproduction number curve fit.
#' @param cicol Color for shading the credible envelope.
#' @param transparency Controls the transparency of the credible envelope.
#' @param epicol The color for the epidemic curve.
#' @param incibars Should the bars of the incidence time series be shown?
#' @param barwidth The barwidth associated to the incidence time series.
#' @param themetype Type of theme for the plot.
#' @param ... Further arguments to be passed to plot.
#'
#' @examples
#' si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
#' epidemic <- episim(serial_interval = si, Rpattern = 2)
#' epifit <- epilps(incidence = epidemic$y, K = 30, serial_interval = si)
#' gridExtra::grid.arrange(plot(epifit, title = TRUE),
#'                         plot(epifit, plotout = "epicurve", epicol = "blue",
#'                         title = TRUE), nrow = 2)
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

plot.epilps <- function(x, plotout = c("rt", "epicurve"), dates = NULL,
                        datelab = c("7d", "1m", "3m", "6m"),
                        overlayEpiestim = FALSE, Rtitle = "", epititle = "",
                        rtcol = "red", cicol = "gray", transparency = 0.5,
                        epicol = "red", incibars = FALSE, barwidth = 0.35,
                        themetype = c("gray","classic","light","dark"), ...) {

  n <- nrow(x$epifit)
  incidence <- x$incidence
  t_start <- seq(2, n - 6)
  t_end <- t_start + 6
  epiestim_fit <- suppressMessages(suppressWarnings(
                   EpiEstim::estimate_R(incidence, method = "non_parametric_si",
                   config = EpiEstim::make_config(list(si_distr =
                            c(0, x$serial_interval), t_start = t_start,
                            t_end = t_end)))))
  Repiestim <- epiestim_fit$R$`Mean(R)`

  if (x$ci_level == 0.95) {
    Repiestim_CIlow <- epiestim_fit$R$`Quantile.0.025(R)`
    Repiestim_CIup <- epiestim_fit$R$`Quantile.0.975(R)`
  } else if (x$ci_level == 0.90){
    Repiestim_CIlow <- epiestim_fit$R$`Quantile.0.05(R)`
    Repiestim_CIup <- epiestim_fit$R$`Quantile.0.95(R)`
  }

  tdom <- seq(8, n)
  R_estim <-   x$epifit[8:n, 1]
  RCI_low <- x$epifit[8:n, 2]
  RCI_up <-  x$epifit[8:n, 3]
  Rlps <- data.frame(tdom = tdom, R_estim = R_estim, RCI_low = RCI_low,
                     RCI_up = RCI_up, Repiestim = Repiestim,
                     Repiestim_CIlow = Repiestim_CIlow,
                     Repiestim_CIup = Repiestim_CIup)

  # Adapt x-axis labels if dates are provided
  if(!is.null(dates)){
    Rlps[, 1] <- dates[8:n]
    datelab <- match.arg(datelab)
    if(datelab == "7d"){
    xlabtype <- eval(parse(
      text = "ggplot2::scale_x_date(date_breaks = '7 days')"))
    } else if(datelab == "1m"){
      xlabtype <- eval(parse(
        text = "ggplot2::scale_x_date(date_breaks = '1 month')"))
    } else if(datelab == "3m"){
      xlabtype <- eval(parse(
        text = "ggplot2::scale_x_date(date_breaks = '3 months')"))
    } else if(datelab == "6m"){
      xlabtype <- eval(parse(
        text = "ggplot2::scale_x_date(date_breaks = '6 months')"))
    }
  } else{
    xlabtype <- eval(parse(text = "ggplot2::xlim(0,n)"))
  }

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


  out_type <- match.arg(plotout)
  if (out_type == "rt") {
    if (overlayEpiestim == FALSE) {
      # -- Plot of Rt with EpiLPS
      plotR_EpiLPS <- ggplot2::ggplot(data = Rlps, ggplot2::aes(x = tdom)) +
                      ggplot2::ggtitle(Rtitle) +
                      ggplot2::geom_ribbon(ggplot2::aes(ymin = RCI_low,
                                                        ymax = RCI_up),
                                           alpha = transparency, fill = cicol) +
                      ggplot2::geom_line(ggplot2::aes(y = R_estim),
                                         color = rtcol, size = 1.1) +
                      ggplot2::geom_hline(yintercept = 1, linetype = "dotted",
                                          size = 1.1) +
                      xlabtype +
                      ggplot2::xlab("Time") + ggplot2::ylab("R") + themeval +
                      ggplot2::theme(
                        plot.title = ggplot2::element_text(size = 15),
                        axis.title.x = ggplot2::element_text(size = 13),
                        axis.title.y = ggplot2::element_text(size = 13),
                        axis.text.x = ggplot2::element_text(size = 13),
                        axis.text.y = ggplot2::element_text(size = 13))
      return(plotR_EpiLPS)
    } else{
      # -- Plot of Rt with EpiLPS and EpiEstim
      colors <- c("EpiLPS" = rtcol, "EpiEstim" = "black")
      linetypes <- c("EpiLPS" = 1, "EpiEstim" = 2)

      plotR_EpiLPS <- ggplot2::ggplot(data = Rlps, ggplot2::aes(x = tdom)) +
                      ggplot2::ggtitle(Rtitle) +
                      ggplot2::geom_ribbon(ggplot2::aes(ymin = RCI_low,
                                                        ymax = RCI_up),
                                          alpha = transparency, fill = cicol) +
                      ggplot2::geom_line(ggplot2::aes(y = R_estim,
                                                      color = "EpiLPS",
                                                      linetype = "EpiLPS"),
                                         size = 1.1) +
                      ggplot2::geom_ribbon(ggplot2::aes(ymin = Repiestim_CIlow,
                                                        ymax = Repiestim_CIup),
                                           alpha = transparency,
                                           color = "black", linetype = "dashed",
                                           fill = NA) +
                      ggplot2::geom_line(ggplot2::aes(y = Repiestim,
                                                      color = "EpiEstim",
                                                      linetype = "EpiEstim"),
                                         size = 1.1) +
                      ggplot2::labs(x = "Time", y = "R", color = "Legend",
                                    linetype = "Legend") +
                      ggplot2::scale_color_manual(values = colors) +
                      ggplot2::scale_linetype_manual(values = linetypes) +
                      ggplot2::geom_hline(yintercept = 1, linetype = "dotted",
                                          size = 1.1) +
                      xlabtype + themeval +
                      ggplot2::theme(
                        plot.title = ggplot2::element_text(size = 15),
                        axis.title.x = ggplot2::element_text(size = 13),
                        axis.title.y = ggplot2::element_text(size = 13),
                        axis.text.x = ggplot2::element_text(size = 13),
                        axis.text.y = ggplot2::element_text(size = 13))
      return(plotR_EpiLPS)
    }

  }else if (out_type == "epicurve"){
    tdom <- seq_len(n)
    mu_estim <-   x$epifit[, 4]
    muCI_low <- x$epifit[, 5]
    muCI_up <-  x$epifit[, 6]
    mulps <- data.frame(tdom = tdom, mu_estim = mu_estim,
                        muCI_low = muCI_low, muCI_up = muCI_up)

    if(incibars == TRUE){
    incidence_bars <- eval(parse(
      text = "ggplot2::geom_bar(stat = 'identity', width = barwidth,
              color = 'black', fill = 'black')"))
    } else{
      incidence_bars <- eval(parse(
      text = "ggplot2::geom_bar(stat = 'identity', width = barwidth,
            color = NA, fill = NA)"))
    }
    if(!is.null(dates)){
      mulps[,1] <- dates
    }else{
      xlabtype <- eval(parse(text = "ggplot2::xlim(0,n+1)"))
    }

    #-- Plot incidence data and smoothed mean
    plot_incidence <- ggplot2::ggplot(data = mulps,
                                      ggplot2::aes(x = tdom, y = mu_estim)) +
                      incidence_bars +
                      ggplot2::ggtitle(epititle) +
                      ggplot2::geom_ribbon(ggplot2::aes(ymin = muCI_low,
                                                        ymax = muCI_up),
                                           alpha = transparency, fill = cicol) +
                      ggplot2::geom_line(ggplot2::aes(y = mu_estim),
                                         color = epicol, size = 1.1) +
                      xlabtype +
                      ggplot2::xlab("Time") +
                      ggplot2::ylab("Incidence") +
                      themeval +
                      ggplot2::theme(
                        plot.title = ggplot2::element_text(size = 15),
                        axis.title.x = ggplot2::element_text(size = 13),
                        axis.title.y = ggplot2::element_text(size = 13),
                        axis.text.x = ggplot2::element_text(size = 13),
                        axis.text.y = ggplot2::element_text(size = 13))
    return(plot_incidence)
  }

}

