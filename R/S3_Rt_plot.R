#' Plot the estimated reproduction number
#'
#' @description This routine can be used to plot the estimated reproduction
#' number based on an object of class \code{Rt}.
#'
#' @usage
#' \method{plot}{Rt}(x, datelab = "7d", cilevel = 0.95, col = "black", cicol = "gray",
#'  xtickangle = 0, legendpos = "right", title = "Estimated R",
#'  addfit = c("none","Cori","WT"), theme = "gray", timecut = 0,...)
#'
#' @param x An object of class \code{Rt}.
#' @param datelab Spacing for the ticks on the x-axis.
#' @param cilevel Level of the credible interval.
#' @param col Color of the fitted \eqn{R_t} curve for LPS.
#' @param cicol Color for shading the credible envelope.
#' @param xtickangle Angle of the x-ticks. Default is 0 (horizontal).
#' @param legendpos Position of the legend.
#' @param title Title of the plot.
#' @param addfit Should an additional \eqn{R_t} fit be added?
#' @param theme Theme, either "gray", "classic", "light", "dark"
#' @param timecut Cut time points on plot.
#' @param ... Further arguments to be passed to plot.
#'
#' @examples
#' si <- c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
#' epidemic <- episim(si = si, Rpattern = 2, endepi = 30)
#' epifit <- estimR(incidence = epidemic$y, K = 30, si = si)
#' plot(epifit)
#'
#' @return A plot of the fitted time-varying reproduction number.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export


plot.Rt <- function(x, datelab = "7d", cilevel = 0.95, col = "black",
                        cicol = "gray", xtickangle = 0, legendpos = "right",
                        title = "Estimated R", addfit = c("none","Cori","WT"),
                        theme = "gray", timecut = 0,...) {
    if(!inherits(x, "Rt"))
      stop("x must be an Rt object")
    if(cilevel != 0.90 & cilevel != 0.95)
      stop("cilevel must be either 0.90 or 0.95")

    # Extract R data frame for LPS method
    LPSRout <- x$RLPS[-(1:(7 + timecut)),]
    Corilayer <- ggplot2::geom_blank()
    Tlen <- nrow(LPSRout)
    xlabtype <- ggplot2::geom_blank()
    if (inherits(LPSRout$Time, "Date"))
      xlabtype <- KerDateTicker(x = LPSRout$Time, spacing = datelab)

    LPSlayer <-
      ggplot2::ggplot(data = LPSRout, ggplot2::aes(x = Time)) +
      ggplot2::ggtitle(title) + ggplot2::xlab("Time") +
      ggplot2::ylab("R") +
      ggplot2::geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.8,
                          color = "gray46") +
      ggplot2::ylim(0, NA) +
      ggplot2::geom_ribbon(ggplot2::aes(
        ymin = eval(parse(text = paste0("`Rq", (1 - cilevel) * 0.5, "`"))),
        ymax = eval(parse(text = paste0("`Rq", 1 - (1 - cilevel) * 0.5, "`")))),
        alpha = 0.5, fill = cicol) + xlabtype +
      ggplot2::geom_line(ggplot2::aes(y = `R`, linetype = "EpiLPS"),
                         color = col, linewidth = 0.8) +
      eval(parse(text = paste0("ggplot2::theme_",theme,"()"))) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = xtickangle,
                                                         vjust = 0.55)) +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = legendpos)

    if(is.data.frame(x$RCori)) {# Compute Cori plot layer
      'Time' <- 'R' <- 'Mean(R)' <- NULL
      'Quantile.0.025(R)' <- 'Quantile.0.975(R)' <- NULL
      'Quantile.0.05(R)' <- 'Quantile.0.95(R)' <- NULL
      if(timecut > 0){
        CoriRout <- x$RCori[-(1:timecut),]
      }else{
      CoriRout <- x$RCori
      }
      CoriRout$Time <- LPSRout$Time
      Corilayer <-  ggplot2::geom_line(data = CoriRout, ggplot2::aes(
        x = Time, y = `Mean(R)`, color = "EpiEstim", linetype = "EpiEstim"),
        linewidth = 0.8, linetype = "twodash")
      CoriCI <- ggplot2::geom_ribbon(data = CoriRout,ggplot2::aes(
        ymin = eval(parse(
          text = paste0("`Quantile.", (1 - cilevel) * 0.5, "(R)`"))),
        ymax = eval(parse(
          text = paste0("`Quantile.", 1 - (1 - cilevel) * 0.5, "(R)`")))),
        alpha = 0.5, fill = "cornflowerblue")
      Coricalled <- TRUE
    } else{
      Coricalled <- FALSE
    }

    if(is.data.frame(x$RWT)) {# Compute Wallinga-Teunis plot layer
      'Time' <- 'R' <- 'Mean(R)' <- NULL
      'Quantile.0.025(R)' <- 'Quantile.0.975(R)' <- NULL
      WTRout <- x$RWT
      WT_tstart <- WTRout$t_end[1]
      WTRout$Time <- x$RLPS$Time[WT_tstart:WTRout$t_end[nrow(WTRout)]]
      if ((timecut + 7) > WT_tstart) {
        WTRout <- WTRout[-(1:(timecut + 8 - WT_tstart)), ]
      }
      WTlayer <-  ggplot2::geom_line(data = WTRout, ggplot2::aes(
        x = Time, y = `Mean(R)`, color = "EpiEstimWT", linetype = "EpiEstimWT"),
        linewidth = 0.8, linetype = "twodash")
      WTCI <- ggplot2::geom_ribbon(data = WTRout,ggplot2::aes(
        ymin = `Quantile.0.025(R)`, ymax = `Quantile.0.975(R)`),
        alpha = 0.5, fill = "coral2")
      WTcalled <- TRUE
    } else{
      WTcalled <- FALSE
    }

    if(match.arg(addfit) == "none"){# Compute plot to output
      plotout <- LPSlayer
    } else if(match.arg(addfit) == "Cori"){
      if(isFALSE(Coricalled))
        stop("Cori method has not been called. Use CoriR = TRUE in the
             epilps() routine to call it.")
      plotout <- LPSlayer + CoriCI + Corilayer +
        ggplot2::scale_color_manual(values = c("EpiEstim" = "cornflowerblue"))
    } else if (match.arg(addfit) == "WT"){
      if(isFALSE(WTcalled))
        stop("Wallinga-Teunis method has not been called. Use WTR = TRUE in the
             epilps() routine to call it.")
      if(cilevel == 0.90)
        warning("The CI level for the Wallinga-Teunis method is at 95%.")
      plotout <- LPSlayer + WTCI + WTlayer +
        ggplot2::scale_color_manual(values = c("EpiEstimWT" = "coral2"),
                                    labels = c("EpiEstimWT"="WT"))
    }

    return(plotout)

  }

