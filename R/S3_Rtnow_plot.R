#' Plot the nowcasted reproduction number
#'
#' @description This routine can be used to plot the nowcasted reproduction
#' number based on an object of class \code{Rtnow}.
#'
#' @usage
#' \method{plot}{Rtnow}(x, datelab = "7d", cilevel = 0.95, xtickangle = 0, legendpos = "right",
#'  nowcastcol = "cornflowerblue", ...)
#'
#' @param x An object of class \code{Rtnow}.
#' @param datelab Spacing for the ticks on the x-axis.
#' @param cilevel Level of the credible interval.
#' @param xtickangle Angle of the x-ticks. Default is 0 (horizontal).
#' @param legendpos Position of the legend.
#' @param nowcastcol Color of the nowcasted \eqn{R_t} curve.
#' @param ... Further arguments to be passed to plot.
#'
#' @return A plot of the nowcasted time-varying reproduction number.
#'
#' @author Bryan Sumalinab (writing) and Oswaldo Gressani (editing).
#'
#' @export


plot.Rtnow <- function(x, datelab = "7d", cilevel = 0.95, xtickangle = 0,
                       legendpos = "right", nowcastcol = "cornflowerblue", ...){
    method <- x$method
    Rnow <- x$Rnow
    tstart <- 7
    Rnow.tstart <- Rnow[-(1:tstart), ]
    Rnow.tstart$status <- factor(Rnow.tstart$status,
                                 levels = c("observed", "nowcasted"),
                                 labels = c("observed", "nowcasted"))
    xlabtype <- KerDateTicker(x = Rnow.tstart$Time, spacing = datelab)
    Rnowcasted <- Rnow.tstart[Rnow.tstart$status == "nowcasted",]
    Time <- Rnow.tstart$Time
    R <- Rnow.tstart$R

    if(cilevel == 0.95){
    CIg1 <-  eval(parse(text = 'ggplot2::geom_ribbon(ggplot2::aes(ymin = Rq0.025, ymax = Rq0.975), alpha = 0.5, fill = "gray")'))
    CIg2 <-  eval(parse(text = 'ggplot2::geom_ribbon(data = Rnowcasted, ggplot2::aes(x = Time, ymin = Rq0.025, ymax = Rq0.975, fill = "nowcasted"), fill = nowcastcol, alpha = 0.3)'))
    } else if (cilevel == 0.90){
    CIg1 <-  eval(parse(text = 'ggplot2::geom_ribbon(ggplot2::aes(ymin = Rq0.05, ymax = Rq0.95), alpha = 0.5, fill = "gray")'))
    CIg2 <-  eval(parse(text = 'ggplot2::geom_ribbon(data = Rnowcasted, ggplot2::aes(x = Time, ymin = Rq0.05, ymax = Rq0.95, fill = "nowcasted"), fill = nowcastcol, alpha = 0.3)'))
    } else{
      stop("cilevel must be either 0.95 or 0.90")
    }

    g1 <- ggplot2::ggplot(data = Rnow.tstart, ggplot2::aes(x = Time, y = R)) +
      CIg1 +
      ggplot2::geom_line(ggplot2::aes(y = R), linewidth = 1, color = "black") +
      xlabtype +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = xtickangle, vjust = 0.55)) +
      ggplot2::ylim(0, NA)

    g2 <- CIg2

    g3 <-  ggplot2::geom_line(data = Rnowcasted, ggplot2::aes(x = Time, y = R,
                              color = "nowcasted"), linewidth = 1)

    plotout <- g1 + g3 + g2 +
      ggplot2::labs(x = "Time", y = "R") +
      ggplot2::geom_hline(yintercept = 1, linetype = "dotted", linewidth = 0.8,
                          color = "gray46") +
      ggplot2::ylab("R") +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = legendpos) +
      ggplot2::ggtitle('Estimated R') +
      ggplot2::scale_color_manual(values = c("nowcasted" = nowcastcol)) +
      ggplot2::scale_fill_manual(values = c("nowcasted" = nowcastcol))

    return(plotout)
  }
