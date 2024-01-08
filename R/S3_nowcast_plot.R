#' Plot of nowcasted cases
#'
#' @description This routine plots the reported cases, the nowcasted cases and
#' associated \eqn{95\%} credible interval for the nowcast. In addition, it
#' can be used to show the contour plot of the estimated delay distribution.
#'
#' @usage
#' \method{plot}{nowcasted}(x, type = c("cases", "delay"), ...)
#'
#' @param x An object of class \code{nowcasted}.
#' @param type Either "cases" (default) to show the nowcasted cases or
#'  "delay" to show the contour plot of the estimated delay distribution.
#' @param ... Further arguments to be passed to plot.
#'
#' @return A plot of the reported and nowcasted cases.
#'
#' @author Bryan Sumalinab (writing) and Oswaldo Gressani (editing).
#'
#' @export

plot.nowcasted <- function(x, type = c("cases", "delay"), ...){
  data <- x$data
  data_CI <- stats::na.omit(x$cases.now)
  D <- max(data$d)
  Date <- data$Date
  Cases <- data$Cases
  Reported <- data$Reported
  y <- x$cases.now$y
  CI95_low <- data_CI$CI95L
  CI95_up <- data_CI$CI95R
  data1 <- data[data$Reported != "Not yet reported", ]
  data1$Reported <- factor(data1$Reported,
    levels = c("Reported", "Nowcast", "95% CI"),
    labels = c("Reported", "Nowcast", "95% CI")
  )
  data1 <- data1[data1$Date > max(data1$Date) - (2 * D), ]
  data1 <- stats::aggregate(Cases ~ Date + Reported, data = data1, FUN = sum)
  data2 <- data.frame("Date" = data_CI$Date, "Reported" = "Nowcast",
                      "Cases" = utils::tail(x$cases.now$y, n = D) -
                        utils::tail(data1$Cases, n = D))
  data3 <- rbind(data1,data2)

  # Plot of nowcasted cases
  caseplot <- ggplot2::ggplot(data = data3,
                              mapping = ggplot2::aes(x = Date, y = Cases,
                                                       fill = Reported)) +
    ggplot2::geom_col(width = 1,
                      position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_fill_manual(values = c("blue",
                                          grDevices::adjustcolor("black",
                               alpha = 0.5), "orange"), name = "",
                               drop = FALSE) +
    ggplot2::geom_crossbar(data = data_CI, mapping = ggplot2::aes(x = Date,
                           y = y, ymin = `CI95_low`, ymax = `CI95_up`),
                  fill = grDevices::adjustcolor("orange", alpha = 0.5),
                  colour = NA, width = 1, inherit.aes = FALSE)

  # Contour plot of delay distribution
  Delay <- x$delay$Delay
  density <- x$delay$density
  delayplot <- ggplot2::ggplot(data = x$delay,
                            mapping = ggplot2::aes(Date, Delay, z = density)) +
    ggplot2::geom_contour() +
    ggplot2::geom_contour_filled()


  if(match.arg(type) == "cases"){# Compute plot to output
    plotout <- caseplot
  } else if(match.arg(type) == "delay"){
    plotout <- delayplot
  }
  return(plotout)
}
