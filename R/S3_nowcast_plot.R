#' Plot of nowcasted cases and delay distribution
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
  data1 <- data[data$Reported != "Not yet reported", ]
  data1 <- data1[data1$Date > max(data1$Date) - (2 * D), ]
  data1 <- stats::aggregate(Cases ~ Date + Reported, data = data1, FUN = sum)
  data2 <- data.frame("Date" = data_CI$Date, "Reported" = "Nowcast",
                      "Cases" = utils::tail(x$cases.now$y, n = D) -
                        utils::tail(data1$Cases, n = D))
  data3 <- rbind(data1,data2)
  CI95_low <- c(rep(NA,nrow(data3)-D), data_CI$CI95L)
  CI95_up <- c(rep(NA,nrow(data3)-D), data_CI$CI95R)
  data3 <- cbind(data3, CI95_low, CI95_up)
  datarep <- data3[data3$Reported == "Reported",]
  datanow <- data3[data3$Reported == "Nowcast",]
  datanow$Cases <- datanow$Cases + utils::tail(datarep$Cases,D)

  layer_rep <- ggplot2::ggplot(data = datarep, ggplot2::aes(x = Date, y = Cases)) +
    ggplot2::geom_col(ggplot2::aes(fill = "Reported"), width = 1,
                      position = ggplot2::position_stack(reverse = TRUE))
  layer_now <- ggplot2::geom_col(data = datanow, ggplot2::aes(fill = "Nowcast"),
                                 width = 1, position = ggplot2::position_stack(reverse = TRUE))
  layer_rep2 <- ggplot2::geom_col(ggplot2::aes(fill = "Reported"), width = 1,
                                  position = ggplot2::position_stack(reverse = TRUE))
  layer_nowCI <- ggplot2::geom_crossbar(data = datanow, mapping = ggplot2::aes(x = Date,
                        y = Cases, ymin = `CI95_low`, ymax = `CI95_up`,
                        fill = "95% CI"), col = NA, width = 1)

  caseplot <- layer_rep + layer_now + layer_rep2 + layer_nowCI +
    ggplot2::scale_fill_manual(values = c(grDevices::adjustcolor("orange", alpha = 0.5),
                                          "gray50", "blue"), name = "")

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
