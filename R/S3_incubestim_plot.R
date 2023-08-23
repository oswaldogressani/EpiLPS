#' Plot the estimated incubation distribution
#'
#' @description This routine can be used to plot the estimated incubation
#' distribution based on an object of class \code{incubestim}.
#'
#' @usage
#' \method{plot}{incubestim}(x, type = c("pdf", "cdf", "incubwin"), ...)
#'
#' @param x An object of class \code{incubestim}.
#' @param type The type of plot. Default is "pdf" for the plot of the density
#' function. Setting it to "cdf" returns the cumulative distribution function
#' and "incubwin" gives a bar plot showing the width of the incubation window.
#' @param ... Further arguments to be passed to plot.
#'
#' @examples
#' set.seed(123)
#' simdat <- incubsim(n = 30, tmax = 20) # Simulate incubation data
#' data <- simdat$Dobsincub              # Incubation bounds
#' incubfit <- estimIncub(x = data, niter = 500, tmax = 20, verbose = TRUE)
#' plot(incubfit)
#'
#' @return
#'\itemize{
#'  \item{The probability density function of the estimated incubation period with
#'   95\% credible envelope.}
#'  \item{The cumulative distribution function of the estimated incubation period with
#'   95\% credible envelope.}
#' \item{A bar plot showing the width of the incubation windows.}
#' }
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

plot.incubestim <- function(x, type = c("pdf", "cdf", "incubwin"), ...) {
  n <- nrow(x$incubdat)

  # Plot incubation windows
  index <- seq_len(n)
  lower <- x$incubdat$tL
  upper <- x$incubdat$tR
  df1 <- data.frame(index = index, lower = lower, upper = upper)

  incub <- ggplot2::ggplot(data = df1, ggplot2::aes(x = index)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = lower, ymax = upper),
      color = "darkslateblue",
      linewidth = 1
    ) +
    ggplot2::xlab("Individual index number") +
    ggplot2::ylab("Incubation bound") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 17),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 9)
    )

  # Plot probability density function with 95%CI
  tdom <- x$tg
  fhat <- x$ftg
  flow <- apply(x$ftgMCMC, 2, stats::quantile, probs = 0.025)
  fup <- apply(x$ftgMCMC, 2, stats::quantile, probs = 0.975)
  df2 <- data.frame(tdom = tdom, fhat = fhat, flow = flow, fup = fup)

  densplot <-  ggplot2::ggplot(data = df2, ggplot2::aes(x = tdom)) +
    ggplot2::xlab("Incubation period (days)") +
    ggplot2::ylab("Probability density function") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = flow, ymax = fup), alpha = 0.5,
                         fill = "gray60") +
    ggplot2::geom_line(ggplot2::aes(y = fhat), color = "dodgerblue2",
                       linewidth = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 17),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 9)
    )

  # Plot cumulative distribution function with 95% CI
  dt <- tdom[2] - tdom[1]
  Fhat <- cumsum(fhat * dt)
  nMCMC <- nrow(x$ftgMCMC)
  FMCMC <- matrix(0, nrow = nMCMC, ncol = ncol(x$ftgMCMC))
  for (j in 1:nMCMC) {
    FMCMC[j, ] <- cumsum(x$ftgMCMC[j, ] * dt)
  }
  Flow <- apply(FMCMC, 2, stats::quantile, probs = 0.025)
  Fup <- apply(FMCMC, 2, stats::quantile, probs = 0.975)
  df3 <- data.frame(tdom = tdom, Fhat = Fhat, Flow = Flow, Fup = Fup)

  cdfplot <- ggplot2::ggplot(data = df3, ggplot2::aes(x = tdom)) +
    ggplot2::xlab("Incubation period (days)") +
    ggplot2::ylab("Cumulative distribution function") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Flow, ymax = Fup), alpha = 0.5,
                         fill = "gray60") +
    ggplot2::geom_line(ggplot2::aes(y = Fhat), color = "brown3",
                       linewidth = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 17),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 9)
    )

  if(match.arg(type) == "pdf"){
    densplot
  } else if(match.arg(type) == "cdf"){
    cdfplot
  } else{
    incub
  }
  # gridExtra::grid.arrange(incub, densplot, cdfplot, nrow = 1)
}




















