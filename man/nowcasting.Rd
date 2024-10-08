% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nowcasting.R
\name{nowcasting}
\alias{nowcasting}
\title{Nowcasting and estimation of occurred but not yet reported events}
\usage{
nowcasting(data, day.effect = TRUE, ref.day = "Monday", verbose = TRUE)
}
\arguments{
\item{data}{A data frame containing the data for each time and delay combination with
the following 6 columns. The first column is a numeric variable associated to the calendar date.
The second column is a numeric variable indicating the delay of reporting. The third column
corresponds to the calendar date of the event (e.g. death) and the fourth column to the
calendar date at which the event of interest was reported. The fifth column indicates the
number of cases for each time and delay combination. Finally, the sixth column indicates whether
the cases are already reported or not yet reported. To see an example of such a data structure
type \code{data("cov19mort2021")} and then \code{head(cov19mort2021)}. This will illustrate the
data structure used for nowcasting the daily number of cases based on mortality data for Belgium
in 2021.}

\item{day.effect}{If TRUE (defaut), include the day of the week effect.}

\item{ref.day}{If \code{day.effect = TRUE}, then the reference category for
the day of the week must be specified. The default is "Monday".}

\item{verbose}{Show summary output of nowcast in console? Default is TRUE.}
}
\value{
A list with the following components:
\itemize{
 \item data: The data frame used as an input.
 \item cases.now: A data frame containing the nowcasting results with \eqn{95\%} prediction intervals.
 \item delay: A data frame containing the two-dimensional delay distribution.
 \item lambda_estim: Estimated penalty parameters of the P-splines model.
 \item phi_estim: Estimated overdispersion parameter from the negative binomial model.
 \item Kt: The number of B-spline basis functions in the time dimension.
 \item Kd: The number of B-spline basis functions in the delay dimension.
 \item xi_mode: The mode of the conditional posterior distribution of the latent vector.
 }
}
\description{
This routine can be used to estimate cases that have not yet been reported also known as
occurred-but-not-yet-reported-events. Daily cases are typically subject to reporting
delays, so that the reported number of infected individuals is not always reflecting the
true epidemic status. Nowcasting aims to correct this underreporting phenomenon by estimating
the number of infections that have occurred but that have not yet been reported. The
latter number is then combined with the already reported cases and interpreted as a nowcast or
prediction for the true epidemic status regarding the number of daily cases. The routine
is anchored around Laplacian-P-splines in an epidemic context (Gressani et al. 2022) and the
detailed methodology can be found in Sumalinab et al. (2023).
}
\examples{
# data("cov19mort2021")
# ncast <- nowcasting(data = cov19mort2021, day.effect = FALSE)
# plot(ncast) # Show nowcasted cases
# plot(ncast, type = "delay") # Show contour of delay distribution

}
\references{
Gressani, O., Wallinga, J., Althaus, C. L., Hens, N. and Faes, C.
 (2022). EpiLPS: A fast and flexible Bayesian tool for estimation of the
 time-varying reproduction number. \emph{Plos Computational Biology},
 \strong{18}(10): e1010618.

Sumalinab, B., Gressani, O., Hens, N. and Faes, C. (2023). Bayesian
 nowcasting with Laplacian-P-splines. MedRxiv preprint.
 \url{https://www.medrxiv.org/content/10.1101/2022.08.26.22279249v2}
}
\author{
Bryan Sumalinab (writing) and Oswaldo Gressani (editing).
}
