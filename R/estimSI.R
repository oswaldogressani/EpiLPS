#' Nonparametric serial interval estimation
#'
#' @description
#' The serial interval (SI) is defined as the time between illness onset in
#' a primary case (infector) and illness onset in a secondary case (infectee).
#' This routine estimates different features of the SI distribution based on
#' transmission pair data. The input \code{x} should be a data frame with two
#' columns containing the left and right boundary of the serial interval window
#' for each transmission pair. A resampling technique (bootstrap) is used to
#' draw samples from a nonparametric estimate of the cumulative distribution
#' function of the serial interval.
#'
#' @usage estimSI(x, B = 5000, verbose = TRUE)
#'
#' @param x A data frame with two columns based on transmission pair data. The
#'  number of rows is the number of transmission pairs (sample size). The first
#'  column should be the left boundary of the SI window and the second column
#'  should be the right boundary of the SI window. NA values are not permitted.
#' @param B The number of bootstrap samples.
#' @param verbose Print basic summary features of the SI distribution?
#'
#' @return A list containing basic summary features with point estimates and
#' interval estimates.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Gressani, O. and Hens, N. (2024). Nonparametric serial interval
#' estimation. \emph{MedRxiv preprint}.
#'
#' @examples
#' data(influenza2009serial)
#' SIfit <- estimSI(x = influenza2009serial)
#'
#' @export

estimSI <- function(x, B = 5000, verbose = TRUE){
  tic <- proc.time()
  sidat <- as.matrix(x)
  sifit <- Rcpp_Kerserialint(sidat, B = B)
  out <- as.data.frame(sifit$estimres)
  toc <- proc.time() - tic

  if(isTRUE(verbose)){
    cat("---------------------------------------------------------------------\n")
    cat("EpiLPS nonparametric serial interval (SI) estimation \n")
    cat("---------------------------------------------------------------------\n")
    cat("Time elapsed (in seconds):", round(toc[3],2), "\n")
    cat("Sample size:", nrow(x), "(Bootstrap:", paste0(B,")"), "\n")
    cat("Mean   SI  :", format(round(out$Estim[1], 2),nsmall = 2), "days (95% CI:",
        paste0(format(round(out$CI95p_l[1],2), nsmall = 2),"-",
               format(round(out$CI95p_r[1],2), nsmall = 2), ")"), "\n")
    cat("Median SI  :", format(round(out$Estim[5], 2), nsmall = 2), "days (95% CI:",
        paste0(format(round(out$CI95p_l[5],2), nsmall = 2),"-",
               format(round(out$CI95p_r[5], 2), nsmall = 2),")"), "\n")
    cat("St.dev SI  :", format(round(out$Estim[2], 2), nsmall = 2), "days (95% CI:",
        paste0(format(round(out$CI95p_l[2],2), nsmall = 2),"-",
               format(round(out$CI95p_r[2], 2), nsmall = 2),")"), "\n")
    cat("---------------------------------------------------------------------\n")

  }

  outlist <- list(estim = out, data = x, B = B)
  return(outlist)
}
