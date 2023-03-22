#' Summarize the estimated reproduction number
#'
#' @description This routine can be used to summarize estimation results for
#' related to the reproduction number.
#'
#' @param object An object of class \code{Rt}.
#' @param ... Further arguments to be passed.
#'
#' @return A summary output.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @export

summary.Rt <- function(object, ...){
  if(!inherits(object, "Rt"))
    stop("object must be an Rt object")
  totdays <- length(object$incidence)
  timing <- object$LPS_elapsed
  method <- object$method
  if(method == "LPSMAP"){
    method <- "Maximum a posteriori (MAP)"
  } else{
    method <- "MCMC (with Langevin diffusion)"
  }
  optim_method <- object$optim_method
  if(optim_method == "NelderMead"){
    optim <- "Nelder-Mead"
  } else if (optim_method == "HillClimb"){
    optim <- "Gradient ascent"
  }
  converged <- object$optimconverged
  meanR <- round(mean(object$RLPS$R[-(1:7)]),3)
  maxR <- round(max(object$RLPS$R[-(1:7)]),3)
  minR <- round(min(object$RLPS$R[-(1:7)]),3)
  #--- Print output
  cat("Estimation of the reproduction number with Laplacian-P-splines \n")
  cat("-------------------------------------------------------------- \n")
  cat("Total number of days:         ", totdays,"\n")
  cat("Routine time (seconds):       ", timing, "\n")
  cat("Method:                       ", method, "\n")
  cat("Hyperparam. optim method:     ", optim, "\n")
  cat("Hyperparam. optim convergence:", converged, "\n")
  cat("Mean reproduction number:     ", format(meanR, nsmall = 3), "\n")
  cat("Min  reproduction number:     ", format(minR, nsmall = 3), "\n")
  cat("Max  reproduction number:     ", format(maxR, nsmall = 3), "\n")
  cat("-------------------------------------------------------------- \n")
}
