#' Adjusting a positive definite matrix
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd


KerAdjustPD <- function(x, eigentol = 1e-06, correct = TRUE){
  if (!is.matrix(x))
    stop("Input argument must be a matrix")
  if (!isSymmetric(x))
    stop("Input matrix must be symmetric")
  if (ncol(x) != nrow(x))
    stop("Input matrix must be square")
  if (eigentol <= 0)
    stop("eigentol must be a positive number")
  if (!is.logical(correct))
    stop("correct must be either T or F")
  eigvals <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  checkPD <- !any(eigvals < eigentol)
  if (correct == TRUE) {
    if (checkPD == FALSE) {
      xcorrect <- x + diag(abs(min(eigvals)) + 1e-04, ncol(x))
    }
    else {
      xcorrect <- x
    }
    list_out <- list(isPD = checkPD, PD = xcorrect)
  }
  else {
    list_out <- list(isPD = checkPD)
  }
  return(list_out)
}
