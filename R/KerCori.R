#' Calling Cori method to estimate reproduction number
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerCori <- function(Dobs, sinter){
  val <- suppressMessages(suppressWarnings(
    EpiEstim::estimate_R(
      incid  = Dobs,
      method = "non_parametric_si",
      config = EpiEstim::make_config(list(si_distr = c(0, sinter))))$R))
  return(val)
}
