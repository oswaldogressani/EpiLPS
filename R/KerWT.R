#' Calling Wallinga-Teunis method to estimate reproduction number
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerWT <- function(Dobs, sinter){
  n <- length(Dobs)
  val <- suppressMessages(suppressWarnings(
    EpiEstim::wallinga_teunis(
      incid  = Dobs,
      method = "non_parametric_si",
      config = EpiEstim::make_config(list(si_distr = c(0, sinter),
                                          t_start = seq(2, n - 6),
                                          t_end = seq(8, n))))$R))
  return(val)
}
