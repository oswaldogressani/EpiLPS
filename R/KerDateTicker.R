#' Construction of date labels for x-axis
#' @param x A vector of class Date.
#' @param spacing The spacing required for the date labels.
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerDateTicker <- function(x, spacing = "7d") {
  # Kernel routine
  # Author: Oswaldo Gressani (oswaldo_gressani@hotmail.fr)
  if (!inherits(x, "Date"))
    stop("x is not of class 'Date'.")
  if (!is.character(spacing))
    stop("spacing must be a character, e.g. 7d for 7 days.")
  spacing_len <- nchar(spacing)
  if (spacing_len > 3 || spacing_len < 2)
    stop("x should be of format #d with # in 1-31 (days) or *m with
       * in 1-12 (months).")
  timescale <- substr(spacing, spacing_len, spacing_len)
  if (timescale == "d") {
    timelab <- "days"
    daynums <- as.numeric(substr(spacing, 1, spacing_len - 1))
    if (daynums < 1 || daynums > 31)
      stop("Number of days should be between 1 and 31.")
    ggtext <- eval(parse(
      text = paste0("ggplot2::scale_x_date(date_breaks = '", daynums, " days')")
    ))
  } else if (timescale == "m") {
    timelab <- "months"
    monthnums <- as.numeric(substr(spacing, 1, spacing_len - 1))
    if (monthnums < 1 || monthnums > 12)
      stop("Number of months should be betwwen 1 and 12.")
    ggtext <- eval(parse(
      text = paste0("ggplot2::scale_x_date(date_breaks = '",
                    monthnums," months')")
    ))
  } else{
    stop("spacing should be either day(s) 'd' or month(s) 'm'.")
  }
  return(ggtext)
}
















