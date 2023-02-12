#' Checking the time series of incidence data
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#' @noRd

KerIncidCheck <- function(Dobs) {
  # Kernel routine
  # Author: Oswaldo Gressani (oswaldo_gressani@hotmail.fr)
  if (!is.numeric(Dobs))
    stop("incidence must be numeric")
  n <- length(Dobs)
  if (n <= 7)
    stop("Incidence data is required for at least 8 days.")
  if (any(is.infinite(Dobs)))
    stop("Incidence data must contain finite values.")
  if (anyNA(Dobs)) {# Rule to deal with NAs
    nreplace <- sum(is.na(Dobs))
    NAloc <- which(is.na(Dobs))
    for (j in 1:nreplace) {
      if (1 < NAloc[j] && NAloc[j] < n) {# NA is interior
        if (!is.na(Dobs[NAloc[j] + 1])) {
          Dobs[NAloc[j]] <-
            round((Dobs[NAloc[j] - 1] + Dobs[NAloc[j] + 1]) * 0.5)
        } else{
          Dobs[NAloc[j]] <- Dobs[NAloc[j] - 1]
        }
      } else if (NAloc[j] == 1) {# NA in first position
        Dobs[1] <- round(mean(Dobs, na.rm = TRUE))
      } else if (NAloc[j] == n) {# NA in last position
        Dobs[n] <- Dobs[NAloc[j] - 1]
      }
    }
    warning("Incidence data contains NA values.")
  }
  if (any(Dobs < 0))
    stop("Incidence data contains negative values.")
  return(Dobs)
}
