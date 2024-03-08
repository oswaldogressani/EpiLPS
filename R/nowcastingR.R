#' Nowcasting the reproduction number
#'
#' @description
#' This routine can be used to nowcast the time-varying reproduction number. Daily cases are typically
#' subject to reporting delays, so that the reported number of infected individuals is not always reflecting the
#' true epidemic status. Nowcasting aims to correct this underreporting phenomenon by estimating
#' the number of infections that have occurred but that have not yet been reported. The
#' latter number is then combined with the already reported cases and interpreted as a nowcast or
#' prediction for the true epidemic status regarding the number of daily cases. The routine
#' is anchored around Laplacian-P-splines in an epidemic context (Gressani et al. 2022) and the
#' detailed methodology can be found in Sumalinab et al. (2023). Two different models can be fitted,
#' named M3 (the default) and M2. M3 uses a joint approach that simultaneously models the delay dimension
#' and the time-varying reproduction number. M2 uses reported cases and a nowcast of the not yet reported
#' cases. See Sumalinab et al. (2023) and the vignette \url{https://epilps.com/NowcastingRt.html} for
#' more details.
#'
#' @usage nowcastingR(data, day.effect = TRUE, ref.day = "Monday", si, method = c("M3", "M2"))
#'
#' @param data A data frame containing the data for each time and delay combination with
#' the following 6 columns. The first column is a numeric variable associated to the calendar date.
#' The second column is a numeric variable indicating the delay of reporting. The third column
#' corresponds to the calendar date of the event (e.g. death) and the fourth column to the
#' calendar date at which the event of interest was reported. The fifth column indicates the
#' number of cases for each time and delay combination. Finally, the sixth column indicates whether
#' the cases are already reported or not yet reported. To see an example of such a data structure
#' type \code{data("cov19incidence2022")} and then \code{head(cov19incidence2022)}. This will illustrate the
#' required data structure for nowcasting the reproduction number based on incidence data for
#' Belgium in 2022.
#' @param day.effect If TRUE (default), include the day of the week effect.
#' @param ref.day If \code{day.effect = TRUE}, then the reference category for
#' the day of the week must be specified. The default is "Monday".
#' @param si The (discrete) serial interval distribution.
#' @param method The model to be fitted, either M3 (default) or M2.
#'
#' @return A list with the following components:
#' \itemize{
#'  \item data: The data frame used as an input.
#'  \item Rnow: A data frame containing the nowcasted reproduction number.
#'  \item lambda_estim: Estimated penalty parameters of the P-splines model.
#'  \item phi_estim: Estimated overdispersion parameter from the negative binomial model.
#'  \item method: The model choice, i.e. either M3 or M2.
#'  }
#'
#' @author Bryan Sumalinab (writing) and Oswaldo Gressani (editing).
#'
#' @references Gressani, O., Wallinga, J., Althaus, C. L., Hens, N. and Faes, C.
#'  (2022). EpiLPS: A fast and flexible Bayesian tool for estimation of the
#'  time-varying reproduction number. \emph{Plos Computational Biology},
#'  \strong{18}(10): e1010618.
#' @references Sumalinab, B., Gressani, O., Hens, N. and Faes, C. (2023). An
#'  efficient approach to nowcasting the time-varying reproduction number. MedRxiv preprint.
#'
#' @examples
#' # data("cov19incidence2022")
#' # si_covid <- c(0.344, 0.316, 0.168, 0.104, 0.068) # serial interval distribution
#' # Sys.setlocale("LC_TIME", "English")              # set system locale to English
#' # Rnowfit <- nowcastingR(data = cov19incidence2022, si = si_covid)
#' # tail(Rnowfit$Rnow)
#' # plot(Rnowfit)
#'
#' @export

nowcastingR <- function(data, day.effect = TRUE, ref.day = "Monday", si, method = c("M3", "M2")){
  modelnum <- match.arg(method)
  if (modelnum == "M2") {
    nowcast_inc <- nowcasting(data = data, day.effect = day.effect, verbose = FALSE)
    epifit_nowcast <- estimR(incidence = nowcast_inc$cases.now$y, si = si)
    Rnow.epi <- epifit_nowcast$RLPS
    D <- max(data$d)    # Maximum delay
    TT <- max(data$t)   # Maximum time
    status <- rep("observed", TT)
    status[seq(TT-(D-1),TT)] <- "nowcasted"
    Rnow <- data.frame("Time" = unique(data$Date),
                       "status" = status,
                       Rnow.epi[,-1])
    phi_estim <- nowcast_inc$phi_estim
    lambda_estim <- nowcast_inc$lambda_estim
    method <- method
  } else if (modelnum == "M3") {
    method <- method
    date.start <- as.Date(min(data$Date)) # Starting date
    date.now <- as.Date(max(data$Date))   # Nowcast date
    max.delay <- max(data$d)              # Maximum delay
    data <- data[order(data$d), ]

    # Hyperparameters for Gamma prior of delta
    a.delta <- 1e-05
    b.delta <- 1e-05
    # Prior for overdispersion parameter phi
    a.phi <- 1e-05
    b.phi <- 1e-05
    nu <- 3             # prior parameter for the penalty
    D <- max(data$d)    # Maximum delay
    TT <- max(data$t)   # Maximum time
    nyr <- which(data$Reported == "Not yet reported") # Rows for not yet reported

    # Time and delay
    t <- unique(data$t)
    d <- unique(data$d)

    # Model matrices
    Kt <- 40 # Number of B-splines for time dimension
    Kd <- 10 # Number of B-splines for delay dimension

    # B-spline basis matrix
    Bt <- Rcpp_KercubicBspline(t, lower = min(t), upper = max(t), K = Kt)
    Bd <- Rcpp_KercubicBspline(d, lower = min(d), upper = max(d), K = Kd)
    B <- kronecker(Bd, Bt)  # Two dimensional B-spline matrix
    y <- data$Cases[-nyr]   # Reported cases
    penorder <- 2           # Penalty order

    # Penalty for column (delay dimension)
    Dd <- diag(Kt)
    for (k in 1:penorder) Dd <- diff(Dd)
    Pd <- t(Dd) %*% Dd
    Pd <- Pd + diag(1e-12, Kt)

    # Penalty for row (time dimension)
    Dt <- diag(Kd)
    for (k in 1:penorder) Dt <- diff(Dt)
    Pt <- t(Dt) %*% Dt
    Pt <- Pt + diag(1e-12, Kd)

    if(day.effect == T){
      # Add day of the week
      data$Day <- weekdays(as.Date(data$Rep.date))
      english_days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
      if (!all(data$Day %in% english_days)) {
        stop("Day names are not in English. Set the R language setting to English.")
      }
      data$Day <- stats::relevel(factor(data$Day),ref.day)

      zeta <- 1e-05 # precision week effects coefficient

      # Design matrix for day of the week effect
      Z <- stats::model.matrix(~ Day, data = data)

      # Global design matrix
      # Global design matrix
      X1 <- cbind(B, Z)
      p <- ncol(Z)-1
      X <- X1[-nyr,]    # Design matrix for reported cases

      # Precision matrix for B-spline parameters
      Pv <- function(v) exp(v[1])*(kronecker(diag(1,Kd),Pd)) +
        exp(v[2])*(kronecker(Pt,diag(1,Kt)))

      # Function to create a block-diagonal matrix from matrices A and B
      block_diagonal <- function(A, B) {
        n_A <- nrow(A)
        m_A <- ncol(A)
        n_B <- nrow(B)
        m_B <- ncol(B)
        n <- n_A + n_B
        m <- m_A + m_B
        result <- matrix(0, n, m)
        result[1:n_A, 1:m_A] <- A
        result[(n_A + 1):n, (m_A + 1):m] <- B
        return(result)
      }

      # Precision matrix for parameter xi
      Qv <- function(v) as.matrix(block_diagonal(Pv(v), diag(zeta, p+1)))
    } else{
      # Full design matrix
      X1 <- B
      X_nyr <- X1[nyr,] # Design matrix for not yet reported cases
      X <- X1[-nyr,]    # Design matrix for reported cases

      # Precision matrix for B-spline parameters
      Pv <- function(v) exp(v[1])*(kronecker(diag(1,Kd),Pd)) +
        exp(v[2])*(kronecker(Pt,diag(1,Kt)))
      Qv <- Pv # Precision matrix for parameter xi
    }

    # Negative binomial GLM with log-link
    mu.nb <- function(xi) exp(as.numeric(X %*% xi))
    var.nb <- function(xi, v) {
      muval <- mu.nb(xi)
      res <- muval + (1 / exp(v[3])) * (muval ^ 2)
      return(res)
    }
    W.nb <- function(xi, v) {
      muval <- exp(as.numeric(X %*% xi))
      varval <- muval + (1 / exp(v[3])) * (muval ^ 2)
      res <- diag(((muval) ^ 2) * (1 / varval))
      return(res)
    }
    D.nb <- function(xi) diag(1/mu.nb(xi))
    M.nb <- function(xi) diag(y - mu.nb(xi))
    V.nb <- function(xi, v){
      muval <- mu.nb(xi)
      varval <-  muval + (1 / exp(v[3])) * (muval ^ 2)
      res <- diag(muval * (1/varval - (muval/(varval^2)) *
                             (1 + 2*muval*(1/exp(v[3])))))
      return(res)
    }
    gamma.nb <- function(xi, v) {
      muval <- mu.nb(xi)
      res <- exp(v[3]) * log(muval / (muval + exp(v[3])))
      return(res)
    }
    bgamma.nb <- function(xi, v) - (exp(v[3])^2) *
      log(exp(v[3])/(exp(v[3]) + mu.nb(xi)))

    # Log conditional posterior of xi given v
    log_pxi <- function(xi, v) {
      value <- (1/exp(v[3])) * sum((y * gamma.nb(xi, v)) - bgamma.nb(xi, v)) -
        .5 * t(xi) %*% Qv(v[1:2]) %*% xi
      return(value)
    }

    Grad.logpxi <- function(xi,v){
      muval <- exp(as.numeric(X %*% xi))
      varval <- muval + (1 / exp(v[3])) * (muval ^ 2)
      W.nbval <- diag(((muval) ^ 2) * (1 / varval))
      D.nbval <- diag(1/muval)
      value <- as.numeric(t(X)%*%W.nbval%*%D.nbval%*%(y - muval) -
                            Qv(v[1:2])%*%xi)
      return(value)
    }

    # Hessian of parameter xi
    Hess.logpxi <- function(xi,v){
      value <- t(X)%*%M.nb(xi)%*%V.nb(xi, v)%*%X -
        t(X)%*%W.nb(xi, v)%*%X - Qv(v[1:2])
      value
    }

    # Initial values for log-penalty and log-overdispersion parameter
    v_init <- c(1, 1, 1)
    xi_init <- as.numeric(Rcpp_KerLaplaceNowcast(xi0=rep(0,dim(X)[2]), v = v_init,
                                                 dimxi = ncol(X), Dlogpxi = Grad.logpxi,
                                                 D2logpxi = Hess.logpxi)$xistar)
    M.nbxi_init <- M.nb(xi = xi_init)
    muval_init <- mu.nb(xi_init)

    log_pvcond <- function(vpar){
      v <- c(vpar[1], -3 ,vpar[2])
      varval <-  muval_init + (1 / exp(v[3])) * (muval_init ^ 2)
      Vnb <- diag(muval_init * (1/varval - (muval_init/(varval^2)) *
                                  (1 + 2*muval_init*(1/exp(v[3])))))
      Wnb <- diag(((muval_init) ^ 2) * (1 / varval))
      gammanb <- exp(v[3]) * log(muval_init / (muval_init + exp(v[3])))
      bgammanb <- (-1) * (exp(v[3])^2) * log(exp(v[3])/(exp(v[3]) + muval_init))

      e1 <- eigen(Pv(v[1:2]),only.values = T)$values
      e2 <- eigen(-t(X)%*%(M.nbxi_init%*%Vnb-Wnb)%*%X + Qv(v[1:2]),
                  only.values = T)$values

      value <- sum((1/exp(v[3]))*((y * gammanb) - bgammanb) +
                     lgamma(y + exp(v[3])) - lgamma(exp(v[3]))) +
        0.5*sum(sapply(e1[e1>0], log)) -
        0.5 * sum((xi_init * Qv(v[1:2])) %*% xi_init) +
        a.phi*v[3] - b.phi*exp(v[3]) - 0.5*sum(sapply(e2[e2>0],log)) +
        0.5*nu*(v[1]+v[2]) - (0.5*nu + a.delta)*(log(b.delta + 0.5*nu*exp(v[1]))+
                                                   log(b.delta + 0.5*nu*exp(v[2])))

      return(value)
    }

    vstar <- stats::optim(par = c(1,1), fn = log_pvcond,method = "Nelder-Mead",
                          control = list(fnscale = -1))$par

    # Conditional posterior mode of v
    v_mode <- c(vstar[1],-3, vstar[2])
    xi_mode <- as.numeric(Rcpp_KerLaplaceNowcast(xi0=xi_init, v = v_mode,
                                                 dimxi = ncol(X), Dlogpxi = Grad.logpxi,
                                                 D2logpxi = Hess.logpxi)$xistar)

    if(day.effect == T){
      m <- dim(X1)[2]
      X1[, (m - (p-1)):m] <- 1/7
      Z[, 2:7] <- 1/7
    } else{X1 <- X1}

    # Mean nowcast
    mu_est <- exp(X1%*%xi_mode)
    mean.nowcast <- data
    mean.nowcast[,"Cases"] <- mu_est
    sum.nowcast <- stats::aggregate(Cases ~ t, data = mean.nowcast, FUN = sum)
    colnames(sum.nowcast) <- c("t", "nowcast")

    # R_t estimation
    Rt_LPS <- function(t) {
      if (t == 1) {
        res <- mu_estim[t]
      } else if (t >= 2 && t <= smax) {
        res <-
          mu_estim[t] * ((sum(rev(mu_estim[1:(smax - 1)][1:(t - 1)]) *
                                si[1:(smax - 1)][1:(t - 1)])) ^ (-1))
      } else if (t > smax && t <= n) {
        res <- mu_estim[t] * (sum(rev(mu_estim[(t - smax):(t - 1)]) *
                                    si) ^ (-1)) * (t > smax && t <= n)
      }
      return(res)
    }

    mu_estim <- sum.nowcast$nowcast
    n <- length(mu_estim)
    smax <- length(si)
    R_estim <- sapply(seq_len(n), Rt_LPS)


    #### Credible Interval for Rt

    # Covariance of xi
    sigma_xi_mode.mu <- solve(-Hess.logpxi(xi = xi_mode, v = v_mode))

    CIRt_LPS.delta <- function(t) {
      log_Rt <- log(Rt_LPS(t))

      mustar <- matrix(mu_est, nrow = TT, ncol = D + 1, byrow = FALSE)
      sumexp.bspline <- matrix(NA, nrow = n, ncol = Kt * Kd)
      dmustar.bspline <- t(mustar[t,]) %*% kronecker(t(Bt[t, ]), Bd)

      for (tt in 1:n) {
        Bstar <- kronecker(t(Bt[tt, ]), Bd)
        sumexp.bspline[tt, ] <- t(mustar[tt,]) %*% kronecker(t(Bt[tt, ]), Bd)
      }

      if(day.effect == TRUE){
        sumexp.week <- matrix(NA, nrow = n, ncol = dim(Z)[2])
        for (l in 1:dim(Z)[2]) {
          for (tt in 1:TT) {
            sumexp.week[tt,l] <- mustar[tt, ] %*% matrix(Z[, l], nrow = TT, ncol = D + 1, byrow = TRUE)[tt, ]
          }
        }
        sumexp <- cbind(sumexp.bspline,sumexp.week)
        dmustar.week <- c()
        for (l in 1:dim(Z)[2]) {
          dmustar.week[l] <- sum(mustar[t, ] %*% matrix(Z[, l], nrow = TT, ncol = D + 1, byrow = TRUE)[t, ])
        }
        dmustar <- c(dmustar.bspline, dmustar.week)
        dhstar_t <- (1 / mu_estim[t]) * c(as.vector(matrix(dmustar.bspline, nrow = Kt, ncol = Kd, byrow = TRUE)),dmustar.week)
      } else{
        sumexp <- sumexp.bspline
        dmustar <- dmustar.bspline
        dhstar_t <- (1 / mu_estim[t]) * as.vector(matrix(dmustar, nrow = Kt, ncol = Kd, byrow = TRUE))
      }

      if (t == 1) {
        dxi <- dhstar_t

      } else if (t == 2) {
        dg.xi_jk <- (-1) * sumexp[1:(smax - 1),][1:(t - 1),] * rev(si[1:(smax - 1)][1:(t - 1)])
        if(day.effect == TRUE){
          dg.xi <- c(as.vector(matrix(dg.xi_jk[1:(Kt*Kd)], nrow = Kt, ncol = Kd, byrow = TRUE)),dg.xi_jk[((Kt*Kd)+1):length(dg.xi_jk)])
        } else{
          dg.xi <- as.vector(matrix(dg.xi_jk, nrow = Kt, ncol = Kd, byrow = TRUE))
        }
        g.xi <- (sum(rev(mu_estim[1:(smax - 1)][1:(t - 1)]) * si[1:(smax - 1)][1:(t - 1)])) ^ (-1)
        dxi <- dhstar_t + g.xi * dg.xi

      } else if (t > 2 && t <= smax) {
        dg.xi_jk <- (-1) * colSums(sumexp[1:(smax - 1),][1:(t - 1),] * rev(si[1:(smax - 1)][1:(t - 1)]))

        if(day.effect == TRUE){
          dg.xi <- c(as.vector(matrix(dg.xi_jk[1:(Kt*Kd)], nrow = Kt, ncol = Kd, byrow = TRUE)),dg.xi_jk[((Kt*Kd)+1):length(dg.xi_jk)])
        } else{
          dg.xi <- as.vector(matrix(dg.xi_jk, nrow = Kt, ncol = Kd, byrow = TRUE))
        }
        g.xi <- (sum(rev(mu_estim[1:(smax - 1)][1:(t - 1)]) * si[1:(smax - 1)][1:(t - 1)])) ^ (-1)
        dxi <- dhstar_t + g.xi * dg.xi

      } else if (t > smax && t <= n) {
        dg.xi_jk <- (-1) * colSums(sumexp[(t - smax):(t - 1), ] * rev(si))

        if(day.effect == TRUE){
          dg.xi <- c(as.vector(matrix(dg.xi_jk[1:(Kt*Kd)], nrow = Kt, ncol = Kd, byrow = TRUE)),dg.xi_jk[((Kt*Kd)+1):length(dg.xi_jk)])
        } else{
          dg.xi <- as.vector(matrix(dg.xi_jk, nrow = Kt, ncol = Kd, byrow = TRUE))
        }
        g.xi <- (sum(rev(mu_estim[(t - smax):(t - 1)]) * si)) ^ (-1)
        dxi <- dhstar_t + g.xi * dg.xi
      }

      VlogRt <- as.numeric(t(dxi) %*% sigma_xi_mode.mu %*% dxi)
      sdlogRt <- sqrt(VlogRt)
      Rtsd <- sqrt((exp((sdlogRt^2)) - 1) * exp(2 * log_Rt + sdlogRt^2))
      quantiles <- c(0.025, 0.05, 0.25, 0.50, 0.75, 0.95, 0.975)

      Rq <- function(quantile, log_Rt, VlogRt) {
        qz <- stats::qnorm(quantile, lower.tail = TRUE)
        exp(log_Rt + qz * sqrt(VlogRt))
      }
      Rq_values <- c(Rtsd,lapply(quantiles, Rq, log_Rt = log_Rt, VlogRt = VlogRt))
      names(Rq_values) <- c("Rsd", "Rq0.025", "Rq0.05","Rq0.25",
                            "Rq0.50", "Rq0.75", "Rq0.95", "Rq0.975")
      return(Rq_values)
    }

    status <- rep("observed", TT)
    status[seq(TT-(D-1),TT)] <- "nowcasted"
    Rquant <- data.frame(do.call(cbind, lapply(data.frame(t(sapply(seq_len(n), CIRt_LPS.delta))), unlist)))
    Rnow <- data.frame("Time" = unique(data$Date), "status" = status, "R" = R_estim, Rquant)

    lambda_estim <- data.frame("lambda" = exp(v_mode[1:2]),
                               "description" = c(" lambda_t (penalty for time)",
                                                 "lambda_d (penalty for delay)"))

    phi_estim <- data.frame("phi" = exp(v_mode[3]),
                            "description" = c("overdispersion parameter"))
  } else {
    stop("Method should either be M3 or M2")
  }

  outputlist <- list(data = data,
                     Rnow = Rnow,
                     lambda_estim = lambda_estim,
                     phi_estim = phi_estim,
                     method = modelnum)
  attr(outputlist, "class") <- "Rtnow"
  outputlist
}
