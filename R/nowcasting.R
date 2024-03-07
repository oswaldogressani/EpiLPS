#' Nowcasting and estimation of occurred but not yet reported events
#'
#' @description
#' This routine can be used to estimate cases that have not yet been reported also known as
#' occurred-but-not-yet-reported-events. Daily cases are typically subject to reporting
#' delays, so that the reported number of infected individuals is not always reflecting the
#' true epidemic status. Nowcasting aims to correct this underreporting phenomenon by estimating
#' the number of infections that have occurred but that have not yet been reported. The
#' latter number is then combined with the already reported cases and interpreted as a nowcast or
#' prediction for the true epidemic status regarding the number of daily cases. The routine
#' is anchored around Laplacian-P-splines in an epidemic context (Gressani et al. 2022) and the
#' detailed methodology can be found in Sumalinab et al. (2023).
#'
#' @param data A data frame containing the data for each time and delay combination with
#' the following 6 columns. The first column is a numeric variable associated to the calendar date.
#' The second column is a numeric variable indicating the delay of reporting. The third column
#' corresponds to the calendar date of the event (e.g. death) and the fourth column to the
#' calendar date at which the event of interest was reported. The fifth column indicates the
#' number of cases for each time and delay combination. Finally, the sixth column indicates whether
#' the cases are already reported or not yet reported. To see an example of such a data structure
#' type \code{data("cov19mort2021")} and then \code{head(cov19mort2021)}. This will illustrate the
#' data structure used for nowcasting the daily number of cases based on mortality data for Belgium
#' in 2021.
#' @param day.effect If TRUE (defaut), include the day of the week effect.
#' @param ref.day If \code{day.effect = TRUE}, then the reference category for
#' the day of the week must be specified. The default is "Monday".
#' @param verbose Show summary output of nowcast in console? Default is TRUE.
#'
#' @return A list with the following components:
#' \itemize{
#'  \item data: The data frame used as an input.
#'  \item cases.now: A data frame containing the nowcasting results with \eqn{95\%} prediction intervals.
#'  \item delay: A data frame containing the two-dimensional delay distribution.
#'  \item lambda_estim: Estimated penalty parameters of the P-splines model.
#'  \item phi_estim: Estimated overdispersion parameter from the negative binomial model.
#'  }
#'
#' @author Bryan Sumalinab (writing) and Oswaldo Gressani (editing).
#'
#' @references Gressani, O., Wallinga, J., Althaus, C. L., Hens, N. and Faes, C.
#'  (2022). EpiLPS: A fast and flexible Bayesian tool for estimation of the
#'  time-varying reproduction number. \emph{Plos Computational Biology},
#'  \strong{18}(10): e1010618.
#' @references Sumalinab, B., Gressani, O., Hens, N. and Faes, C. (2023). Bayesian
#'  nowcasting with Laplacian-P-splines. MedRxiv preprint.
#'  \url{https://www.medrxiv.org/content/10.1101/2022.08.26.22279249v2}
#'
#' @examples
#' # data("cov19mort2021")
#' # ncast <- nowcasting(data = cov19mort2021, day.effect = FALSE)
#' # plot(ncast) # Show nowcasted cases
#' # plot(ncast, type = "delay") # Show contour of delay distribution
#'
#' @export


nowcasting <- function(data, day.effect = TRUE, ref.day = "Monday", verbose = TRUE){
  tic <- proc.time()
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
    english_days <- c("Sunday", "Monday", "Tuesday", "Wednesday",
                      "Thursday", "Friday", "Saturday")
    if (!all(data$Day %in% english_days)) {
      stop("Day names are not in English. Set R language settings to English.")
    }
    data$Day <- stats::relevel(factor(data$Day),ref.day)

    zeta <- 1e-05 # precision week effects coefficient

    # Design matrix for day of the week effect
    Z <- stats::model.matrix(~ Day, data = data)

    # Global design matrix
    X1 <- cbind(B, Z)
    p <- ncol(Z)-1
    X_nyr <- X1[nyr,] # Design matrix for not yet reported cases
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
                 control = list(fnscale = -1, reltol = 1e-3))$par

  # Conditional posterior mode of v
  v_mode <- c(vstar[1],-3, vstar[2])
  xi_mode <- as.numeric(Rcpp_KerLaplaceNowcast(xi0=xi_init, v = v_mode,
                                        dimxi = ncol(X), Dlogpxi = Grad.logpxi,
                                        D2logpxi = Hess.logpxi)$xistar)

  # Nowcast for not yet reported
  mu_nyr <- exp(X_nyr%*%xi_mode)
  nowcast <- data
  nowcast[nyr,"Cases"] <- mu_nyr

  # Nowcasted cases (reported + nowcast)
  cases.now <- stats::aggregate(Cases ~ t, data = nowcast,
                         FUN = function(x) ceiling(sum(x)))
  colnames(cases.now) <- c("t", "y")
  cases.now$Date <- unique(nowcast$Date)
  cases.now$CI95L <- NA
  cases.now$CI95R <- NA

  ##### Prediction Interval
  # Time that has not yet reported cases is t = T-(D-1),...,T.
  Date <- data$Date
  Cases <- data$Cases
  Reported <- data$Reported

  t.now <- max(data$t)-(D-1)
  nowcast.nyr <- subset(nowcast, t >= t.now, select = c(Date, Cases))

  # Summarize the data and nowcast results to be used for summarizing prediction
  # interval in next lines of code
  data1 <- stats::aggregate(Cases ~ Date + Reported, data = data, FUN = sum)
  data1$Reported <- factor(data1$Reported,
                           levels = c("Reported", "Not yet reported",
                                      "Nowcast"),
                           labels = c("Reported", "Not yet reported",
                                      "Nowcast"))
  data1$Date <- as.Date(data1$Date)
  data1 <- as.data.frame(data1)

  data2 <- stats::aggregate(Cases ~ Date, data = nowcast.nyr, FUN = sum)
  data2$Date <- as.Date(data2$Date)
  data2 <- as.data.frame(data2)

  # Covariance of xi
  sigma_xi <- -solve(Hess.logpxi(xi = xi_mode, v = v_mode))

  logmu <- c()
  logmu.var <- c()
  for(i in 1:dim(X_nyr)[1]){
    logmu[i] <- X_nyr[i,]%*%xi_mode
    logmu.var[i] <- X_nyr[i,]%*%sigma_xi%*%X_nyr[i,]
  }

  # Generate negative binomial samples
  r.nb <- list()
  N <- 1000
  for (i in 1:length(logmu)) {
    rn <- stats::rnorm(N, mean = logmu[i], sd = sqrt(logmu.var[i]))
    mu <- exp(rn)
    r.nb[[i]] <- stats::rnbinom(n = N, mu = mu, size = exp(v_mode[3]))
  }

  # Data for not yet reported cases
  data_nyr <- data[nyr,]
  data_nyr$nowcast <- exp(logmu)
  data_CI <- stats::aggregate(nowcast ~ t, data = data_nyr, FUN = sum)

  T.now <- as.numeric(date.now - date.start) + 1
  days <- seq(T.now-(D-1),T.now) # dates that have nowcast

  # Sum of negative binomial samples for each t (days)
  for (i in 1:dim(data_CI)[1]) {
    r.nb.sum <- Reduce("+",r.nb[which(data_nyr$t == days[i])])
    data_CI[i,"lower_nyr"] <- stats::quantile(r.nb.sum,probs = 0.025)
    data_CI[i,"upper_nyr"] <- stats::quantile(r.nb.sum,probs = 0.975)
  }

  CI_nyr <- data.frame(data2,data_CI[,c("t","nowcast","lower_nyr","upper_nyr")])
  data_rep_cases <- subset(data1, Date %in% CI_nyr$Date & Reported == "Reported")
  data_rep_cases <- as.data.frame(data_rep_cases)
  CI_nyr$Rep_Cases <- data_rep_cases$Cases

  CI_nyr$CI95L <- ceiling(CI_nyr$lower_nyr + CI_nyr$Rep_Cases)
  CI_nyr$CI95R <- ceiling(CI_nyr$upper_nyr + CI_nyr$Rep_Cases)

  cases.now[(nrow(cases.now) - max.delay + 1):nrow(cases.now),
            c("CI95L", "CI95R")] <- CI_nyr[, c("CI95L", "CI95R")]
  cases.now$status <- ifelse(is.na(cases.now$CI95L), "observed", "nowcasted")

  # Delay density
  mu_hat <- exp(X1[,1:dim(B)[2]]%*%xi_mode[1:dim(B)[2]])
  nowcast2 <- data
  nowcast2[,"Cases"] <- mu_hat

  cases_matrix <- as.data.frame(matrix(nowcast2$Cases,nrow = TT, ncol = D+1,
                                       byrow = F))
  delaydist <- t(apply(cases_matrix, MARGIN = 1, function(i) i/sum(i)))
  data_delay <- data.frame("Date" = as.Date(data$Date),
                           "Delay" = data$d,"density" = as.numeric(delaydist))

  lambda_estim <- data.frame("lambda" = exp(v_mode[1:2]),
                             "description" = c(" lambda_t (penalty for time)",
                                               "lambda_d (penalty for delay)"))

  phi_estim <- data.frame("phi" = exp(v_mode[3]),
                          "description" = c("overdispersion parameter"))

  toc <- proc.time() - tic

  if(isTRUE(verbose)){
    outprint <- cases.now[cases.now$status=="nowcasted",]
    outprint2 <- matrix(0, nrow = nrow(outprint), ncol = ncol(outprint)-1)
    colnames(outprint2) <- c("Time", "Date", "Nowcast", "CI95L", "CI95R")
    outprint2 <- as.data.frame(outprint2)
    outprint2[,1] <- outprint$t
    outprint2[,2] <- as.Date(outprint$Date)
    outprint2[,3] <- outprint$y
    outprint2[,4] <- outprint$CI95L
    outprint2[,5] <- outprint$CI95R
    cat("Nowcast results: \n")
    cat("-------------------------------------------------------\n")
    print(outprint2)
    cat("-------------------------------------------------------\n")
    cat("Time elapsed: ", paste0(round(toc[3],3),"s"))
  }

  outputlist <- list(data = data,
                 cases.now = cases.now,
                 delay = data_delay,
                 lambda_estim = lambda_estim,
                 phi_estim = phi_estim)

  attr(outputlist, "class") <- "nowcasted"
  outputlist
}
