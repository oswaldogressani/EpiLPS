#' Routine to measure the performance of estimR and estimRmcmc
#'
#' @description
#' This routine can be used to check the 'statistical performance' of the
#' \code{estimR()} and \code{estimRmcmc()} routines to estimate the reproduction
#' number \eqn{R_t}. It simulates epidemics using the \code{episim()} function
#' and computes the Bias, MSE, coverage probability (CP) and width of \eqn{90\%} and
#' \eqn{95\%} credible intervals for \eqn{R_t} averaged over days \eqn{t=8,...,T},
#' where \eqn{T} is the total number of days of the simulated epidemics. As such,
#' it can be used to reproduce part of the results in Gressani et al. (2022)
#' Table 1 and Table 2, respectively. Small differences in results are due to a
#' restructuring of the code since version 1.0.6. If strict reproducible results
#' are required, please refer to version 1.0.6 of the EpiLPS package or visit
#' the GitHub repository \url{https://github.com/oswaldogressani/EpiLPS-ArticleCode}.
#'
#' @usage perfRestim(nsim = 100, scenario = 1, days = 40, K = 40,
#'  method = c("LPSMAP", "LPSMALA"), mcmciter = 3000, burnin = 1000,
#'  si = c("flu", "sars", "mers"), seed = 1325, overdisp = 1000)
#'
#' @param nsim Total number of simulated epidemics.
#' @param scenario The scenario to be used in \code{episim()}.
#' @param days Number of days for the simulated epidemics.
#' @param K Number of B-splines basis function in the P-spline model.
#' @param method The method for LPS, either LPSMAP or LPSMALA.
#' @param mcmciter Number of MCMC samples for method LPSMALA.
#' @param burnin Burn-in for method LPSMALA.
#' @param si The discrete serial interval distribution. Possible specifications
#' are "flu", "sars" or "mers".
#' @param seed A seed for reproducibility.
#' @param overdisp The value of the overdispersion parameter for the negative
#' binomial model in the \code{episim()} routine.
#'
#' @return A list with the following components:
#' \itemize{
#'  \item LPS: Results for the LPS approach.
#'  \item EpiEstim: Results for the EpiEstim approach with weekly sliding
#'  windows.
#'  \item inciplot: The simulated incidence time series.
#'  \item Rlpsplot: Estimated \eqn{R_t} trajectories with LPS.
#'  \item Repiestimplot: Estimated \eqn{R_t} trajectories with EpiEstim.
#' }
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Gressani, O., Wallinga, J., Althaus, C. L., Hens, N. and Faes, C.
#'  (2022). EpiLPS: A fast and flexible Bayesian tool for estimation of the
#'  time-varying reproduction number. \emph{Plos Computational Biology},
#'  \strong{18}(10): e1010618.
#'
#' @examples
#' # # FLU serial interval (Scenarios 1-4)
#' # S1 <- perfRestim(si = "flu", scenario = 1, seed = 1325)
#' # S1mcmc <- perfRestim(si = "flu", scenario = 1, seed = 1325, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S1$inciplot, S1$Rlpsplot, S1$Repiestimplot, nrow = 1))
#' # S2 <- perfRestim(si = "flu", scenario = 2, seed = 1123)
#' # S2mcmc <- perfRestim(si = "flu", scenario = 2, seed = 1123, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S2$inciplot, S2$Rlpsplot, S2$Repiestimplot, nrow = 1))
#' # S3 <- perfRestim(si = "flu", scenario = 3, seed = 1314)
#' # S3mcmc <- perfRestim(si = "flu", scenario = 3, seed = 1314, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S3$inciplot, S3$Rlpsplot, S3$Repiestimplot, nrow = 1))
#' # S4 <- perfRestim(si = "flu", scenario = 4, seed = 1966)
#' # S4mcmc <- perfRestim(si = "flu", scenario = 4, seed = 1966, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S4$inciplot, S4$Rlpsplot, S4$Repiestimplot, nrow = 1))
#' #
#' # # SARS serial interval (Scenarios 5-8)
#' # S5 <- perfRestim(si = "sars", scenario = 1, seed = 1998, overdisp = 5)
#' # S5mcmc <- perfRestim(si = "sars", scenario = 1, seed = 1998, overdisp = 5, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S5$inciplot, S5$Rlpsplot, S5$Repiestimplot, nrow = 1))
#' # S6 <- perfRestim(si = "sars", scenario = 2, seed = 1870, overdisp = 5)
#' # S6mcmc <- perfRestim(si = "sars", scenario = 2, seed = 1870, overdisp = 5, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S6$inciplot, S6$Rlpsplot, S6$Repiestimplot, nrow = 1))
#' # S7 <- perfRestim(si = "sars", scenario = 3, seed = 115,  overdisp = 5)
#' # S7mcmc <- perfRestim(si = "sars", scenario = 3, seed = 115,  overdisp = 5, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S7$inciplot, S7$Rlpsplot, S7$Repiestimplot, nrow = 1))
#' # S8 <- perfRestim(si = "sars", scenario = 4, seed = 1464, overdisp = 5)
#' # S8mcmc <- perfRestim(si = "sars", scenario = 4, seed = 1464, overdisp = 5, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S8$inciplot, S8$Rlpsplot, S8$Repiestimplot, nrow = 1))
#' #
#' # # MERS serial interval (Scenario 9)
#' # S9 <- perfRestim(si = "mers", scenario = 5, days = 60,  seed = 1905, overdisp = 50)
#' # S9mcmc <- perfRestim(si = "mers", scenario = 5, days = 60,
#' # seed = 1905, overdisp = 50, method = "LPSMALA")
#' # suppressWarnings(gridExtra::grid.arrange(S9$inciplot, S9$Rlpsplot,
#' # S9$Repiestimplot, nrow = 1))
#' #
#' # #(Partially recovering Table 2 and Table 3 of Gressani et al. 2022)
#' # simsummary <- matrix(0, nrow = 36, ncol = 7)
#' # colnames(simsummary) <- c("Method", "Bias", "MSE", "CP90%", "CP95%",
#' #                           "CIwidth90%", "CIwidth95%")
#' # simsummary <- as.data.frame(simsummary)
#' #
#' # # Scenario 1
#' # simsummary[1,] <- c(rownames(S1$LPS),S1$LPS)
#' # simsummary[2,] <- c(rownames(S1mcmc$LPS),S1mcmc$LPS)
#' # simsummary[3,] <- c(rownames(S1$EpiEstim),S1$EpiEstim)
#' # simsummary[4,] <- rep("--",7)
#' # # Scenario 2
#' # simsummary[5,] <- c(rownames(S2$LPS),S2$LPS)
#' # simsummary[6,] <- c(rownames(S2mcmc$LPS),S2mcmc$LPS)
#' # simsummary[7,] <- c(rownames(S2$EpiEstim),S2$EpiEstim)
#' # simsummary[8,] <- rep("--",7)
#' # # Scenario 3
#' # simsummary[9,] <- c(rownames(S3$LPS),S3$LPS)
#' # simsummary[10,] <- c(rownames(S3mcmc$LPS),S3mcmc$LPS)
#' # simsummary[11,] <- c(rownames(S3$EpiEstim),S3$EpiEstim)
#' # simsummary[12,] <- rep("--",7)
#' # # Scenario 4
#' # simsummary[13,] <- c(rownames(S4$LPS),S4$LPS)
#' # simsummary[14,] <- c(rownames(S4mcmc$LPS),S4mcmc$LPS)
#' # simsummary[15,] <- c(rownames(S4$EpiEstim),S4$EpiEstim)
#' # simsummary[16,] <- rep("--",7)
#' # # Scenario 5
#' # simsummary[17,] <- c(rownames(S5$LPS),S5$LPS)
#' # simsummary[18,] <- c(rownames(S5mcmc$LPS),S5mcmc$LPS)
#' # simsummary[19,] <- c(rownames(S5$EpiEstim),S5$EpiEstim)
#' # simsummary[20,] <- rep("--",7)
#' # # Scenario 6
#' # simsummary[21,] <- c(rownames(S6$LPS),S6$LPS)
#' # simsummary[22,] <- c(rownames(S6mcmc$LPS),S6mcmc$LPS)
#' # simsummary[23,] <- c(rownames(S6$EpiEstim),S6$EpiEstim)
#' # simsummary[24,] <- rep("--",7)
#' # # Scenario 7
#' # simsummary[25,] <- c(rownames(S7$LPS),S7$LPS)
#' # simsummary[26,] <- c(rownames(S7mcmc$LPS),S7mcmc$LPS)
#' # simsummary[27,] <- c(rownames(S7$EpiEstim),S7$EpiEstim)
#' # simsummary[28,] <- rep("--",7)
#' # # Scenario 8
#' # simsummary[29,] <- c(rownames(S8$LPS),S8$LPS)
#' # simsummary[30,] <- c(rownames(S8mcmc$LPS),S8mcmc$LPS)
#' # simsummary[31,] <- c(rownames(S8$EpiEstim),S8$EpiEstim)
#' # simsummary[32,] <- rep("--",7)
#' # # Scenario 9
#' # simsummary[33,] <- c(rownames(S9$LPS),S9$LPS)
#' # simsummary[34,] <- c(rownames(S9mcmc$LPS),S9mcmc$LPS)
#' # simsummary[35,] <- c(rownames(S9$EpiEstim),S9$EpiEstim)
#' # simsummary[36,] <- rep("--",7)
#' # simsummary
#'
#' @export

perfRestim <- function(nsim = 100, scenario = 1, days = 40, K = 40,
                       method = c("LPSMAP", "LPSMALA"), mcmciter = 3000,
                       burnin = 1000, si = c("flu", "sars", "mers"), seed = 1325,
                       overdisp = 1000){

  set.seed(seed)
  if(match.arg(si) == "flu"){
    scen_num <- scenario
    serial_interval <- c(0.233,0.359,0.198,0.103,0.053,0.027,0.014,0.007,0.003,
                         0.002,0.001)
  } else if(match.arg(si) == "sars"){
    scen_num <- scenario + 4
    serial_interval <- c(0.001,0.012,0.043,0.078,0.104,0.117,0.116,0.108,0.094,
                         0.078,0.063, 0.049,0.038,0.028,0.021,0.015,0.011,
                         0.008,0.005,0.004, 0.003,0.002, 0.001,0.001)
  } else if (match.arg(si) == "mers"){
    scen_num <- 9
    si_mers <- round(EpiEstim::discr_si(k = seq(1,20), mu = 6.8, sigma=4.1),3)
    serial_interval <- si_mers/sum(si_mers)
  }
  if(scenario == 1){
    Rconst <- 1.3
  }

  sim_incid <- matrix(0, nrow = nsim, ncol = days)
  sim_muy <- matrix(0, nrow = nsim, ncol = days)
  RLPS_estim <- matrix(0, nrow = nsim, ncol = days - 7)
  REpiEstim_estim <- matrix(0, nrow = nsim, ncol = days - 7)
  BiasLPS <- matrix(0, nrow = nsim, ncol = days - 7)
  BiasEpiEstim <- matrix(0, nrow = nsim, ncol = days - 7)
  MSELPS <- matrix(0, nrow = nsim, ncol = days - 7)
  MSEEpiEstim <- matrix(0, nrow = nsim, ncol = days - 7)
  CPLPS90 <- matrix(0, nrow = nsim, ncol = days - 7)
  CPEpiEstim90 <- matrix(0, nrow = nsim, ncol = days - 7)
  CPLPS95 <- matrix(0, nrow = nsim, ncol = days - 7)
  CPEpiEstim95 <- matrix(0, nrow = nsim, ncol = days - 7)
  ciwidthLPS90 <- matrix(0, nrow = nsim, ncol = days - 7)
  ciwidthEpiEstim90 <- matrix(0, nrow = nsim, ncol = days - 7)
  ciwidthLPS95 <- matrix(0, nrow = nsim, ncol = days - 7)
  ciwidthEpiEstim95 <- matrix(0, nrow = nsim, ncol = days - 7)
  hyperoptim_convergence <- c()

  for(s in 1:nsim){
    epidemic <- episim(si = serial_interval, endepi = days,
                       Rpattern = scenario, Rconst = Rconst, dist = "negbin",
                       overdisp = overdisp)
    sim_incid[s,] <- epidemic$y
    sim_muy[s,] <- epidemic$mu_y
  }

  Rtarget <- sapply(seq_len(days), epidemic$Rtrue)[-(1:7)]


  cat(paste0("Simulating and fitting ",nsim, " epidemics \n"))
  progbar <- utils::txtProgressBar(min = 1, max = nsim, initial = 1,
                                   style = 3, char ="~")
  for(s in 1:nsim) {

    if(match.arg(method)=="LPSMAP"){
      #-- Estimation with estimR
      epilps_fit <- estimR(incidence = sim_incid[s, ], si = serial_interval,
                           K = K, CoriR = TRUE)
      hyperoptim_convergence[s] <- epilps_fit$optimconverged
    } else if (match.arg(method)=="LPSMALA"){
      #-- Estimation with estimRmcmc
      epilps_fit <- estimRmcmc(incidence = sim_incid[s, ], si = serial_interval,
                               K = K, CoriR = TRUE, niter = mcmciter,
                               burnin = burnin)
      hyperoptim_convergence[s] <- epilps_fit$optimconverged
    }

    #-- Estimation with EpiEstim using 7 days windows
    epiestim_fit <- epilps_fit$RCori

    # EpiLPS results
    RLPS_estim[s, ] <- epilps_fit$RLPS$R[-(1:7)]
    BiasLPS[s, ] <- RLPS_estim[s, ] - Rtarget
    MSELPS[s, ] <- (BiasLPS[s, ] ^ 2)
    if(match.arg(method)=="LPSMAP"){
    CPLPS90[s, ] <- (epilps_fit$RLPS$Rq0.05[-(1:7)] <= Rtarget) &
      (Rtarget <= epilps_fit$RLPS$Rq0.95[-(1:7)])
    CPLPS95[s, ] <- (epilps_fit$RLPS$Rq0.025[-(1:7)] <= Rtarget) &
      (Rtarget <= epilps_fit$RLPS$Rq0.975[-(1:7)])
    ciwidthLPS90[s,] <- (epilps_fit$RLPS$Rq0.95-epilps_fit$RLPS$Rq0.05)[-(1:7)]
    ciwidthLPS95[s,] <- (epilps_fit$RLPS$Rq0.975-epilps_fit$RLPS$Rq0.025)[-(1:7)]
    } else if (match.arg(method) == "LPSMALA"){
      CPLPS90[s, ] <- (epilps_fit$HPD90_Rt[,1][-(1:7)] <= Rtarget) &
        (Rtarget <= epilps_fit$HPD90_Rt[,2][-(1:7)])
      CPLPS95[s, ] <- (epilps_fit$HPD95_Rt[,1][-(1:7)] <= Rtarget) &
        (Rtarget <= epilps_fit$HPD95_Rt[,2][-(1:7)])
      ciwidthLPS90[s,] <- (epilps_fit$HPD90_Rt[,2]-epilps_fit$HPD90_Rt[,1])[-(1:7)]
      ciwidthLPS95[s,] <- (epilps_fit$HPD95_Rt[,2]-epilps_fit$HPD95_Rt[,1])[-(1:7)]
    }

    # EpiEstim results
    REpiEstim_estim[s, ] <- epiestim_fit$`Mean(R)`
    BiasEpiEstim[s, ] <- REpiEstim_estim[s, ] - Rtarget
    MSEEpiEstim[s, ] <- (BiasEpiEstim[s, ] ^ 2)
    CPEpiEstim90[s, ] <- (epiestim_fit$`Quantile.0.05(R)` <= Rtarget) &
      (Rtarget <= epiestim_fit$`Quantile.0.95(R)`)
    CPEpiEstim95[s, ] <- (epiestim_fit$`Quantile.0.025(R)`<= Rtarget) &
      (Rtarget <= epiestim_fit$`Quantile.0.975(R)`)
    ciwidthEpiEstim90[s,] <- (epiestim_fit$`Quantile.0.95(R)`-
                                epiestim_fit$`Quantile.0.05(R)`)
    ciwidthEpiEstim95[s,] <- (epiestim_fit$`Quantile.0.975(R)`-
                                epiestim_fit$`Quantile.0.025(R)`)
    utils::setTxtProgressBar(progbar, s)
  }
  close(progbar)

  # LPS summary
  BiasLPS <- mean(colMeans(BiasLPS))
  MSELPS <- mean(colMeans(MSELPS))

  CPLPS90 <- mean(colMeans(CPLPS90))
  CPLPS95 <- mean(colMeans(CPLPS95))

  meanciwidthLPS90 <- mean(colMeans(ciwidthLPS90))
  meanciwidthLPS95 <- mean(colMeans(ciwidthLPS95))

  summaryLPS <- matrix(0, nrow = 1, ncol = 6)
  colnames(summaryLPS) <- c("Bias","MSE","CP90%","CP95%","CIwidth90%",
                            "CIwidth95%")
  rownames(summaryLPS) <- epilps_fit$method
  summaryLPS[1,] <- round(c(BiasLPS, MSELPS, CPLPS90 * 100, CPLPS95 * 100,
                            meanciwidthLPS90, meanciwidthLPS95), 3)

  # EpiEstim summary
  BiasEpiEstim <- mean(colMeans(BiasEpiEstim), na.rm = TRUE)
  MSEEpiEstim <- mean(colMeans(MSEEpiEstim), na.rm =  TRUE)

  CPEpiEstim90 <- mean(colMeans(CPEpiEstim90), na.rm = TRUE)
  CPEpiEstim95 <- mean(colMeans(CPEpiEstim95), na.rm = TRUE)

  meanciwidthEpiEstim90 <- mean(colMeans(ciwidthEpiEstim90), na.rm = TRUE)
  meanciwidthEpiEstim95 <- mean(colMeans(ciwidthEpiEstim95), na.rm = TRUE)

  summaryEpiEstim <- matrix(0, nrow = 1, ncol = 6)
  colnames(summaryEpiEstim) <- c("Bias","MSE","CP90%","CP95%","CIwidth90%",
                                 "CIwidth95%")
  rownames(summaryEpiEstim) <- "EpiEstim (7d windows)"
  summaryEpiEstim[1,] <- round(c(BiasEpiEstim, MSEEpiEstim, CPEpiEstim90 * 100,
                                 CPEpiEstim95 * 100, meanciwidthEpiEstim90,
                                 meanciwidthEpiEstim95), 3)

  # Plot 1: Incidence data
  Days <- rep(seq_len(days), nsim)
  Incidence <- as.vector(t(sim_incid))
  incid_datframe <- data.frame(Days = Days, Incidence = Incidence)

  inciplot <- ggplot2::ggplot(data = incid_datframe,
                              ggplot2::aes(x = Days, y = Incidence)) +
    ggplot2::geom_point(shape = 8) +
    ggplot2::ggtitle(paste0("Scenario ",scen_num)) +
    ggplot2::xlab("Time (days)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14)
    )

  # Plot 2: EpiLPS trajectories
  colnames(RLPS_estim) <- paste0("Day ", seq(8, days))
  colnames(REpiEstim_estim) <- paste0("Day ", seq(8, days))
  Rtruth <- NULL
  epilpsmedian <- as.numeric(apply(RLPS_estim, 2, "median"))
  epiestimedian <- as.numeric(apply(REpiEstim_estim, 2, "median"))
  Repilps_datframe <- suppressWarnings(
    data.frame(Days = rep(seq(epiestim_fit$t_end[1], days), nsim),
               Rtruth = Rtarget, t(RLPS_estim),
               epilpsmedian = epilpsmedian,
               epiestimedian = epiestimedian))
  datcolnames <- colnames(Repilps_datframe)
  Repilps_datframe <- base::subset(Repilps_datframe, Days > 7)
  myepilpscol <- grDevices::rgb(0, 211, 255, maxColorValue = 255)
  medepiescol <- grDevices::rgb(212, 0, 52, maxColorValue = 255)
  medepilpscol <- grDevices::rgb(0, 69, 245, maxColorValue = 255)
  # see https://color.adobe.com/create/color-wheel for color definition
  colors <- c("Target R" = "black", "EpiLPS" = myepilpscol,
              "EpiLPS median" = medepilpscol, "EpiEstim median" = medepiescol)
  linetypes <- c("Target R" = 1, "EpiLPS" = 1,
                 "EpiLPS median" = 4, "EpiEstim median" = 3)

  Rlpsplot <- ggplot2::ggplot(data = Repilps_datframe, ggplot2::aes(x = Days)) +
    ggplot2::ylim(0,max(rbind(RLPS_estim[,
                        which(colnames(RLPS_estim)=="Day 8"):ncol(RLPS_estim)],
        REpiEstim_estim[,
              which(colnames(REpiEstim_estim)=="Day 8"):ncol(REpiEstim_estim)]),
                      na.rm = TRUE)+ 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = Rtruth, color = "Target R",
                                    linetype = "Target R"), size = 1.2) +
    ggplot2::labs(x = "Time (days)", y = "R", color = "Legend",
                  linetype = "Legend")  +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_linetype_manual(values = linetypes)

  for(s in 3:(nsim+2)){
    Rlpsplot <- Rlpsplot + eval(parse(
      text = paste("ggplot2::geom_line(ggplot2::aes(y=",datcolnames[s],
                   "),colour='#00D3FF')",sep = "")))
  }

  Rlpsplot <- Rlpsplot +
    ggplot2::geom_line(ggplot2::aes(y = Rtruth), size = 1.2) +
    ggplot2::geom_line(ggplot2::aes(y = epilpsmedian, color = "EpiLPS median",
                                    linetype = "EpiLPS median"), size = 1.1) +
    ggplot2::geom_line(ggplot2::aes(y = epiestimedian,
                                    color = "EpiEstim median",
                                    linetype = "EpiEstim median"),
                       size = 1.1) +
    ggplot2::ggtitle(paste0(method, " trajectories")) +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 9)
    )

  # Plot 3: EpiEstim trajectories
  Repiestim_datframe <- suppressWarnings(
    data.frame(Days = rep(seq(epiestim_fit$t_end[1], days), nsim),
               Rtruth = Rtarget, t(REpiEstim_estim),
               epilpsmedian = epilpsmedian,
               epiestimedian = epiestimedian))
  datcolnames <- colnames(Repiestim_datframe)
  Repiestim_datframe <- base::subset(Repiestim_datframe, Days > 7)
  myepiestimcol <- grDevices::rgb(42, 211, 134, maxColorValue = 255)
  colors <- c("Target R" = "black", "EpiEstim" = myepiestimcol,
              "EpiLPS median" = medepilpscol, "EpiEstim median" = medepiescol)
  linetypes <- c("Target R" = 1, "EpiEstim" = 1,
                 "EpiLPS median" = 4, "EpiEstim median" = 3)

  Repiestimplot <- ggplot2::ggplot(data = Repiestim_datframe,
                                   ggplot2::aes(x = Days)) +
    ggplot2::ylim(0,max(rbind(RLPS_estim[,
              which(colnames(RLPS_estim)=="Day 8"):ncol(RLPS_estim)],
                            REpiEstim_estim[,
              which(colnames(REpiEstim_estim)=="Day 8"):ncol(REpiEstim_estim)]),
                      na.rm = TRUE)+ 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = Rtruth, color = "Target R",
                                    linetype = "Target R"), size = 1.2) +
    ggplot2::labs(x = "Time (days)", y = "R", color = "Legend",
                  linetype = "Legend")  +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_linetype_manual(values = linetypes)

  for(s in 3:(nsim+2)){
    Repiestimplot <- Repiestimplot + eval(parse(
      text = paste("ggplot2::geom_line(ggplot2::aes(y=",datcolnames[s],
                   "),colour='#2AD386')",sep = "")))
  }

  Repiestimplot <- Repiestimplot +
    ggplot2::geom_line(ggplot2::aes(y = Rtruth), size = 1.2) +
    ggplot2::geom_line(ggplot2::aes(y = epilpsmedian, color = "EpiLPS median",
                                    linetype = "EpiLPS median"), size = 1.1) +
    ggplot2::geom_line(ggplot2::aes(y = epiestimedian,
                                    color = "EpiEstim median",
                                    linetype = "EpiEstim median"),
                       size = 1.1) +
    ggplot2::ggtitle("EpiEstim 7d windows trajectories") +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 9)
    )

  convergence_message <- paste0("Algorithm for hyperparameter optimization converged for ",
                                sum(hyperoptim_convergence)/
                                  nsim * 100,"% of the simulated epidemics.")
  outlist <- list(LPS = summaryLPS, EpiEstim = summaryEpiEstim,
                  convergence_message = convergence_message,
                  inciplot = inciplot, Rlpsplot = Rlpsplot,
                  Repiestimplot = Repiestimplot)

  return(outlist)
}

