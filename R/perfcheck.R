#' Check the statistical performance of EpiLPS with simulations
#'
#' @description
#' The \code{perfcheck()} routine can be used to check the performance of EpiLPS
#' in various epidemic scenarios. The user can choose between 4 scenarios,
#' each scenario corresponding to a different data generating process for the
#' incidence data with a specific target dynamics for the reproduction number.
#' The aim of these simulations is to assess how close EpiLPS can reproduce
#' the target reproduction number curve. Different metrics are given as outputs
#' and comparisons with the \code{estimate_R()} routine of the
#' EpiEstim package (Cori et al. 2013) is also shown.
#'
#' @usage perfcheck(S = 10, serial_interval, scenario = 3, K = 30, method = "LPSMAP",
#'           slidewindow = 6, ci_level = 0.95,
#'           themetype = c("classic","gray","light","dark"), seed = 123)
#'
#' @param S The total number of replications.
#' @param scenario The scenario (1,2,3 or 4).
#' @param serial_interval The serial interval distribution.
#' @param method Either LPSMAP (fully sampling-free) or LPSMALA (MCMC-based).
#' @param K Number of (cubic) B-splines in the basis.
#' @param slidewindow The sliding window for EpiEstim (defaults to 1 week).
#' @param ci_level Level of the credible intervals to be computed.
#' @param themetype What theme should be use for plotting the R curves?
#' @param seed A seed for reproducibility.
#'
#' @return An object of class \code{perfcheck} containing a table of summary
#'  statistics for the EpiLPS and EpiEstim routines.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr}
#'
#' @references Cori, A., Ferguson, N.M., Fraser, C., Cauchemez, S. (2013).
#'  A new framework and software to estimate time-varying reproduction numbers
#'  during epidemics. \emph{American Journal of Epidemiology},
#'  \strong{178}(9):1505-1512.
#'
#' @examples
#' simex <- perfcheck(S = 10, serial_interval = c(0.2, 0.4, 0.2, 0.1, 0.1),
#'                    scenario = 3, ci_level = 0.95,  seed = 1234, themetype = "gray")
#'
#' @export

perfcheck <- function(S = 10, serial_interval, scenario = 3, K = 30,
                      method = "LPSMAP", slidewindow = 6, ci_level = 0.95,
                      themetype = c("classic","gray","light","dark"),
                      seed = 123){

  #-- Declare elements to host simulation results
  epidays <- 50
  Repilps <- matrix(0, nrow = S, ncol = epidays)
  colnames(Repilps) <- paste0("Day ", seq_len(epidays))
  Repiestim <- matrix(0, nrow = S, ncol = epidays - (slidewindow + 1))
  colnames(Repiestim) <- paste0("Day ", seq(slidewindow + 2, epidays))
  epilpsCI <- matrix(0, nrow = S, ncol = epidays)
  colnames(epilpsCI) <- paste0("Day ", seq_len(epidays))
  epiestimCI <- matrix(0, nrow = S, ncol = epidays - (slidewindow + 1))
  colnames(epiestimCI) <- paste0("Day ", seq(slidewindow + 2, epidays))
  sim_incid <- matrix(0, nrow = S, ncol = epidays)
  ciwidth_epilps   <- matrix(0, nrow = S, ncol = epidays)
  ciwidth_epiestim <- matrix(0, nrow = S, ncol = epidays - (slidewindow + 1))
  epilps_timing <- c()
  epiestim_timing <- c()
  set.seed(seed)

  #-- Progress bar
  progbar <- progress::progress_bar$new(
    format = crayon::white$yellow("Simulation in progress [:elapsed :spin] [:bar] :percent"),
    total = S,
    clear = FALSE
  )

  #-- Start loop
  for(s in 1:S) {

    #-- Simulate epidemic with episim() routine
    epidemic <- episim(serial_interval = serial_interval, endepi = epidays,
                       Rpattern = scenario)
    sim_incid[s,] <- epidemic$y
    incidence <- epidemic$y
    n <- epidays
    p <- epidemic$serial_interval
    Rtarget <- sapply(seq_len(epidays), epidemic$Rtrue)

    #-- Estimation with epilps
    epilps_fit <- epilps(incidence = incidence, K = K, method = method,
                         serial_interval = p, ci_level = ci_level,
                         verbose = FALSE, progmala = FALSE, tictoc = TRUE)
    epilps_timing[s] <- epilps_fit$elapsed
    ciwidth_epilps[s,] <- (epilps_fit$epifit[, 4] - epilps_fit$epifit[, 3])

    #-- Estimation with epiestim
    t_start <- seq(2, epidays - slidewindow)
    t_end <- t_start + slidewindow
    epiestim_tic <- proc.time()
    epiestim_fit <- suppressMessages(suppressWarnings(
      EpiEstim::estimate_R(incidence, method = "non_parametric_si",
                 config = EpiEstim::make_config(list(si_distr = c(0, p),
                                                t_start = t_start,
                                                t_end = t_end)))))
    epiestim_toc <- proc.time() - epiestim_tic
    epiestim_timing[s] <- round(epiestim_toc[3], 2)
    if (ci_level == 0.95){
      ciwidth_epiestim[s, ] <- epiestim_fit$R$`Quantile.0.975(R)` -
        epiestim_fit$R$`Quantile.0.025(R)`
    } else if (ci_level == 0.90){
      ciwidth_epiestim[s, ] <- epiestim_fit$R$`Quantile.0.95(R)` -
        epiestim_fit$R$`Quantile.0.05(R)`
    }

    # Store point estimation
    Repilps[s, ]   <- epilps_fit$epifit$R_estim
    Repiestim[s, ] <- epiestim_fit$R$`Mean(R)`

    # Check if credible interval contains the truth for epilps
    for (j in 1:epidays) {
      epilpsCI[s, j] <- Rtarget[j] >= epilps_fit$epifit[j, 3] &&
        Rtarget[j] <= epilps_fit$epifit[j, 4]
    }

    # Check if credible interval contains the truth for epiestim
    for (j in 1:length(t_end)) {
      if (ci_level == 0.95) {
        epiestimCI[s, j] <- Rtarget[t_end[j]] >=
          (epiestim_fit$R$`Quantile.0.025(R)`)[j] &&
          Rtarget[t_end[j]] <=
          (epiestim_fit$R$`Quantile.0.975(R)`)[j]
      } else if (ci_level == 0.90) {
        epiestimCI[s, j] <- Rtarget[t_end[j]] >=
          (epiestim_fit$R$`Quantile.0.05(R)`)[j] &&
          Rtarget[t_end[j]] <=
          (epiestim_fit$R$`Quantile.0.95(R)`)[j]
      }
    }
    progbar$tick()
  }

  Repilps <- Repilps[, t_end[1]:epidays]
  epilpsCI <- epilpsCI[, t_end[1]:epidays]
  meanciwidth_epilps   <- colMeans(ciwidth_epilps[, t_end[1]:epidays])
  meanciwidth_epiestim <- colMeans(ciwidth_epiestim)

  #-- Compute metrics starting from day t_end[1] to 50
  Rtruth <- Rtarget[t_end[1]:epidays]

  # Bias
  Bias_epilps <- colMeans(Repilps - matrix(rep(Rtruth,S), nrow = S, byrow = T))
  Bias_epiestim <- colMeans(Repiestim - matrix(rep(Rtruth,S),
                                               nrow = S, byrow = T))
  # MSE
  MSE_epilps <- colMeans((Repilps -
                            matrix(rep(Rtruth,S), nrow = S, byrow = T)) ^ 2)
  MSE_epiestim <- colMeans((Repiestim -
                              matrix(rep(Rtruth,S), nrow = S, byrow = T)) ^ 2)

  #-- Coverage of credible interval at day t
  coverage_epilps   <- round(colMeans(epilpsCI) * 100, 2)
  coverage_epiestim <- round(colMeans(epiestimCI) * 100, 2)

  #-- Summarize performance metrics
  perf_metrics <- matrix(0, nrow = length(t_end), ncol = 6)
  colnames(perf_metrics) <- c("Bias (EpiLPS)", "Bias (EpiEstim)",
                              "MSE (EpiLPS)", "MSE (EpiEstim)",
                              paste0("CP",ci_level * 100,"% (EpiLPS)"),
                              paste0("CP",ci_level * 100,"% (EpiEstim)"))
  rownames(perf_metrics) <- paste0("Day ", t_end)
  perf_metrics[, 1] <- Bias_epilps
  perf_metrics[, 2] <- Bias_epiestim
  perf_metrics[, 3] <- MSE_epilps
  perf_metrics[, 4] <- MSE_epiestim
  perf_metrics[, 5] <- coverage_epilps
  perf_metrics[, 6] <- coverage_epiestim

  simul_summary <- round(perf_metrics[seq(1, length(t_end), by = 1), ], 3)

  #-- Plot results with ggplot2

  themetype <- match.arg(themetype)
  if (themetype == "classic") {
    themeval <- eval(parse(text = "ggplot2::theme_classic()"))
  } else if (themetype == "gray") {
    themeval <- eval(parse(text = "ggplot2::theme_gray()"))
  } else if (themetype == "light") {
    themeval <- eval(parse(text = "ggplot2::theme_light()"))
  } else if (themetype == "dark") {
    themeval <- eval(parse(text = "ggplot2::theme_dark()"))
  }

  # Plot 1: Incidence data
  Days <- rep(seq_len(n), S)
  Incidence <- as.vector(t(sim_incid))
  incid_datframe <- data.frame(Days = Days, Incidence = Incidence)

  inciplot <- ggplot2::ggplot(data = incid_datframe,
                              ggplot2::aes(x = Days, y = Incidence)) +
    ggplot2::geom_point(shape = 8) +
    ggplot2::xlab("Time (days)") +
    themeval +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14)
    )

  # Plot 2: EpiLPS trajectories
  epilpsmedian <- as.numeric(apply(Repilps, 2, "median"))
  epiestimedian <- as.numeric(apply(Repiestim, 2, "median"))
  Repilps_datframe <- suppressWarnings(
    data.frame(Days = rep(seq(t_end[1], epidays), S),
               Rtruth = Rtruth, t(Repilps),
               epilpsmedian = epilpsmedian,
               epiestimedian = epiestimedian))
  datcolnames <- colnames(Repilps_datframe)
  # define colors
  myepilpscol <- grDevices::rgb(0, 211, 255, maxColorValue = 255)
  medepiescol <- grDevices::rgb(212, 0, 52, maxColorValue = 255)
  medepilpscol <- grDevices::rgb(0, 69, 245, maxColorValue = 255)
  # see https://color.adobe.com/create/color-wheel for color definition
  colors <- c("Target R" = "black", "EpiLPS" = myepilpscol,
              "EpiLPS median" = medepilpscol, "EpiEstim median" = medepiescol)
  linetypes <- c("Target R" = 1, "EpiLPS" = 1,
                 "EpiLPS median" = 4, "EpiEstim median" = 3)

  Rlpsplot <- ggplot2::ggplot(data = Repilps_datframe, ggplot2::aes(x = Days)) +
    ggplot2::ylim(0, max(Repilps) + 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = Rtruth, color = "Target R",
                                    linetype = "Target R"), size = 1.2) +
    ggplot2::labs(x = "Time (days)", y = "R", color = "Legend",
                  linetype = "Legend")  +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_linetype_manual(values = linetypes)

  for(s in 3:(S+2)){
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
                       size = 1.1) + themeval +
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
    data.frame(Days = rep(seq(t_end[1], epidays), S),
               Rtruth = Rtruth, t(Repiestim),
               epilpsmedian = epilpsmedian,
               epiestimedian = epiestimedian))
  datcolnames <- colnames(Repiestim_datframe)
  # define colors
  myepiestimcol <- grDevices::rgb(42, 211, 134, maxColorValue = 255)
  # see https://color.adobe.com/create/color-wheel for color definition
  colors <- c("Target R" = "black", "EpiEstim" = myepiestimcol,
              "EpiLPS median" = medepilpscol, "EpiEstim median" = medepiescol)
  linetypes <- c("Target R" = 1, "EpiEstim" = 1,
                 "EpiLPS median" = 4, "EpiEstim median" = 3)

  Repiesplot <- ggplot2::ggplot(data = Repiestim_datframe,
                                ggplot2::aes(x = Days)) +
    ggplot2::ylim(0, max(Repiestim) + 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = Rtruth, color = "Target R",
                                    linetype = "Target R"), size = 1.2) +
    ggplot2::labs(x = "Time (days)", y = "R", color = "Legend",
                  linetype = "Legend")  +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_linetype_manual(values = linetypes)

  for(s in 3:(S+2)){
    Repiesplot <- Repiesplot + eval(parse(
      text = paste("ggplot2::geom_line(ggplot2::aes(y=",datcolnames[s],
                   "),colour='#2AD386')",sep = "")))
  }


  Repiesplot <- Repiesplot +
    ggplot2::geom_line(ggplot2::aes(y = Rtruth), size = 1.2) +
    ggplot2::geom_line(ggplot2::aes(y = epilpsmedian, color = "EpiLPS median",
                                    linetype = "EpiLPS median"), size = 1.1) +
    ggplot2::geom_line(ggplot2::aes(y = epiestimedian,
                                    color = "EpiEstim median",
                                    linetype = "EpiEstim median"),
                       size = 1.1) + themeval +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank()) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 9)
    )


  summary_plot <- gridExtra::grid.arrange(inciplot, Rlpsplot,
                                          Repiesplot, nrow = 1)


  outputlist <- list(simul_summary = simul_summary,
                     plot_summary = summary_plot,
                     ciwidth_epilps = ciwidth_epilps,
                     ciwidth_epiestim = ciwidth_epiestim)

  cat("Comparing ",method," vs EpiEstim in S=",S,
      " replications (epidemic T=50 days). \n", sep ="")
  cat("Mean Bias on days ",t_end[1],"-",epidays, ":\n", sep = "")
  cat("-- EpiLPS mean Bias: ", round(mean(Bias_epilps),5), "\n", sep = "")
  cat("-- EpiEstim mean Bias: ", round(mean(Bias_epiestim, na.rm = TRUE),5),
      "\n", sep = "")
  cat("Mean MSE on days ",t_end[1],"-",epidays, ":\n", sep = "")
  cat("-- EpiLPS mean MSE:   ", round(mean(MSE_epilps),5), "\n", sep = "")
  cat("-- EpiEstim mean MSE: ", round(mean(MSE_epiestim, na.rm = TRUE),5),
      "\n", sep = "")
  cat("Mean credible interval coverage on days ",t_end[1],"-",epidays,
      " (nominal level: ",ci_level * 100," %)",":\n", sep = "")
  cat("-- EpiLPS mean coverage:   ", round(mean(coverage_epilps),5),
      "\n", sep = "")
  cat("-- EpiEstim mean coverage: ", round(mean(coverage_epiestim,
                                                na.rm = TRUE),5),
      "\n", sep = "")
  cat("-- EpiLPS mean CI width: ", round(mean(meanciwidth_epilps),2),
      "\n", sep = "")
  cat("-- EpiEstim mean CI width: ", round(mean(meanciwidth_epiestim,
                                                na.rm = TRUE),2),
      "\n", sep = "")


  attr(outputlist, "class") <- "perfcheck"
  outputlist

}
