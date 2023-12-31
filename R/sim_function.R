#' Cointegration test for settings of large N and T
#'
#' Runs a simulation on the H0 for the Bykhovskaya-Gorin test for cointegration and returns an empirical p-value. Paper can be found at: https://doi.org/10.48550/arXiv.2202.07150
#' @param N  The number of time series used in simulations.
#' @param tau  The length of the time series used in simulations.
#' @param stat_value The test statistic value for which the p-value is calculated.
#' @param k The number of lags that we wish to employ in the vector autoregression. The default value is k = 1.
#' @param r The number of largest eigenvalues used in the test. The default value is r = 1.
#' @param fin_sample_corr A boolean variable indicating whether we wish to employ finite sample correction on our test statistics. The default value is fin_sample_corr = FALSE.
#' @param sim_num The number of simulations that the function conducts for H0. The default value is sim_num = 1000.
#' @examples
#' result <- sim_function(N=90, tau=501, stat_value=-0.27,k=1,r=1,sim_num=2000)
#' @return A list that contains the simulation values, the empirical percentage (realizations larger than the test statistic provided by the user) and a histogram.
#' @export
#'

sim_function <- function(N=NULL, tau=NULL,stat_value=NULL,k=1,r=1, fin_sample_corr = FALSE, sim_num=1000){

  # Stopping conditions
  error_indicator <- check_input_simfun(N,tau,stat_value,k,r,fin_sample_corr,sim_num)
  if(error_indicator == 1){
    stop()
  }

  print("This function should only be used for quick approximate assessments, as precise computations of the statistics need much larger numbers of simulations.")

  ### this is for progress bar
  options(width = 80)
  n <- sim_num
  ############################

  # extract parameters based on data input
  t = tau-1


  #simulation loop
  stat_vec <- matrix(0,sim_num,1)
  for (j in 1:sim_num){
    X_tilde <- matrix(0, N, tau)
    dX <- matrix(rnorm(N * tau), N, tau)

    for (i in 2:tau){
      if (i == 2){
        X_tilde[,2] <- dX[,1]
      }
      if (i>2){
        X_tilde[, i] <- rowSums(dX[, 1:(i - 1)])
      }
    }

    data_sim <- t(X_tilde) #change to N,T layout because that is the format that largevar() asks for: input is going to be X_tilde

    output <- largevar_scel(data_sim,k,r,fin_sample_corr)
    stat_vec[j,1]<-output

    # #progress report on simulation
    # if(j%%100==0){
    #   base::print(paste("Simulation loop is at",j))
    # }

    ii <- j
    extra <- nchar('||100%')
    width <- options()$width
    step <- round(ii / n * (width - extra))
    text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
                       strrep(' ', width - step - extra), round(ii / n * 100))
    cat(text, " \r")
    flush.console()
   }

  x <- stat_value

  values <- stat_vec[,1]
  percentage <- (length(values[values > x]))/sim_num

  plot <- hist(values,breaks=2*(ceiling(log2(length(values)))+1))
  abline(v=x,col="red",lwd=3)

  list <- list("sim_results"=stat_vec,"empirical_percentage"=percentage,plot)
  new("simfun_output", list)
}
