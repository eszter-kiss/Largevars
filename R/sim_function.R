#' Cointegration test for settings of large N and T
#'
#' Runs a simulation on the H0 for the Bykhovskaya-Gorin test for cointegration and returns an empirical p-value. Paper can be found at: https://doi.org/10.48550/arXiv.2202.07150
#' @param data a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @param fin_sample_corr A boolean variable indicating whether we wish to employ finite sample correction on our test statistic. Default is false.
#' @param sim_num The number of simulations that the function conducts for H_0 (the function is slow!). Default is 1000.
#' @examples
#' result <- sim_function(data,k=1,r=1,sim_num=2000)
#' @return A list that contains the simulation values, the empirical percentage (realizations larger than the test statistic for the original data) and a histogram.
#' @export
#'

sim_function <- function(data=NULL,k=1,r=1, fin_sample_corr = FALSE, sim_num=1000){

  # Stopping conditions
  error_indicator <- check_input_simfun(data,k,r,fin_sample_corr, sim_num)
  if(error_indicator == 1){
    stop()
  }

  print("This function should only be used for quick approximate assessments, as precise computations of the statistics need much larger numbers of simulations.")

  # extract parameters based on data input
  ss = dim(data)
  tau = ss[1]
  t = tau-1
  N =ss[2]

  #simulation loop
  stat_vec <- matrix(0,sim_num,1)
  for (j in 1:sim_num){
    X_tilde <- matrix(0, N, tau)
    dX <- matrix(rnorm(N * tau), N, tau)

    for (i in 2:tau){
      if (i == 2){
        X_tilde[,2] <- dX[,1]-1/tau* rowSums(dX[,]) #rowSums does not work when we only have one column so separate case
      }
      if (i>2){
        X_tilde[, i] <- rowSums(dX[, 1:(i - 1)]) - ((i - 1)/tau) * rowSums(dX[,])  # detrend the data and do time shift
      }
    }

    data_sim <- t(X_tilde) #change to N,T layout because that is the format that largevar() asks for: input is going to be X_tilde

    output <- largevar_scel(data_sim,k,r,fin_sample_corr)
    stat_vec[j,1]<-output

    #progress report on simulation
    if(j%%100==0){
      base::print(paste("Simulation loop is at",j))
    }
  }

  x <- largevar_scel(data,k,r,fin_sample_corr)

  values <- stat_vec[,1]
  percentage <- (length(values[values > x]))/sim_num

  plot <- hist(values,breaks=3*(ceiling(log2(length(values)))+1))
  abline(v=x,col="red",lwd=3)

  list <- list("sim_results"=stat_vec,"empirical_percentage"=percentage,plot)
  new("simfun_output", list)
}
