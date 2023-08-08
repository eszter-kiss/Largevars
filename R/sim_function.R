#' Cointegration test for settings of large N and T
#'
#' Runs a simulation on the H0 for the Bykhovskaya-Gorin test for cointegration and returns an empirical p-value. Paper can be found at: https://doi.org/10.48550/arXiv.2202.07150
#' @param data a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @param sim_num The number of simulations that the function conducts for H_0 (the function is slow!). Default is 500.
#' @examples
#' result <- sim_function(data,k=1,r=1,sim_num=20)
#' @return A list that contains the simulation values, the empirical percentage (realizations larger than the test statistic for the original data) and a histogram.
#' @export
#'

sim_function <- function(data,k=1,r=1,sim_num=500){

  # Stopping conditions
  #check if data has been provided
  myArgs <- match.call()
  dataInArgs <- ("data" %in% names(myArgs))
  stopifnot("`data` is a mandatory input" = dataInArgs==TRUE)

  #check correctness of input
  stopifnot("`data` must be a numeric matrix." = is.matrix(data)&is.numeric(data))

  # extract parameters based on data input
  ss = dim(data)
  tau = ss[1]
  t = tau-1
  N =ss[2]

  #check correctness of input
  stopifnot("`k` must be a positive integer." = k%%1==0&k>0)
  stopifnot("`r` must be a positive integer." = r%%1==0&r>0)

  #limit max number of k
  stopifnot("`k` must be such that k+1 < T/N holds." = k<t/N-1)

  #limit max number of r
  stopifnot("`r` must be less than or equal to the number of variables in your dataset." = r<=N)

  #check if sim_num is a positive integer
  stopifnot("`sim_num` must be a positive integer." = sim_num%%1==0&sim_num>0)

  #simulation loop
  stat_vec <- matrix(0,sim_num,1)
  for (j in 1:sim_num){
    X_tilde <- matrix(0, N, tau)
    dX <- matrix(rnorm(N * tau), N, tau)
    # ts.plot(dX[5,])

    for (i in 2:tau){
      if (i == 2){
        X_tilde[,2] <- dX[,1]-1/tau* rowSums(dX[,]) #rowSums does not work when we only have one column so separate case
      }
      if (i>2){
        X_tilde[, i] <- rowSums(dX[, 1:(i - 1)]) - ((i - 1)/tau) * rowSums(dX[,])  # detrend the data and do time shift
      }
    }
    #ts.plot(X_tilde[87,])

    data_sim <- t(X_tilde) #change to N,T layout because that is the format that largevar() asks for: input is going to be X_tilde

    output <- Largevars:::largevar_scel(data_sim,k,r)
    stat_vec[j,1]<-output

    #progress report on simulation
    if(j%%10==0){
      base::print(paste("Simulation loop is at",j))
    }
  }

  x <- Largevars:::largevar_scel(data,k,r)

  values <- stat_vec[,1]
  percentage <- (length(values[values > x]))/sim_num

  plot <- hist(values,breaks=3*(ceiling(log2(length(ev_values)))+1))
  abline(v=x,col="red",lwd=3)

  list <- list("sim_results"=stat_vec,"empirical_percentage"=percentage,plot)
  return(list)
}
