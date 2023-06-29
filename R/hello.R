#' Cointegration test for settings of large N and T
#'
#' Runs the Bykhovskaya-Gorin test for cointegration. Paper can be found at: https://doi.org/10.48550/arXiv.2202.07150
#' @param data a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @param fin_sample_corr A boolean variable indicating whether we wish to employ finite sample correction on our test statistic. Default is false
#' @param plot_output A boolean variable indicating whether we wish to generate the distribution of the eigenvalues (default: TRUE)
#' @examples
#' result <- largevar(data=data,k=1,r=1,plot_output=FALSE);
#' @return A list that contains the test statistic, the upper and lower bounds of the support of the measure of the Wachter distribution, the eigenvalues of our sample matrix, and a table of significance for the test statistic.
#' @export
#'

largevar <- function(data,k=1,r=1, fin_sample_corr = FALSE, plot_output=TRUE){
  #Stopping conditions
  #check if data has been provided (because unlike for the other inputs, this does not have a default)
  myArgs <- match.call()
  dataInArgs <- ("data" %in% names(myArgs))
  stopifnot("`data` is a mandatory input" = dataInArgs==TRUE)

  #parameters extracted from data input that we need
  ss = dim(data) # MATLAB CODE: ss=size(data)
  tau = ss[1] # MATLAB CODE:  tau=ss(1,1);
  t = tau-1 #MATLAB CODE: T=tau-1; %ret_1 is used as X0
  N =ss[2]

  #check correctness of input
  stopifnot("`k` must be a positive integer." = k%%1==0&k>0)
  stopifnot("`r` must be a positive integer." = r%%1==0&r>0)

  #check correctness of input
  stopifnot("`data` must be a numeric matrix." = is.matrix(data)&is.numeric(data))
  stopifnot("`plot_output` must be a boolean." = plot_output==FALSE|plot_output==TRUE)
  stopifnot("`fin_sample_corr` must be a boolean." = fin_sample_corr==FALSE|fin_sample_corr==TRUE)

  #limit max number of k
  stopifnot("`k` must be such that k+1 < T/N holds." = k<t/N-1)

  #limit max number of r
  stopifnot("`r` must be less than or equal to the number of variables in your dataset." = r<=N)

  X_tilde=matrix(nrow=N,ncol=t)
  dX=matrix(nrow=N,ncol=t)

  for (i in 1:t){
    X_tilde[,i]=data[i,]-((i-1)/t)*(data[tau,]-data[1,]); # Step 1, VAR(1) paper: De-trend data, time shift
    # and here we change from the N,T layout in input to T,N layout in our objects
    dX[,i]=(data[i+1,]-data[i,]) #difference the data
  }

  if (k==1){
    # Step 2, VAR(1) paper: De-mean data
    R0=matrix(0, N, t)
    R1=matrix(0,N,t)
    meanvec_tilde <- apply(X_tilde, 1 , mean)
    R1 <- apply(X_tilde, 2, function(x) x-meanvec_tilde)
    meanvec_d <- apply(dX, 1 , mean)
    R0 <- apply(dX,2,function(x) x-meanvec_d);

    # Step 3 part 1, VAR(1) paper: Calculate the squared sample canonical correlations between R0 and R1
    S00=R0%*%t(R0) #in MATLAB CODE (for sp100) THIS IS S11
    S11=R1%*%t(R1) #in MATLAB CODE THIS IS S00
    S01=R0%*%t(R1) #in MATLAB CODE THIS IS S10
    S10=R1%*%t(R0) #in MATLAB CODE THIS IS S01

    # Step 3 part 2, VAR(1) paper: Calculate the eigenvalues of S10 S00^-1 S01 S11^-1 matrix
    can_corr_mat <- S10%*%solve(S00)%*%S01%*%solve(S11)
    ev_values <- eigen(can_corr_mat)$values #the function returns the eigenvalues in descending order in a vector object

    # Step 4: form the test statistic
    loglambda <- log(rep(1,length(ev_values))-ev_values)
    NT <- sum(loglambda[c(1:r)]) # test statistic
    if (fin_sample_corr == FALSE){
      p <- 2
      q <- t/N - 1
    } else{
      p <- 2-2/N
      q <- t/N - 1 - 2/N
    }

    lambda_m <- 1/((p+q)^2)*((sqrt(p*(p+q-1))-sqrt(q))^2)
    lambda_p <- 1/((p+q)^2)*((sqrt(p*(p+q-1))+sqrt(q))^2)
    c_1 <- log(1-lambda_p)
    c_2 <- -((2^(2/3)*lambda_p^(2/3))/(((1-lambda_p)^(1/3))*((lambda_p-lambda_m)^(1/3))))*((p+q)^(-2/3))

    LR_nt <- (NT-r*c_1)/((N^(-2/3))*c_2)

    NT_1 <- sum(loglambda[c(1:1)]) # test statistic
    NT_2 <- sum(loglambda[c(1:2)]) # test statistic
    NT_3 <- sum(loglambda[c(1:3)]) # test statistic

    # tab contains the quantiles for $\sum_{i=1}^r a_i
    # this section appends that table by the test statistic calculated for r=1,2,3

    LR_nt_1 <- (NT_1-1*c_1)/((N^(-2/3))*c_2)
    LR_nt_2 <- (NT_2-2*c_1)/((N^(-2/3))*c_2)
    LR_nt_3 <- (NT_3-3*c_1)/((N^(-2/3))*c_2)

    tests <- c(LR_nt_1,LR_nt_2,LR_nt_3)
    table <- cbind(sum_airy_quantiles, tests)
    colnames(table)[dim(table)[2]] <- "Test statistic"

    decision <- 'If the test statistic is larger than the quantile, we reject the null of no cointegration at the (1-α) level.'
    list_table <- list("significance_table"=table,"Statistical decision"=decision)

    #for r=1 we have more precise statistics so we enter an if loop based on whether r=1 was chosen
    # we round down to two decimals because we have right sided test statistic and two decimals is the precision we have in our table:
    ## if negative: ceiling(abs(a)*100)/100*sign(a)
    ## otherwise: trunc(a*100)/100
    if (r==1){
      # round down to 2 decimals
      airy_1_tab <- as.matrix(airy_1_tab)
      if (sign(LR_nt)==-1){
        val_stat <- (ceiling(abs(LR_nt)*100)/100)*sign(LR_nt)
      }else{#0 or 1
        val_stat <- trunc(LR_nt*100)/100
      }

      coord <- which.min(abs(airy_1_tab[,1] - val_stat)) # this will work fine because we already rounded to two decimals so either it is an exact match, or it is bigger than 3.59 or smaller than -3.89 in which case the ends of the empirical support will be matched to our stat.
      min_coord <- min(abs(airy_1_tab[,1] - val_stat)) #we use this to throw a warning when the statistic matches to the bounds of our empirical support but is way off in distance that p value is actually larger/smaller (similarly to KPSS test in R)
      p_value <- 1-airy_1_tab[coord,2]

      #if loop: if the test statistic value matches to the boundary (-3.89 or 3.59) and its difference is larger than 0.01 (min_coord), we let the user know
      # in practice, enough to check whether min_coord is greater than 0.01 because of our previous rounding
      if (min_coord > 0.01){
        if (sign(val_stat)==-1){
          decision_2 <- paste('The p-value for the test statistic in the r=1 case is greater than ',p_value)
        }
        #
        if (sign(val_stat)==1){
          decision_2 <- paste('The p-value for the test statistic in the r=1 case is less than ',p_value)
        }
      }else{
        decision_2 <- paste('The p-value for the test statistic in the r=1 case is',p_value)
      }

      list_2 <- list("p_value"=p_value, "Precise decision"=decision_2)
      list_table <- append(list_table,list_2)
    }

    list <- list("statistic"=LR_nt, "measure_upper_bound"=lambda_p, "measure_lower_bound"=lambda_m, "eigenvalues"=ev_values,"significance_test"=list_table,"k"=k,"r"=r)

    # Plot the function
    if (plot_output==TRUE){
      my_function <- function(x){(p+q)*sqrt(pmax(0,(lambda_p-x)*(x-lambda_m)))/x/(1-x)/2/pi}
      plot <- hist(ev_values, breaks = 20, probability = TRUE, col = "lightblue", border = "white", main = paste("VAR(",k,") Eigenvalues"),xlab = "Eigenvalues", ylab = "Frequency" ,xlim=c(0,1))
      curve(my_function,add = TRUE)

      list <- append(list,plot)
    }
    return(list)

  } else { #if k not 1
    #Create cyclic lag matrix
    m <- matrix(1,nrow=t-1,t-1)
    m <- m - lower.tri(m,diag=FALSE)-upper.tri(m,diag = FALSE)
    m <- rbind(rep(0,t-1),m)
    m <- cbind(m,rep(0,t))
    m[1,t] <- 1

    #Create variable matrices for regressions

    Z1=matrix(1,nrow=N*(k-1)+1,ncol=t); #vector of regressors (dX_{t-1},...,dX_{t-k+1},dX_{t-k},1)

    Zk <- X_tilde
    for (i in 1:(k-1)){ #shouldnt this go until k?
      Zk <- Zk%*%t(m)
    }

    cyclic_lag <- dX
    for (j in 1:(k-1)){
      cyclic_lag <- cyclic_lag%*%t(m) #creates appropriately lagged dX because m is idempotent so raising it to a power doesn't transform it and R executes the powering first as opposed to multiplying dX by t(m) j amount of times
      Z1[(1+N*(j-1)):(N*j),]<-cyclic_lag
    }

    # stacked regressions
    M11=Z1%*%t(Z1)/t;
    R0=dX-(dX%*%t(Z1)/t)%*%solve(M11)%*%Z1
    Rk=Zk-(Zk%*%t(Z1)/t)%*%solve(M11)%*%Z1
    S00=R0%*%t(R0)
    Skk=Rk%*%t(Rk)
    S0k=R0%*%t(Rk)
    Sk0=Rk%*%t(R0)
    ev_values=eigen(solve(Skk)%*%Sk0%*%solve(S00)%*%S0k)$values

    #test
    loglambda <- log(rep(1,length(ev_values))-ev_values)
    NT <- sum(loglambda[c(1:r)]) # test statistic

    if (fin_sample_corr == FALSE){
      p <- 2
      q <- t/N - k
    } else{
      p <- 2-2/N
      q <- t/N - k - 2/N
    }

    lambda_m <- 1/((p+q)^2)*((sqrt(p*(p+q-1))-sqrt(q))^2)
    lambda_p <- 1/((p+q)^2)*((sqrt(p*(p+q-1))+sqrt(q))^2)
    c_1 <- log(1-lambda_p)
    c_2 <- -((2^(2/3)*lambda_p^(2/3))/(((1-lambda_p)^(1/3))*((lambda_p-lambda_m)^(1/3))))*((p+q)^(-2/3))

    LR_nt <- (NT-r*c_1)/((N^(-2/3))*c_2)

    NT_1 <- sum(loglambda[c(1:1)]) # test statistic
    NT_2 <- sum(loglambda[c(1:2)]) # test statistic
    NT_3 <- sum(loglambda[c(1:3)]) # test statistic

    # tab contains the quantiles for $\sum_{i=1}^r a_i
    # this section appends that table by the test statistic calculated for r=1,2,3

    LR_nt_1 <- (NT_1-1*c_1)/((N^(-2/3))*c_2)
    LR_nt_2 <- (NT_2-2*c_1)/((N^(-2/3))*c_2)
    LR_nt_3 <- (NT_3-3*c_1)/((N^(-2/3))*c_2)

    tests <- c(LR_nt_1,LR_nt_2,LR_nt_3)
    table <- cbind(sum_airy_quantiles, tests)
    colnames(table)[dim(table)[2]] <- "Test statistic"

    decision <- 'If the test statistic is larger than the quantile, we reject the null of no cointegration at the (1-α) level.'
    list_table <- list("significance_table"=table,"Statistical decision"=decision)

    list <- list("statistic"=LR_nt, "measure_upper_bound"=lambda_p, "measure_lower_bound"=lambda_m, "eigenvalues"=ev_values,"significance_test"=list_table,"k"=k,"r"=r)

    # Plot the function
    if (plot_output==TRUE){
      my_function <- function(x){(p+q)*sqrt(pmax(0,(lambda_p-x)*(x-lambda_m)))/x/(1-x)/2/pi}
      plot <- hist(ev_values, breaks = 20, probability = TRUE, col = "lightblue", border = "white", main = paste("VAR(",k,") Eigenvalues"),xlab = "Eigenvalues", ylab = "Frequency" ,xlim=c(0,1))
      curve(my_function,add = TRUE)
      list <- append(list,plot)
    }
    return(list)
  }
}

#' Cointegration test for settings of large N and T
#'
#' Runs a simulation on the H0 for the Bykhovskaya-Gorin test for cointegration. Paper can be found at: https://doi.org/10.48550/arXiv.2202.07150
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
  # Stopping conditions (same exiting conditions as largevar() function in package for applicable inputs)
  #check if data has been provided
  myArgs <- match.call()
  dataInArgs <- ("data" %in% names(myArgs))
  stopifnot("`data` is a mandatory input" = dataInArgs==TRUE)

  ss = dim(data)
  tau = ss[1]
  t = tau-1
  N =ss[2]

  #check correctness of input
  stopifnot("`k` must be a positive integer." = k%%1==0&k>0)
  stopifnot("`r` must be one of the following positive integers: 1, 2, 3." = r==1|r==2|r==3)

  #check correctness of input
  stopifnot("`data` must be a numeric matrix." = is.matrix(data)&is.numeric(data))

  #limit max number of k

  stopifnot("`k` must be less than the number of time observations after differencing." = k<t)

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

    output <- largevars:::largevar(data_sim,k,r,plot_output=FALSE)
    stat_vec[j,1]<-output$statistic

    #progress report on simulation
    if(j%%10==0){
      base::print(paste("Simulation loop is at",j))
    }
  }
  stat_data <- largevars:::largevar(data,k,r,plot_output=FALSE)
  x <- stat_data$statistic

  values <- stat_vec[,1]
  percentage <- length(values[values > x])/sim_num

  plot <- hist(stat_vec[,1],breaks=2*(log2(sim_num)+1))
  abline(v=0.3,col="red",lwd=3)

  list <- list("sim_results"=stat_vec,"empirical_percentage"=percentage,plot)
  return(list)
}
