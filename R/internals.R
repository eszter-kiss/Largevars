#' Cointegration test for settings of large N and T
#'
#' This is the "skeleton" version of the largevar function in the package. It is called within the sim_function function to make runtime faster. For the actual cointegration test, use the largevar function.
#' @param data a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @return The test statistic.
#' @keywords internal

largevar_scel <- function(data,k=1,r=1){

  #parameters extracted based on data input that we need
  ss = dim(data)
  tau = ss[1]
  t = tau-1
  N =ss[2]

  # create matrix objects to store our transformed data in
  X_tilde=matrix(nrow=N,ncol=t)
  dX=matrix(nrow=N,ncol=t)

  for (i in 1:t){
    X_tilde[,i]=data[i,]-((i-1)/t)*(data[tau,]-data[1,]); # Step 1: De-trend data, time shift
    # and here we change from the (N,T) layout in input to (T,N) layout in our objects
    dX[,i]=(data[i+1,]-data[i,]) #difference the data
  }

  # Based on whether k=1 or k>1 our code is a little different for conducting the test, hence the "if" structure
  if (k==1){
    # Step 2: De-mean data
    R0=matrix(0, N, t)
    R1=matrix(0,N,t)
    meanvec_tilde <- apply(X_tilde, 1 , mean)
    R1 <- apply(X_tilde, 2, function(x) x-meanvec_tilde)
    meanvec_d <- apply(dX, 1 , mean)
    R0 <- apply(dX,2,function(x) x-meanvec_d)

    # Step 3 part 1: Calculate the squared sample canonical correlations between R0 and R1
    S00=R0%*%t(R0)
    S11=R1%*%t(R1)
    S01=R0%*%t(R1)
    S10=R1%*%t(R0)

    # Step 3 part 2: Calculate the eigenvalues of S10 S00^-1 S01 S11^-1 matrix
    can_corr_mat <- S10%*%solve(S00)%*%S01%*%solve(S11)
    ev_values <- eigen(can_corr_mat)$values #the function returns the eigenvalues in descending order in a vector object

    # Step 4: form the test statistic
    loglambda <- log(rep(1,length(ev_values))-ev_values)
    NT <- sum(loglambda[c(1:r)])

    p <- 2
    q <- t/N - 1
    lambda_m <- 1/((p+q)^2)*((sqrt(p*(p+q-1))-sqrt(q))^2)
    lambda_p <- 1/((p+q)^2)*((sqrt(p*(p+q-1))+sqrt(q))^2)
    c_1 <- log(1-lambda_p)
    c_2 <- -((2^(2/3)*lambda_p^(2/3))/(((1-lambda_p)^(1/3))*((lambda_p-lambda_m)^(1/3))))*((p+q)^(-2/3))

    # test statistic
    LR_nt <- (NT-r*c_1)/((N^(-2/3))*c_2)

    return(LR_nt)
  } else { #if k not 1

    #Create cyclic lag matrix
    m <- matrix(1,nrow=t-1,t-1)
    m <- m - lower.tri(m,diag=FALSE)-upper.tri(m,diag = FALSE)
    m <- rbind(rep(0,t-1),m)
    m <- cbind(m,rep(0,t))
    m[1,t] <- 1

    #Create variable matrices for regressions
    Z1=matrix(1,nrow=N*(k-1)+1,ncol=t)
    Zk <- X_tilde

    # Creating X_{t-k} based on VECM form
    for (i in 1:(k-1)){
      Zk <- Zk%*%t(m)
    }

    # Creating the lags with the cyclic lag operator
    cyclic_lag <- dX
    for (j in 1:(k-1)){
      cyclic_lag <- cyclic_lag%*%t(m)
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
    #eigenvalues
    ev_values=eigen(solve(Skk)%*%Sk0%*%solve(S00)%*%S0k)$values

    #test
    loglambda <- log(rep(1,length(ev_values))-ev_values)
    NT <- sum(loglambda[c(1:r)])

    p <- 2
    q <- t/N - k
    lambda_m <- 1/((p+q)^2)*((sqrt(p*(p+q-1))-sqrt(q))^2)
    lambda_p <- 1/((p+q)^2)*((sqrt(p*(p+q-1))+sqrt(q))^2)
    c_1 <- log(1-lambda_p)
    c_2 <- -((2^(2/3)*lambda_p^(2/3))/(((1-lambda_p)^(1/3))*((lambda_p-lambda_m)^(1/3))))*((p+q)^(-2/3))

    #test statistic
    LR_nt <- (NT-r*c_1)/((N^(-2/3))*c_2)

    return(LR_nt)
  }
}

#' Define a class for the largervar() output list
#' # Define the class

setClass("stat_test", representation("list"))

# Define the custom show method
setMethod("show", "stat_test", function(object) {
  cat("Output for the largevars function", "\n")
  cat("===================================", "\n")
  cat("Table of statistics\n")
  cat("\n")  # Print double dashed line
  print(object$significance_test$significance_table)
  cat("\n")
  cat(object$significance_test$`Statistical decision`,"\n")
  cat("============================================================================", "\n")
  cat("Test statistic:", object$statistic ,"\n")
  cat(object$significance_test$text,"\n")
  cat("Decision about the H0: ", object$significance_test$boolean_decision,"\n")
})
