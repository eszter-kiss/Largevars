#' Cointegration test for settings of large N and T
#'
#' Runs the Bykhovskaya-Gorin test for cointegration. Paper can be found at: https://doi.org/10.48550/arXiv.2202.07150
#' @param data a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @param fin_sample_corr A boolean variable indicating whether we wish to employ finite sample correction on our test statistic. Default is false
#' @param plot_output A boolean variable indicating whether we wish to generate the distribution of the eigenvalues (default: TRUE)
#' @param significance_level Specify the significance level at which the decision about the H0 should be made. For r=1 this can be any level of significance. For r=2 and r=3, the significance level input will be rounded up to the nearest of the following: 0.1, 0.05, 0.025, 0.01. If the significance level is larger than 0.1, the decision will be made at the 10% level. For r>3 only the test statistic is returned. For an empirical p-value for r>3 use the sim_function fun. in the package.
#' @examples
#' result <- largevar(data=data,k=1,r=1,fin_sample_corr=FALSE, plot_output=FALSE,significance_level=0.05);
#' @return A list that contains the test statistic, the upper and lower bounds of the support of the measure of the Wachter distribution, the eigenvalues of our sample matrix, and a table of significance for the test statistic.
#' @export
#'

largevar <- function(data,k=1,r=1, fin_sample_corr = FALSE, plot_output=TRUE, significance_level = 0.05){

#Stopping conditions
  #check if data has been provided (because this input does not have a default)
  myArgs <- match.call()
  dataInArgs <- ("data" %in% names(myArgs))
  stopifnot("`data` is a mandatory input" = dataInArgs==TRUE)

  #check correctness of input
  stopifnot("`data` must be a numeric matrix." = is.matrix(data)&is.numeric(data))

  #parameters extracted based on data input that we need
  ss = dim(data)
  tau = ss[1]
  t = tau-1
  N =ss[2]

  #check correctness of input
  stopifnot("`k` must be a positive integer." = k%%1==0&k>0)
  stopifnot("`r` must be a positive integer." = r%%1==0&r>0)
  stopifnot("`significance_level` must be a real number strictly between 0 and 1." = significance_level<1&significance_level>0)
  stopifnot("`plot_output` must be a boolean." = plot_output==FALSE|plot_output==TRUE)
  stopifnot("`fin_sample_corr` must be a boolean." = fin_sample_corr==FALSE|fin_sample_corr==TRUE)

  #limit max number of k
  stopifnot("`k` must be such that k+1 < T/N holds." = k<t/N-1)

  #limit max number of r
  stopifnot("`r` must be less than or equal to the number of variables in your dataset." = r<=N)

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

      #based on whether we have finite sample correction request by the user, the parameters p and q are corrected by a 2/N factor in this "if" function
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

      # test statistic
      LR_nt <- (NT-r*c_1)/((N^(-2/3))*c_2)

      # these are calculated in addition to the above test statistic based on the user's input in order to merge them into the statistical table in the output to help the user make a decision on H0

      NT_1 <- sum(loglambda[c(1:1)])
      LR_nt_1 <- (NT_1-1*c_1)/((N^(-2/3))*c_2) # test statistic for r=1
      NT_2 <- sum(loglambda[c(1:2)])
      LR_nt_2 <- (NT_2-2*c_1)/((N^(-2/3))*c_2) # test statistic for r=2
      NT_3 <- sum(loglambda[c(1:3)])
      LR_nt_3 <- (NT_3-3*c_1)/((N^(-2/3))*c_2) # test statistic for r=3

      # this section appends the statistical table for r=1,2,3 by the test statistics calculated for r=1,2,3
      tests <- c(LR_nt_1,LR_nt_2,LR_nt_3)
      table <- cbind(sum_airy_quantiles[,1:4], tests)
      colnames(table)[dim(table)[2]] <- "Test statistic"

      #guide for how to make the decision on H0
      decision <- 'If the test statistic is larger than the quantile, we reject the null of no cointegration at the α level.'
      # makes intermediary list to merge in to the final list that is the output of the function
      list_table <- list("significance_table"=table,"Statistical decision"=decision)

    #for r=1 we have more precise statistics so we enter an if loop based on whether r=1 was chosen
    # we round down to two decimals because we have right sided test statistic and two decimals is the precision we have in our table (airy_1_tab)
    #rounding is as follows:
      ## if negative: ceiling(abs(a)*100)/100*sign(a)
      ## otherwise: trunc(a*100)/100
    if (r==1){
      # transforms table in package into matrix format
      airy_1_tab <- as.matrix(airy_1_tab)

      if (sign(LR_nt)==-1){
        val_stat <- (ceiling(abs(LR_nt)*100)/100)*sign(LR_nt)
      }else{#0 or 1
        val_stat <- trunc(LR_nt*100)/100
      }

      coord <- which.min(abs(airy_1_tab[,1] - val_stat)) # this will work fine because we already rounded to two decimals so either it is an exact match, or it is bigger than 3.59 or smaller than -3.89 in which case the ends of the empirical support will be matched to our stat.
      min_coord <- min(abs(airy_1_tab[,1] - val_stat)) # we use this to know when the statistic matches to the bounds of our empirical support but is actually far away from it so the p-value is actually larger/smaller (similarly to KPSS test in R) than what the output will indicate
      p_value <- 1-airy_1_tab[coord,2]

      #if loop: if the test statistic value matches to the boundary (-3.89 or 3.59) and its difference is larger than 0.01 (min_coord), we let the user know in the decision_2 text
      # enough to check whether min_coord is greater than 0.01 because of our previous rounding to two decimals
      if (min_coord > 0.01){
        if (sign(val_stat)==-1){ # if the sign of the test stat. is negative, it matched to the left boundary
          decision_2 <- paste('The p-value for the test statistic in the r=1 case is greater than ',p_value)
          # Accept (0) or reject (1) based on variable `significance_level`

          decision_3 <- 0 # if the statistic is so small that it matches to the left boundary, the p-value of that statistic is giant, we accept
        }
        #
        if (sign(val_stat)==1){ # if the sign of the test stat. is positive, it matched to the right boundary
          decision_2 <- paste('The p-value for the test statistic in the r=1 case is less than ',p_value)
          # Accept (0) or reject (1) based on variable `significance_level`
          # if significance level specified by the user is smaller than the smallest p-value that we know the statistic for based on the airy_1_tab then technically we don't know what the decision would be at the very small level of significance we specified so we return the decision based on the smallest p-value

          decision_3 <- 1 # if the statistic is so big that it matches to the right boundary, the p-value of that statistic is tiny, we reject
        }
      }else{
        decision_2 <- paste('The p-value for the test statistic in the r=1 case is',p_value)
        # Accept (0) or reject (1) based on variable `significance_level`
        if (p_value >= significance_level){
          decision_3 <- 0
        }else{
          decision_3 <- 1
        }

      }
      #intermediary list to append our final, output list with
      list_2 <- list("p_value"=p_value, "Precise decision"=decision_2, "boolean_decision" = decision_3)
      list_table <- append(list_table,list_2)
    }


      if (r==2|r==3){ #we still have a significance table we can go by when r=2 or r=3, hence this loop

        if (significance_level > 0.05){ #round the significance level according to the description document (significance_level var. description in ?largevar)
          round_s_f <- 0.1
        }else{
          if (significance_level > 0.025){
            round_s_f <- 0.05
          }else{
            if (significance_level > 0.01){
              round_s_f <- 0.025
            }else{
              round_s_f <- 0.01
            }
          }
        }

        if (r==2){

          if (round_s_f == 0.1){
            if (sum_airy_quantiles[2,1] >= LR_nt){
              decision_3 <- 0
            }else{
              decision_3 <- 1
            }
          }

          if (round_s_f == 0.05){
            if (sum_airy_quantiles[2,2] >= LR_nt){
              decision_3 <- 0
            }else{
              decision_3 <- 1
            }
          }

          if (round_s_f == 0.025){
            if (sum_airy_quantiles[2,3] >= LR_nt){
              decision_3 <- 0
            }else{
              decision_3 <- 1
            }
          }

          if (round_s_f == 0.01){
            if (sum_airy_quantiles[2,4] >= LR_nt){
              decision_3 <- 0
            }else{
              decision_3 <- 1
            }
          }


        }else{#r==3

          if (round_s_f == 0.1){
            if (sum_airy_quantiles[3,1] >= LR_nt){
              decision_3 <- 0
            }else{
              decision_3 <- 1
            }
          }

          if (round_s_f == 0.05){
            if (sum_airy_quantiles[3,2] >= LR_nt){
              decision_3 <- 0
            }else{
              decision_3 <- 1
            }
          }

          if (round_s_f == 0.025){
            if (sum_airy_quantiles[3,3] >= LR_nt){
              decision_3 <- 0
            }else{
              decision_3 <- 1
            }
          }

          if (round_s_f == 0.01){
            if (sum_airy_quantiles[3,4] >= LR_nt){
              decision_3 <- 0
            }else{
              decision_3 <- 1
            }
          }
        }

        list_2 <- list("significance_level"=round_s_f,"boolean_decision" = decision_3)
        list_table <- append(list_table,list_2)
      }

    list <- list("statistic"=LR_nt, "measure_upper_bound"=lambda_p, "measure_lower_bound"=lambda_m, "eigenvalues"=ev_values,"significance_test"=list_table,"k"=k,"r"=r)

    # Plot the function
    if (plot_output==TRUE){
      my_function <- function(x){(p+q)*sqrt(pmax(0,(lambda_p-x)*(x-lambda_m)))/x/(1-x)/2/pi}
      plot <- hist(ev_values, breaks = 2.5*(log2(length(ev_values))+1), probability = TRUE, col = "lightblue", border = "white", main = paste("VAR(",k,") Eigenvalues"),xlab = "Eigenvalues", ylab = "Frequency" ,xlim=c(0,1))
      curve(my_function,add = TRUE)

      # this will guarantee that if plot_output=TRUE by user then the histogram plot will automatically pop up in "Plots" after running function
      list <- append(list,plot)
    }

    # function output
    return(list)

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
    NT <- sum(loglambda[c(1:r)]) # test statistic

    #finite sample correction in case user specified it, otherwise no correction for p and q
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

    #test statistic
    LR_nt <- (NT-r*c_1)/((N^(-2/3))*c_2)

    NT_1 <- sum(loglambda[c(1:1)])
    LR_nt_1 <- (NT_1-1*c_1)/((N^(-2/3))*c_2) # test statistic
    NT_2 <- sum(loglambda[c(1:2)])
    LR_nt_2 <- (NT_2-2*c_1)/((N^(-2/3))*c_2) # test statistic
    NT_3 <- sum(loglambda[c(1:3)])
    LR_nt_3 <- (NT_3-3*c_1)/((N^(-2/3))*c_2) # test statistic

    # appends the statistical table with the above calculated test statistics
    tests <- c(LR_nt_1,LR_nt_2,LR_nt_3)
    table <- cbind(sum_airy_quantiles[,1:4], tests)
    colnames(table)[dim(table)[2]] <- "Test statistic"

    # decision texts
    decision <- 'If the test statistic is larger than the quantile, we reject the null of no cointegration at the α level.'
    list_table <- list("significance_table"=table,"Statistical decision"=decision)

    if (r==1){
      # transforms table in package into matrix format
      airy_1_tab <- as.matrix(airy_1_tab)

      if (sign(LR_nt)==-1){
        val_stat <- (ceiling(abs(LR_nt)*100)/100)*sign(LR_nt)
      }else{#0 or 1
        val_stat <- trunc(LR_nt*100)/100
      }

      coord <- which.min(abs(airy_1_tab[,1] - val_stat)) # this will work fine because we already rounded to two decimals so either it is an exact match, or it is bigger than 3.59 or smaller than -3.89 in which case the ends of the empirical support will be matched to our stat.
      min_coord <- min(abs(airy_1_tab[,1] - val_stat)) # we use this to know when the statistic matches to the bounds of our empirical support but is actually far away from it so the p-value is actually larger/smaller (similarly to KPSS test in R) than what the output will indicate
      p_value <- 1-airy_1_tab[coord,2]

      #if loop: if the test statistic value matches to the boundary (-3.89 or 3.59) and its difference is larger than 0.01 (min_coord), we let the user know in the decision_2 text
      # enough to check whether min_coord is greater than 0.01 because of our previous rounding to two decimals
      if (min_coord > 0.01){
        if (sign(val_stat)==-1){ # if the sign of the test stat. is negative, it matched to the left boundary
          decision_2 <- paste('The p-value for the test statistic in the r=1 case is greater than ',p_value)
          # Accept (0) or reject (1) based on variable `significance_level`

          decision_3 <- 0 # if the statistic is so small that it matches to the left boundary, the p-value of that statistic is giant, we accept
        }
        #
        if (sign(val_stat)==1){ # if the sign of the test stat. is positive, it matched to the right boundary
          decision_2 <- paste('The p-value for the test statistic in the r=1 case is less than ',p_value)
          # Accept (0) or reject (1) based on variable `significance_level`
          # if significance level specified by the user is smaller than the smallest p-value that we know the statistic for based on the airy_1_tab then technically we don't know what the decision would be at the very small level of significance we specified so we return the decision based on the smallest p-value

          decision_3 <- 1 # if the statistic is so big that it matches to the right boundary, the p-value of that statistic is tiny, we reject
        }
      }else{
        decision_2 <- paste('The p-value for the test statistic in the r=1 case is',p_value)
        # Accept (0) or reject (1) based on variable `significance_level`
        if (p_value >= significance_level){
          decision_3 <- 0
        }else{
          decision_3 <- 1
        }

      }
      #intermediary list to append our final, output list with
      list_2 <- list("p_value"=p_value, "Precise decision"=decision_2, "boolean_decision" = decision_3)
      list_table <- append(list_table,list_2)
    }

    if (r==2|r==3){ #we still have a significance table we can go by when r=2 or r=3, hence this loop
      if (significance_level > 0.05){ #round the significance level according to the description document (significance_level var. description in ?largevar)
        round_s_f <- 0.1
      }else{
        if (significance_level > 0.025){
          round_s_f <- 0.05
        }else{
          if (significance_level > 0.01){
            round_s_f <- 0.025
          }else{
            round_s_f <- 0.01
          }
        }
      }

      if (r==2){

        if (round_s_f == 0.1){
          if (sum_airy_quantiles[2,1] >= LR_nt){
            decision_3 <- 0
          }else{
            decision_3 <- 1
          }
        }

        if (round_s_f == 0.05){
          if (sum_airy_quantiles[2,2] >= LR_nt){
            decision_3 <- 0
          }else{
            decision_3 <- 1
          }
        }

        if (round_s_f == 0.025){
          if (sum_airy_quantiles[2,3] >= LR_nt){
            decision_3 <- 0
          }else{
            decision_3 <- 1
          }
        }

        if (round_s_f == 0.01){
          if (sum_airy_quantiles[2,4] >= LR_nt){
            decision_3 <- 0
          }else{
            decision_3 <- 1
          }
        }


      }else{#r==3

        if (round_s_f == 0.1){
          if (sum_airy_quantiles[3,1] >= LR_nt){
            decision_3 <- 0
          }else{
            decision_3 <- 1
          }
        }

        if (round_s_f == 0.05){
          if (sum_airy_quantiles[3,2] >= LR_nt){
            decision_3 <- 0
          }else{
            decision_3 <- 1
          }
        }

        if (round_s_f == 0.025){
          if (sum_airy_quantiles[3,3] >= LR_nt){
            decision_3 <- 0
          }else{
            decision_3 <- 1
          }
        }

        if (round_s_f == 0.01){
          if (sum_airy_quantiles[3,4] >= LR_nt){
            decision_3 <- 0
          }else{
            decision_3 <- 1
          }
        }

      }

      list_2 <- list("significance_level"=round_s_f,"boolean_decision" = decision_3)
      list_table <- append(list_table,list_2)
    }

    list <- list("statistic"=LR_nt, "measure_upper_bound"=lambda_p, "measure_lower_bound"=lambda_m, "eigenvalues"=ev_values,"significance_test"=list_table,"k"=k,"r"=r)

    # Plot the function if user asked for it
    if (plot_output==TRUE){
      my_function <- function(x){(p+q)*sqrt(pmax(0,(lambda_p-x)*(x-lambda_m)))/x/(1-x)/2/pi}
      plot <- hist(ev_values, breaks = 2.5*(log2(length(ev_values))+1), probability = TRUE, col = "lightblue", border = "white", main = paste("VAR(",k,") Eigenvalues"),xlab = "Eigenvalues", ylab = "Frequency" ,xlim=c(0,1))
      curve(my_function,add = TRUE)
      #this guarantees that the histogram will pop up on its own in "Plots" after running the test
      list <- append(list,plot)
    }

    #the output of the function
    return(list)
  }
}

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

  plot <- hist(values,breaks=2.5*(log2(sim_num)+1))
  abline(v=x,col="red",lwd=3)

  list <- list("sim_results"=stat_vec,"empirical_percentage"=percentage,plot)
  return(list)
}


#' Cointegration test for settings of large N and T
#'
#' This is the "skeleton" version of the largevar function in the package. It is called within the sim_function function to make runtime faster. For the actual cointegration test, use the largevar function.
#' @param data a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @return The test statistic.
#' @export
#'

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
