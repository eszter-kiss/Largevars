#' Cointegration test for settings of large N and T
#'
#' This is an internal function that checks the validity of the inputs of the largevar function.
#' @param data a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @param fin_sample_corr A boolean variable indicating whether we wish to employ finite sample correction on our test statistic. Default is false
#' @param plot_output A boolean variable indicating whether we wish to generate the distribution of the eigenvalues (default: TRUE)
#' @param significance_level Specify the significance level at which the decision about the H0 should be made. For r=1 this can be any level of significance. For r=2 and r=3, the significance level input will be rounded up to the nearest of the following: 0.1, 0.05, 0.025, 0.01. If the significance level is larger than 0.1, the decision will be made at the 10% level. For r>3 only the test statistic is returned. For an empirical p-value for r>3 use the sim_function fun. in the package.
#' @keywords internal

check_input_largevar <- function(data,k,r, fin_sample_corr, plot_output, significance_level){

  if(is.null(data)){ # data must be existing input (default is NULL if nothing is input)
    print("`data` is a mandatory input")
    return(1) # arg. 1
  }else if((is.matrix(data)==FALSE)|(is.numeric(data)==FALSE)){
    print("`data` must be a numeric matrix.")
    return(1) # arg. 2

  }else if((is.numeric(k)==FALSE)|(length(k) == 1)==FALSE){
    print("`k` must be a number.")
    return(1) # arg. 3
  }else if(((k%%1==0)==FALSE)|((k>0)==FALSE)){
    print("`k` must be a positive integer.")
    return(1) # arg. 4
  }else if  (k >= ((dim(data)[1]-1)/dim(data)[2])-1){
    print("`k` must be such that k+1 < T/N holds.")
    return(1) # arg. 5

  }else if((is.numeric(r)==FALSE)|(length(r) == 1)==FALSE){
    print("`r` must be a number.")
    return(1) # arg. 6
  }else if(((r%%1==0)==FALSE)|((r>0)==FALSE)){
    print("`r` must be a positive integer.")
    return(1) # arg. 7
  }else if  ( r>dim(data)[2] ){
    print("`r` must be less than or equal to the number of variables in your dataset.")
    return(1) # arg. 8

  }else if((is.numeric(significance_level)==FALSE)|(length(significance_level) == 1)==FALSE){
    print("`significance_level` must be a number.")
    return(1) # arg. 9
  }else if( (significance_level<0)|(significance_level>1) ){
    print("`significance_level` must be a real number strictly between 0 and 1.")
    return(1) # arg. 10

  }else if  ( (isTRUE(plot_output)==FALSE) & (isFALSE(plot_output)==FALSE) ){
    print("`plot_output` must be a boolean.")
    return(1) # arg. 11

  }else if  ( (isTRUE(fin_sample_corr)==FALSE) & (isFALSE(fin_sample_corr)==FALSE) ){
    print("`fin_sample_corr` must be a boolean.")
    return(1) # arg. 12

  } else if (r>10){
    print("Warning: test statistic percentiles are only available for `r` <= 10.")
    return(0) # arg. 13
  }  else{ #if none of the above fails
    return(0)
  }
}
