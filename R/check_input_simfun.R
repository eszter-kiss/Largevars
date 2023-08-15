#' Cointegration test for settings of large N and T
#'
#' This is an internal function that checks the validity of the inputs of the sim_function function.
#' @param data a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @param sim_num The number of simulation we wish to run.
#' @keywords internal

check_input_simfun <- function(data,k,r,sim_num){

  if(is.null(data)){ # data must be existing input (default is NULL if nothing is input)
    print("`data` is a mandatory input")
    return(1)
  }else if((is.matrix(data)==FALSE)|(is.numeric(data)==FALSE)){
    print("`data` must be a numeric matrix.")
    return(1)

  }else if((is.numeric(k)==FALSE)|(length(k) == 1)==FALSE){
    print("`k` must be a number.")
    return(1)
  }else if(((k%%1==0)==FALSE)|((k>0)==FALSE)){
    print("`k` must be a positive integer.")
    return(1)
  }else if  ( k >= ((dim(data)[1]-1)/(dim(data)[2]-1)) ){
    print("`k` must be such that k+1 < T/N holds.")
    return(1)

  }else if((is.numeric(r)==FALSE)|(length(r) == 1)==FALSE){
    print("`r` must be a number.")
    return(1)
  }else if(((r%%1==0)==FALSE)|((r>0)==FALSE)){
    print("`r` must be a positive integer.")
    return(1)
  }else if  ( r>dim(data)[2] ){
    print("`r` must be less than or equal to the number of variables in your dataset.")
    return(1)

  }else if((is.numeric(sim_num)==FALSE)|(length(sim_num) == 1)==FALSE){
    print("`sim_num` must be a number.")
    return(1) # arg. 9

  }else if(((sim_num%%1==0)==FALSE)|((sim_num>0)==FALSE)){
    print("`sim_num` must be a positive integer.")
    return(1)

  }else if  ( (isTRUE(fin_sample_corr)==FALSE) & (isFALSE(fin_sample_corr)==FALSE) ){
    print("`fin_sample_corr` must be a boolean.")
    return(1)

  } else if (sim_num>500){
    print("Warning: simulation may run for several minutes")
    return(0)
  }  else{
    return(0)
  }
}
