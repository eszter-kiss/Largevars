#' Cointegration test for settings of large N and T
#'
#' This is an internal function that checks the validity of the inputs of the sim_function function.
#' @param N a number representing the number of time series we want to simulate in the system
#' @param tau a number representing the length of the time series we want to simulate in the system
#' @param stat_value the test statistic value we want to calculate p-value based on
#' @param k The number of lags we wish to employ in the VECM form (default: k=1)
#' @param r The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)
#' @param sim_num The number of simulation we wish to run.
#' @keywords internal

check_input_simfun <- function(N,tau,stat_value,k,r,fin_sample_corr,sim_num){
  if(is.null(N)){ # (default is NULL if nothing is input)
    print("`N` is a mandatory input")
    return(1)
  }else if((is.numeric(N)==FALSE)|(length(N) == 1)==FALSE){
    print("`N` must be a number.")
    return(1)
  }else if(((N%%1==0)==FALSE)|((N>0)==FALSE)){
    print("`N` must be a positive integer.")
    return(1)

    }else if(is.null(tau)){ # (default is NULL if nothing is input)
      print("`tau` is a mandatory input")
      return(1)
    }else if((is.numeric(tau)==FALSE)|(length(tau) == 1)==FALSE){
      print("`tau` must be a number.")
      return(1)
    }else if(((tau%%1==0)==FALSE)|((tau>0)==FALSE)){
      print("`tau` must be a positive integer.")
      return(1)

    }else if(is.null(stat_value)){ # (default is NULL if nothing is input)
      print("`tau` is a mandatory input")
      return(1)
    }else if((is.numeric(stat_value)==FALSE)|(length(stat_value) == 1)==FALSE){
      print("`stat_value` must be a number.")
      return(1)

  }else if((is.numeric(k)==FALSE)|(length(k) == 1)==FALSE){
    print("`k` must be a number.")
    return(1)
  }else if(((k%%1==0)==FALSE)|((k>0)==FALSE)){
    print("`k` must be a positive integer.")
    return(1)
  }else if  ( k >= ((tau-1)/N)-1){
    print("`k` must be such that k+1 < T/N holds.")
    return(1)

  }else if((is.numeric(r)==FALSE)|(length(r) == 1)==FALSE){
    print("`r` must be a number.")
    return(1)
  }else if(((r%%1==0)==FALSE)|((r>0)==FALSE)){
    print("`r` must be a positive integer.")
    return(1)
  }else if  ( r>N){
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
