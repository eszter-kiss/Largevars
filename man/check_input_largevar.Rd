% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_input_largevar.R
\name{check_input_largevar}
\alias{check_input_largevar}
\title{Input checker for largevar function}
\usage{
check_input_largevar(
  data,
  k,
  r,
  fin_sample_corr,
  plot_output,
  significance_level
)
}
\arguments{
\item{data}{a numeric matrix where columns contain the individual time series that will be examined for presence of cointegrating relationships}

\item{k}{The number of lags we wish to employ in the VECM form (default: k=1)}

\item{r}{The number of cointegrating relationships we impose on the H1 hypothesis (default: r=1)}

\item{fin_sample_corr}{A boolean variable indicating whether we wish to employ finite sample correction on our test statistic. Default is false}

\item{plot_output}{A boolean variable indicating whether we wish to generate the distribution of the eigenvalues (default: TRUE)}

\item{significance_level}{Specify the significance level at which the decision about the H0 should be made. For r=1 this can be any level of significance. For r=2 and r=3, the significance level input will be rounded up to the nearest of the following: 0.1, 0.05, 0.025, 0.01. If the significance level is larger than 0.1, the decision will be made at the 10\% level. For r>3 only the test statistic is returned. For an empirical p-value for r>3 use the sim_function fun. in the package.}
}
\value{
Nothing (or warning message) if all inputs are correct, and an error message otherwise.
}
\description{
This is an internal function that checks the validity of the inputs of the largevar function.
}
\keyword{internal}
