## This script serves to run largevar() as it was at the time of development.
## The main output values are saved into an .rds file. These values are used
## in subsequent testing files for largevar(). To make sure that during package
## building the original .rds file with values is not overwritten (and so the
## testing function for largevar stays meaningful), this script is fully
## commented out and distributed this way.

# library(Largevars)

# # Load the sample data
# data("s_p100_price")
#
# # Transform data according to researcher needs
# dataSP <- log(s_p100_price[, seq(2, dim(s_p100_price)[2])])
# dataSP <- as.matrix(dataSP)
#
# # Run the function
# result <- largevar(data = dataSP, k = 1, r = 1, fin_sample_corr = FALSE,
#                    plot_output = TRUE, significance_level = 0.05)
#
# # Extract relevant components
# expected_result <- list(
#   statistic = result$statistic,
#   p_value = result$significance_test$p_value,
#   boolean_decision = result$significance_test$boolean_decision,
#   eigenvalues = result$eigenvalues
# )
#
# # Save the expected result
# saveRDS(expected_result, file = "tests/testthat/testdata/expected_largevar_result.rds")
