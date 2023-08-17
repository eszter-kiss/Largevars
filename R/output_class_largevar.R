# Define the class for the output list of largevar function
setClass("stat_test", representation("list"))

# Define the custom show method
setMethod("show", "stat_test", function(object) {
  cat("Output for the largevars function", "\n")
  cat("===================================", "\n")
  cat("Row of statistics\n")
  cat("\n")  # Print double dashed line
  print(object$significance_test$significance_row)
  cat("\n")
  cat(object$significance_test$`Statistical decision`,"\n")
  cat("============================================================================", "\n")
  cat("Test statistic:", object$statistic ,"\n")
  cat(object$significance_test$text,"\n")
  cat("Decision about the H0: ", object$significance_test$boolean_decision,"\n")
})
