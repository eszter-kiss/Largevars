# custom_format <- function(x) {
#   formatted <- character(length(x))  # Initialize an empty character vector to store formatted numbers
#   formatted <- character(length(table[,1:3]))  # Initialize an empty character vector to store formatted numbers
#   x <- table[,1:3]
#   for (i in seq_along(x)) {  # Loop through each element in the input vector 'x'
#     if (grepl("^-?\\d+\\.\\d0$", x[i])) {  # Check if the number matches the pattern xy.z0
#       formatted[i] <- x[i]  # Keep xy.z0 as is
#     } else if (grepl("^-?\\d+\\.\\d\\d$", x[i])) {  # Check if the number matches the pattern xy.zw
#       formatted[i] <- sub("^(-?\\d+\\.\\d)\\d$", "\\1", x[i])  # Remove the last digit 'w'
#     } else if (grepl("^-?\\d+\\.\\d$", x[i])) {  # Check if the number is of the form xy.z
#       formatted[i] <- paste0(x[i], "0")  # Append '0' to the end
#     } else if (grepl("^-?\\d+$", x[i])) {  # Check if the number is an integer (xy)
#       formatted[i] <- paste0(x[i], ".0")  # Append '.0' to the integer
#     } else {  # For other cases, keep the original format
#       formatted[i] <- x[i]
#     }
#   }
#
#   formatted  # Return the vector with formatted numbers
# }
#
#   formatted  # Return the vector with formatted numbers
#   table[,1:3]
# }

# Define the class for the output list of largevar function
setClass("stat_test", representation("list"))

# Define the custom show method
setMethod("show", "stat_test", function(object) {
  cat("Output for the largevars function", "\n")
  cat("===================================", "\n")
  cat("Cointegration test for high-dimensional VAR(k)\n")
  cat("\n")  # Print double dashed line
  cat("T = ", object$t,", N = ", object$N,"\n")
  cat("\n")  # Print double dashed line
  print(object$significance_test$significance_row)
  cat("\n")
  cat(object$significance_test$`Statistical decision`,"\n")
  cat("============================================================================", "\n")
  cat("Test statistic:", object$statistic ,"\n")
  cat(object$significance_test$text,"\n")
  cat("Decision about the H0: ", object$significance_test$boolean_decision,"\n")
})
