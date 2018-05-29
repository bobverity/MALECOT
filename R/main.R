#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib MALECOT
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
#' Dummy function
#'
#' Simple test function that demonstrates some of the features of the
#' bellsandwhistles package.
#'
#' @param x vector of values
#'
#' @export
#' @examples
#' # Find square of first 100 values
#' dummy1(1:100)

dummy1 <- function(x = 1:5) {

  # print message to console
  message("running R dummy1 function")

  # get arguments in list form
  args <- list(x = x)

  # run C++ function with these arguments
  output_raw <- dummy1_cpp(args)

  # some optional processing of output
  message("processing output")
  ret <- list(x = output_raw$x_squared)

  # return
  return(ret)
}

