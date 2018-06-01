
#------------------------------------------------
# replace NULL value with default
# (not exported)
#' @noRd
define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}

# -----------------------------------
# user_yes_no
# ask user a yes/no question. Return TRUE/FALSE
# (not exported)
#' @noRd
user_yes_no <- function(x="continue? (Y/N): ") {
  userChoice <- NA
  while (!userChoice %in% c("Y", "y" ,"N", "n")) {
    userChoice <- readline(x)
  }
  return(userChoice %in% c("Y", "y"))
}

# -----------------------------------
# draw from Dirichlet distribution
# (not exported)
#' @noRd
rdirichlet <- function (alpha_vec) {
  Y <- rgamma(length(alpha_vec), shape = alpha_vec, scale = 1)
  output <- Y/sum(Y)
  return(output)
}

#------------------------------------------------
# forty_winks
# short sleep to help console output catch up
# (not exported)
#' @noRd
forty_winks <- function() {
  Sys.sleep(0.01)
}

# -----------------------------------
# mat_to_Rcpp
# takes matrix as input, converts to list format for use within Rcpp code
# (not exported)
#' @noRd
mat_to_rcpp <- function(x) {
  return(split(x, f=1:nrow(x)))
}

# -----------------------------------
# Rcpp_to_mat
# Takes list format returned from Rcpp and converts to matrix.
# (not exported)
#' @noRd
rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}
