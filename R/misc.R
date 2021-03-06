
#------------------------------------------------
#' @title Import file
#'
#' @description Import file from the inst/extdata folder of this package
#' 
#' @param name name of file
#'
#' @export

malecot_file <- function(name) {
  
  # load file from inst/extdata folder
  name_full <- system.file("extdata", name, package = 'MALECOT', mustWork = TRUE)
  ret <- readRDS(name_full)
  
  # return
  return(ret)
}

#------------------------------------------------
# replace NULL value with default
#' @noRd
define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}

#------------------------------------------------
# simple zero-padding function. Not robust to e.g. negative numbers
#' @noRd
zero_pad_simple <- function(x, n = 3) {
  ret <- mapply(function(x) {
                  paste0(paste0(rep(0,n-nchar(x)), collapse = ""), x, collapse = "")
                }, x)
  return(ret)
}

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE
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
#' @noRd
rdirichlet <- function (alpha_vec) {
  Y <- rgamma(length(alpha_vec), shape = alpha_vec, scale = 1)
  output <- Y/sum(Y)
  return(output)
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
mat_to_rcpp <- function(x) {
  return(split(x, f=1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}

#------------------------------------------------
# calls C++ implementation of the Hungarian algorithm for finding best matching
# in a linear sum assigment problem. This is function is used in testing.
#' @noRd
call_hungarian <- function(x) {
  args <- list(cost_mat = mat_to_rcpp(x))
  call_hungarian_cpp(args)
}

#------------------------------------------------
# return 95% quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs=c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
# sum logged values without underflow, i.e. do log(sum(exp(x)))
#' @noRd
log_sum <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  x_max <- max(x, na.rm = TRUE)
  ret <- x_max + log(sum(exp(x-x_max)))
  return(ret)
}

#------------------------------------------------
# return p-value of Geweke's diagnostic convergence statistic, estimated from package coda
#' @noRd
geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# check convergence on values x[1:n]
#' @noRd
test_convergence <- function(x, n) {
  # fail if n = 1
  if (n == 1) {
    return(FALSE)
  }
  
  # fail if ESS too small
  ESS <- coda::effectiveSize(x[1:n])
  if (ESS < 10) {
    return(FALSE)
  }
  
  # fail if geweke p-value < threshold
  g <- geweke_pvalue(mcmc(x[1:n]))
  ret <- (g > 0.01)
  if (is.na(ret)) {
    ret <- FALSE;
  }
  return(ret)
}

#------------------------------------------------
# update progress bar
# (not exported)
#' @noRd
update_progress <- function(pb_list, name, i, max_i) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i==max_i) {
    close(pb_list[[name]])
  }
}


##########################################################################################################
# MISC CLASSES

#------------------------------------------------
# Overload print function for malecot_qmatrix
#' @noRd
print.malecot_qmatrix <- function(x, ...) {
  print(unclass(x))
  invisible(x)
}

#------------------------------------------------
# Overload print function for malecot_loglike_intervals
#' @noRd
print.malecot_loglike_intervals <- function(x, ...) {
  print(unclass(x))
  invisible(x)
}

#------------------------------------------------
# Overload print function for malecot_COI_intervals
#' @noRd
print.malecot_COI_intervals <- function(x, ...) {
  print(unclass(x))
  invisible(x)
}

#------------------------------------------------
# Overload print function for malecot_p_intervals
#' @noRd
print.malecot_p_intervals <- function(x, ...) {
  print(unclass(x))
  invisible(x)
}

#------------------------------------------------
# Overload print function for malecot_e_intervals
#' @noRd
print.malecot_e_intervals <- function(x, ...) {
  print(unclass(x))
  invisible(x)
}

#------------------------------------------------
# Overload print function for malecot_COI_mean_intervals
#' @noRd
print.malecot_COI_mean_intervals <- function(x, ...) {
  print(unclass(x))
  invisible(x)
}

#------------------------------------------------
# Overload print function for malecot_GTI_path
#' @noRd
print.malecot_GTI_path <- function(x, ...) {
  class(x) <- "data.frame"
  print(x)
  invisible(x)
}

#------------------------------------------------
# Overload print function for class malecot_GTI_logevidence
#' @noRd
print.malecot_GTI_logevidence <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Overload print function for class malecot_GTI_posterior
#' @noRd
print.malecot_GTI_posterior <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Overload print function for class malecot_GTI_logevidence_model
#' @noRd
print.malecot_GTI_logevidence_model <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}

#------------------------------------------------
# Overload print function for class malecot_GTI_posterior_model
#' @noRd
print.malecot_GTI_posterior_model <- function(x, ...) {
  print(as.data.frame(unclass(x)))
  invisible(x)
}
