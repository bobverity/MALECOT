
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
# calls C++ implementation of the Hungarian algorithm for binding best matching
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
  quantile(x, probs=c(0.025, 0.5, 0.975))
}

#------------------------------------------------
# return p-value of Geweke's diagnostic convergence statistic, estimated from package coda
#' @noRd
geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# check that geweke p-value non-significant on values x[1:n]
#' @noRd
test_convergence <- function(x, n) {
  if (n==1) {
    return(FALSE)
  }
  g <- geweke_pvalue(mcmc(x[1:n]))
  ret <- (g>0.01)
  if (is.na(ret)) {
    ret <- FALSE;
  }
  return(ret)
}

#------------------------------------------------
# update progress bar
#' @noRd
update_progress <- function(args, type, i, max_i) {
  
  # split by type
  if (type==1) { # scaffold progress bar
    setTxtProgressBar(args$pb_scaf, i)
    if (i==max_i) {
      close(args$pb_scaf)
    }
  } else if (type==2) { # burn-in iterations progress bar
    setTxtProgressBar(args$pb_burnin, i)
    if (i==max_i) {
      close(args$pb_burnin)
    }
  } else if (type==3) { # sampling iterations progress bar
    setTxtProgressBar(args$pb_samples, i)
    if (i==max_i) {
      close(args$pb_samples)
    }
  }
}

##########################################################################################################
# misc classes

#------------------------------------------------
#' @title Overload print function for malecot_loglike_quantiles
#'
#' @description Overload print function for malecot_loglike_quantiles
#'
#' @details TODO
#'
#' @param x TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

print.malecot_loglike_quantiles <- function(x, ...) {
  
  # print raw list
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}


#------------------------------------------------
#' @title Overload summary function for malecot_loglike_quantiles
#'
#' @description Overload summary function for malecot_loglike_quantiles
#'
#' @details TODO
#'
#' @param object TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

summary.malecot_loglike_quantiles <- function(object, ...) {
  
  # print raw summary
  summary(unclass(object))
  
  # return invisibly
  invisible(object)
}

#------------------------------------------------
#' @title Determine if object is of class malecot_loglike_quantiles
#'
#' @description Determine if object is of class malecot_loglike_quantiles
#'
#' @details TODO
#'
#' @param x TODO
#'
#' @export
#' @examples
#' # TODO

is.malecot_loglike_quantiles <- function(x) {
  inherits(x, "malecot_loglike_quantiles")
}

#------------------------------------------------
#' @title Overload print function for malecot_m_quantiles
#'
#' @description Overload print function for malecot_m_quantiles
#'
#' @details TODO
#'
#' @param x TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

print.malecot_m_quantiles <- function(x, ...) {
  
  # print raw list
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}

 
#------------------------------------------------
#' @title Overload summary function for malecot_m_quantiles
#'
#' @description Overload summary function for malecot_m_quantiles
#'
#' @details TODO
#'
#' @param object TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

summary.malecot_m_quantiles <- function(object, ...) {
  
  # print raw summary
  summary(unclass(object))
  
  # return invisibly
  invisible(object)
}

#------------------------------------------------
#' @title Determine if object is of class malecot_m_quantiles
#'
#' @description Determine if object is of class malecot_m_quantiles
#'
#' @details TODO
#'
#' @param x TODO
#'
#' @export
#' @examples
#' # TODO

is.malecot_m_quantiles <- function(x) {
  inherits(x, "malecot_m_quantiles")
}

#------------------------------------------------
#' @title Overload print function for malecot_p_quantiles
#'
#' @description Overload print function for malecot_p_quantiles
#'
#' @details TODO
#'
#' @param x TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

print.malecot_p_quantiles <- function(x, ...) {
  
  # print raw list
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}


#------------------------------------------------
#' @title Overload summary function for malecot_p_quantiles
#'
#' @description Overload summary function for malecot_p_quantiles
#'
#' @details TODO
#'
#' @param object TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

summary.malecot_p_quantiles <- function(object, ...) {
  
  # print raw summary
  summary(unclass(object))
  
  # return invisibly
  invisible(object)
}

#------------------------------------------------
#' @title Determine if object is of class malecot_p_quantiles
#'
#' @description Determine if object is of class malecot_p_quantiles
#'
#' @details TODO
#'
#' @param x TODO
#'
#' @export
#' @examples
#' # TODO

is.malecot_p_quantiles <- function(x) {
  inherits(x, "malecot_p_quantiles")
}

#------------------------------------------------
#' @title Overload print function for malecot_q_matrix
#'
#' @description Overload print function for malecot_q_matrix
#'
#' @details TODO
#'
#' @param x TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

print.malecot_q_matrix <- function(x, ...) {
  
  # print raw list
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}


#------------------------------------------------
#' @title Overload summary function for malecot_q_matrix
#'
#' @description Overload summary function for malecot_q_matrix
#'
#' @details TODO
#'
#' @param object TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

summary.malecot_q_matrix <- function(object, ...) {
  
  # print raw summary
  summary(unclass(object))
  
  # return invisibly
  invisible(object)
}

#------------------------------------------------
#' @title Determine if object is of class malecot_q_matrix
#'
#' @description Determine if object is of class malecot_q_matrix
#'
#' @details TODO
#'
#' @param x TODO
#'
#' @export
#' @examples
#' # TODO

is.malecot_q_matrix <- function(x) {
  inherits(x, "malecot_q_matrix")
}

#------------------------------------------------
#' @title Overload print function for malecot_GTI_path
#'
#' @description Overload print function for malecot_GTI_path
#'
#' @details TODO
#'
#' @param x TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

print.malecot_GTI_path <- function(x, ...) {
  
  # print raw list
  class(x) <- "data.frame"
  print(x)
  
  # return invisibly
  invisible(x)
}


#------------------------------------------------
#' @title Overload summary function for malecot_GTI_path
#'
#' @description Overload summary function for malecot_GTI_path
#'
#' @details TODO
#'
#' @param object TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

summary.malecot_GTI_path <- function(object, ...) {
  
  # print raw summary
  class(x) <- "data.frame"
  summary(x)
  
  # return invisibly
  invisible(object)
}

#------------------------------------------------
#' @title Determine if object is of class malecot_GTI_path
#'
#' @description Determine if object is of class malecot_GTI_path
#'
#' @details TODO
#'
#' @param x TODO
#'
#' @export
#' @examples
#' # TODO

is.malecot_GTI_path <- function(x) {
  inherits(x, "malecot_GTI_path")
}


