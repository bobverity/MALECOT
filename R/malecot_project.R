
#------------------------------------------------
#' @title Define empty malecot_project object
#'
#' @description Define empty malecot_project object
#'
#' @details TODO
#'
#' @export
#' @examples
#' # TODO
#'
malecot_project <- function() {
  
  # create some empty data frames for storing results
  GTI_logevidence_model <- data.frame(set = numeric(),
                                      name = character(),
                                      mean = numeric(),
                                      SE = numeric())
  class(GTI_logevidence_model) <- "maverick_GTI_logevidence_model"
  
  GTI_posterior_model <- data.frame(set = numeric(),
                                    name = character(),
                                    Q2.5 = numeric(),
                                    Q50 = numeric(),
                                    Q97.5 = numeric())
  class(GTI_posterior_model) <- "maverick_GTI_posterior_model"
  
  # initialise project with default values
  ret <- list(data = NULL,
              data_processed = NULL,
              parameter_sets = NULL,
              active_set = 0,
              output = list(single_set = list(),
                            all_sets = list(GTI_logevidence_model = GTI_logevidence_model,
                                            GTI_posterior_model = GTI_posterior_model)
              )
  )
  
  # create class and return invisibly
  class(ret) <- "malecot_project"
  invisible(ret)
}

#------------------------------------------------
#' @title Custom print function for class malecot_project
#'   
#' @description Custom print function for class malecot_project, printing a 
#'   summary of the key elements (also equivalent to \code{summary(x)}). To do 
#'   an ordinary \code{print()} of all elements of the project, use the 
#'   \code{print_full()} function.
#'   
#' @param x object of class \code{malecot_project}
#' @param ... other arguments (ignored)
#'   
#' @export

print.malecot_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Ordinary print function for class malecot_project
#'
#' @description Calling \code{print()} on an object of class malecot_project
#'   results in custom output. This function therefore stands in for the base
#'   \code{print()} function, and is equivalent to running
#'   \code{print(unclass(x))}.
#'
#' @param x object of class \code{malecot_project}
#' @param ... other arguments passed to \code{print()}
#'
#' @export

print_full <- function(x, ...) {
  
  # check inputs
  assert_custom_class(x, "malecot_project")
  
  # print un-classed object
  print(unclass(x), ...)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @title Print summary for class malecot_project
#'   
#' @description Overload summary function for class malecot_project
#'   
#' @param object object of class \code{malecot_project}
#' @param ... other arguments (ignored)
#'   
#' @export

summary.malecot_project <- function(object, ...) {
  
  # print data summary
  cat("DATA:\n")
  if (is.null(object$data)) {
    cat("   (none loaded)\n")
  } else {
    
    # extract data properties
    n <- object$data_processed$n
    L <- object$data_processed$L
    pop <- object$data_processed$pop
    data_format <- object$data_processed$data_format
    name <- object$data_processed$name
    
    if (!is.null(name)) {
      cat(sprintf("   '%s'\n", name))
    }
    cat(sprintf("   data format = %s\n", data_format))
    cat(sprintf("   samples = %s\n", n))
    cat(sprintf("   loci = %s\n", L))
    cat(sprintf("   pops = %s\n", length(unique(pop))))
    if (data_format == "multiallelic") {
      allele_range <- range(object$data_processed$alleles)
      if (allele_range[1] == allele_range[2]) {
        cat(sprintf("   alleles = %s at all loci\n", allele_range[1]))
      } else {
        cat(sprintf("   alleles = between %s and %s\n", allele_range[1], allele_range[2]))
      }
    }
    
    switch(data_format,
           "biallelic" = {
             n_missing <- sum(object$data_processed$data == 0)
           },
           "multiallelic" = {
             n_missing <- sum(object$data_processed$data$haplotype == 0)
           })
    cat(sprintf("   missing data = %s of %s gene copies (%s%%)\n", n_missing, n*L, round(n_missing/(n*L)*100)))
  }
  cat("\n")
  
  
  # print parameter sets summary
  cat("PARAMETER SETS:\n")
  if (length(object$parameter_sets)==0) {
    cat("   (none defined)\n")
  } else {
    # print names of all sets
    s <- object$active_set
    for (i in 1:length(object$parameter_sets)) {
      
      # star next to active set
      if (i==s) {
        cat(" * ")
      } else {
        cat("   ")
      }
      
      # print name of set
      cat(sprintf("SET%s: %s\n", i, object$parameter_sets[[i]]$name))
    }
    cat("\n")
    
    # print details of active set
    name <- object$parameter_sets[[s]]$name
    lambda <- object$parameter_sets[[s]]$lambda
    lambda_range <- range(unlist(lambda))
    COI_model <- object$parameter_sets[[s]]$COI_model
    COI_max <- object$parameter_sets[[s]]$COI_max
    COI_manual <- object$parameter_sets[[s]]$COI_manual
    estimate_COI_mean <- object$parameter_sets[[s]]$estimate_COI_mean
    COI_mean <- object$parameter_sets[[s]]$COI_mean
    COI_dispersion <- object$parameter_sets[[s]]$COI_dispersion
    estimate_error <- object$parameter_sets[[s]]$estimate_error
    e1 <- object$parameter_sets[[s]]$e1
    e2 <- object$parameter_sets[[s]]$e2
    e1_max <- object$parameter_sets[[s]]$e1_max
    e2_max <- object$parameter_sets[[s]]$e2_max
    
    cat(sprintf("ACTIVE SET: SET%s\n", s))
    if (lambda_range[1] == lambda_range[2]) {
      cat(sprintf("   lambda = %s\n", lambda_range[1]))
    } else {
      cat(sprintf("   lambda = between %s and %s\n", lambda_range[1], lambda_range[2]))
    }
    cat(sprintf("   COI model = %s\n", COI_model))
    cat(sprintf("   COI max = %s\n", COI_max))
    if (any(COI_manual != -1)) {
      cat(sprintf("   COI specified manually for %s samples\n", sum(COI_manual != -1)))
    }
    cat(sprintf("   estimate COI mean = %s\n", estimate_COI_mean))
    if (!estimate_COI_mean) {
      cat(sprintf("   COI mean = %s\n", COI_mean))
    }
    cat(sprintf("   COI dispersion = %s\n", COI_dispersion))
    cat(sprintf("   estimate error = %s\n", estimate_error))
    if (estimate_error) {
      cat(sprintf("   e1_max = %s\n", e1_max))
      cat(sprintf("   e2_max = %s\n", e2_max))
    } else {
      cat(sprintf("   e1 = %s\n", e1))
      cat(sprintf("   e2 = %s\n", e2))
    }
    
  }
  cat("\n")
  
}

#------------------------------------------------
#' @title Determine if object is of class malecot_project
#'
#' @description Determine if object is of class malecot_project.
#'
#' @details TODO
#'
#' @param x TODO
#'
#' @export
#' @examples
#' # TODO
#'
is.malecot_project <- function(x) {
  inherits(x, "malecot_project")
}
