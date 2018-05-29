
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

  # initialise project with default values
  ret <- list(data = list(),
              parameter_sets = list(),
              active_set = 0,
              output = list())

  # create class and return invisibly
  class(ret) <- "malecot_project"
  invisible(ret)
}

#------------------------------------------------
# overload print() function to print without class specifier
#' @export
print.malecot_project <- function(x, ...) {

  # print raw list
  print(unclass(x))

  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function.
#' @export
summary.malecot_project <- function(object, ...) {

  proj <- object

  # print data summary
  cat("### DATA\n")
  if (length(proj$data)==0) {
    cat("   none loaded\n")
  } else {
    data <- proj$data$raw
    dataFormat <- proj$data$dataFormat
    missingData <- proj$data$missingData

    if (dataFormat=="biallelic") {
      nsamp <- nrow(data)
      nlocus <- ncol(data)
      nmissing <- sum(data == missingData)
    } else {
      nsamp <- length(unique(data[,1]))
      nlocus <- length(unique(data[,2]))
      nmissing <- sum(data[,3] == missingData)
    }

    cat(paste0("   format:\t", dataFormat, "\n"))
    cat(paste0("   samples:\t", nsamp,"\n"))
    cat(paste0("   loci:\t\t", nlocus, "\n"))
    cat(paste0("   missing:\t", nmissing, "/", nsamp*nlocus, " (", round(100*nmissing/(nsamp*nlocus),1), "%)\n"))
  }
  cat("\n")

  # print parameter sets summary
  cat("### PARAMETER SETS\n")
  s <- proj$activeSet
  if (length(proj$parameterSets)==0) {
    cat("   none defined\n")
  } else {
    cat(paste0("   sets:\t\t", length(proj$parameterSets), "\n"))
    cat(paste0("   active set:\t", s, "\n"))

    cat("\n   ## CURRENT SET\n")
    if (!is.null(proj$parameterSets[[s]]$setDescription)) {
      if (proj$parameterSets[[s]]$setDescription != "") {
        cat(paste0("   ", proj$parameterSets[[s]]$setDescription, "\n"))
      }
    }
    cat(paste0("   K_range =\t\t", paste(proj$parameterSets[[s]]$K_range, collapse=", "), "\n"))
    cat(paste0("   lambda =\t\t", proj$parameterSets[[s]]$lambda, "\n"))
    cat(paste0("   COI_model =\t", proj$parameterSets[[s]]$COI_model, "\n"))
    cat(paste0("   COI_max =\t\t", proj$parameterSets[[s]]$COI_max, "\n"))
    cat(paste0("   COI_dispersion =\t", proj$parameterSets[[s]]$COI_dispersion, "\n"))
    cat(paste0("   e1 =\t\t\t", proj$parameterSets[[s]]$e1, "\n"))
    cat(paste0("   e2 =\t\t\t", proj$parameterSets[[s]]$e2, "\n"))
    cat(paste0("   estimateError =\t", proj$parameterSets[[s]]$estimateError, "\n"))
    cat(paste0("   e1_max =\t\t", proj$parameterSets[[s]]$e1_max, "\n"))
    cat(paste0("   e2_max =\t\t", proj$parameterSets[[s]]$e2_max, "\n"))
  }
  cat("\n")

  # print output summary
  cat("### OUTPUT\n")
  if (length(proj$output)==0) {
    cat("   none\n")
  } else {
    if (length(proj$output[[s]])==0) {
      cat("   none\n")
    } else {
      cat("   THERE IS SOME OUTPUT\n")
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
