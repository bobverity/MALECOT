
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib MALECOT
#' @import assertthat
#' @import parallel
#' @import coda
#' @importFrom Rcpp sourceCpp
#' @import stats
#' @import utils
NULL

#------------------------------------------------
#' @title Bind data to project
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param proj TODO
#' @param data TODO
#' @param data_format TODO
#' @param name TODO
#' @param missing_data TODO
#' @param check_delete_output TODO
#'
#' @export
#' @examples
#' # TODO

bind_data <- function(proj, data, data_format = NULL, name = NULL, missing_data = -1, check_delete_output = TRUE) {

  # check inputs
  assert_malecot_project(proj)
  if (!is.null(data_format)) {
    assert_in(data_format, c("biallelic", "multiallelic"))
  }

  # try to work out data_format automatically if unspecified
  if (is.null(data_format)) {
    all_biallelic_format <- all(apply(data, 1, function(x){
      x %in% c(0, 0.5, 1, missing_data)
      }))
    data_format <- ifelse(all_biallelic_format, "biallelic", "multiallelic")
    message(sprintf("%s format detected", data_format))
  }

  # reformat missing data as -1
  data[data==missing_data] <- -1

  # perform checks on data
  # TODO - more checks on data format to ensure correct
  if (data_format=="multiallelic") {
    assert_that(ncol(data)==3)
    if (any( data[,3]<=0 & data[,3]!=missing_data )) {
      stop("for the multi-allelic format haplotypes must be coded as positive integers")
    }
    n <- length(unique(data[,1]))
    L <- length(unique(data[,2]))
  } else {
    n <- nrow(data)
    L <- ncol(data)
  }

  # check whether there is data loaded already
  if (length(proj$data)!=0) {
    # return existing project if user not happy to continue
    if (check_delete_output) {
      user_continue <- user_yes_no("All existing output for this project will be lost. Continue? (Y/N): ")
      if (!user_continue) {
        message("returning without modification\n")
        return(proj)
      }
    }

    # drop all existing output
    message("overwriting data and deleting old output\n")
    proj$output <- replicate(length(proj$parameter_sets), NULL)
  }

  # update project with new data
  proj$data <- list(raw = as.matrix(data),
                    name = name,
                    n = n,
                    L = L,
                    missing_data = missing_data,
                    data_format = data_format)

  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Create new MALECOT parameter set
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param proj the current MALECOT project
#' @param set_description brief description of this parameter set
#' @param K_range TODO
#' @param lambda TODO
#' @param COI_model TODO
#' @param COI_max TODO
#' @param COI_dispersion TODO
#' @param estimate_error TODO
#' @param e1 TODO
#' @param e2 TODO
#' @param e1_max TODO
#' @param e2_max TODO
#'
#' @export
#' @examples
#' # TODO

new_set <- function(proj, set_description = NULL, K_range = 1:3, lambda = 1.0, COI_model = "poisson", COI_max = 20, COI_dispersion = 1, estimate_error = FALSE, e1 = 0, e2 = 0, e1_max = 0.2, e2_max = 0.2) {

  # check inputs
  assert_malecot_project(proj)
  assert_pos_int(K_range, zero_allowed = FALSE)
  assert_bounded(lambda, left = 0, right = 100)
  assert_in(COI_model, c("uniform", "poisson", "nb"))
  assert_pos_int(COI_max, zero_allowed = FALSE)
  if (COI_model=="nb") {
    assert_gr(COI_dispersion, 1)
  }
  assert_bounded(e1_max)
  assert_bounded(e2_max)
  if (estimate_error) {
    assert_bounded(e1, right = e1_max)
    assert_bounded(e2, right = e2_max)
  }

  # count current parameter sets and add one
  s <- length(proj$parameter_sets) + 1

  # make new set active
  proj$active_set <- s

  # create new parameter set
  proj$parameter_sets[[s]] <- list(set_description = set_description,
                                   K_range = K_range,
                                   lambda = lambda,
                                   COI_model = COI_model,
                                   COI_max = COI_max,
                                   COI_dispersion = COI_dispersion,
                                   e1 = e1,
                                   e2 = e2,
                                   estimate_error = estimate_error,
                                   e1_max = e1_max,
                                   e2_max = e2_max)

  # add matching output object
  proj$output[s] <- list(NULL)

  # rename sets
  set_names <- paste0("set", 1:s)
  names(proj$parameter_sets) <- set_names
  names(proj$output) <- set_names

  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Delete MALECOT parameter set
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param proj the current MALECOT project
#' @param set numerical index of which set to delete
#' @param check_delete_output whether to perform check before deleting output corresponding to this parameter set
#'
#' @export
#' @examples
#' # TODO

delete_set <- function(proj, set = NULL, check_delete_output = TRUE) {

  # target set is active set by default
  s <- define_default(set, proj$active_set)

  # check inputs
  assert_malecot_project(proj)
  assert_pos_int(s, zero_allowed = TRUE)
  assert_bounded(s, left = 1, right = length(proj$parameter_sets), inclusive_left = TRUE, inclusive_right = TRUE)

  # look for output associated with this set
  if (!is.null(proj$output[[s]])) {

    # return existing project if user not happy to continue
    if (check_delete_output) {
      user_continue <- user_yes_no("All existing output associated with this parameter set will be lost. Continue? (Y/N): ")
      if (!user_continue) {
        message("returning without modification\n")
        return(proj)
      }
    }
    message("deleting old output\n")
  }

  # drop target from parameter sets and output
  proj$parameter_sets[[s]] <- NULL
  proj$output[[s]] <- NULL
  n_s <- length(proj$parameter_sets)

  # move active_set down 1
  proj$active_set <- proj$active_set - 1

  # rename parameter sets and output
  if (n_s>0) {
    set_names <- paste0("set", 1:n_s)
    names(proj$parameter_sets) <- set_names
    names(proj$output) <- set_names
  }

  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Example MCMC
#'
#' @description Demonstrates good algorithmic and coding practices for a simple normal mixture model. Implements advanced MCMC methods, including Metropolis-coupling over temperature rungs, a split-merge proposal, and an additional form of proposal referred to here as a "scaffold" proposal. Uses Rcpp and the package parallel for speed and efficiency.
#'
#' @details TODO
#'
#' @param proj TODO
#' @param K the number of mixture components
#' @param burnin the number of burn-in iterations
#' @param samples the number of sampling iterations
#' @param rungs the number of temperature rungs
#' @param auto_converge whether burn-in should be diagnosed automatically
#' @param coupling_on whether to implement Metropolis-coupling over temperature rungs
#' @param scaffold_on whether to use scaffolds to improve mixing
#' @param scaffold_n the number of scaffolds to use
#' @param split_merge_on whether to implement a split-merge proposal
#' @param solve_label_switching_on whether to implement Stevens' solution to the label-switching problem
#' @param precision the level of precision that allele frequencies are represented at
#' @param GTI_pow the power used in the generalised thermodynamic integration method for estimating K
#' @param silent supresses output messages
#' @param output_format choose from a range of output formats
#' @param cluster pass in a cluster environment
#'
#' @export
#' @examples
#' # TODO

run_mcmc <- function(proj, K = NULL, burnin = 1e2, samples = 1e3, rungs = 1, auto_converge = TRUE, coupling_on = TRUE, scaffold_on = TRUE, scaffold_n = 10, split_merge_on = TRUE, solve_label_switching_on = TRUE, precision = 0.01, GTI_pow = 2, silent = FALSE, output_format = 1, cluster = NULL) {
  
  # ---------- check inputs ----------
  
  # check input arguments
  assert_malecot_project(proj)
  assert_pos_int(burnin)
  assert_greq(burnin, 10)
  assert_pos_int(samples, zero_allowed = FALSE)
  assert_pos_int(rungs, zero_allowed = FALSE)
  assert_pos_int(scaffold_n, zero_allowed = FALSE)
  assert_int(GTI_pow)
  assert_bounded(GTI_pow, left = 1, right = 10)
  assert_in(output_format, 1:2)
  
  # get active set and check non-zero
  s <- proj$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # check that precision is either zero or leads to simple sequence
  if (precision!=0) {
    if(!((1/precision)%%1)==0) {
      stop("1/precision must be an integer")
    }
  }
  
  # default value and checks on K
  K <- define_default(K, proj$parameter_sets[[s]]$K_range)
  assert_in(K, proj$parameter_sets[[s]]$K_range)
  
  
  # ---------- process data ----------
  
  # if bi-allelic
  data <- proj$data$raw
  L <- proj$data$L
  data_format <- proj$data$data_format
  missing_data <- proj$data$missing_data
  if (data_format == "biallelic") {
  
    # convert data to simple integer format for use in Rcpp function
    data[data == 0.5] <- 2
    data[data == 0] <- 3
    data[data == missing_data] <- 0
    
    # number of haplotypes at each locus
    Jl <- rep(2,L)
  }
  
  # if multi-allelic
  if (data_format == "multiallelic") {
    
    # convert data to simple integer format for use in Rcpp function
    # TODO - check can convert data from text to numeric
    Jl <- rep(0,L)	# number of haplotypes at each locus
    for (l in 1:L) {
      u <- sort(unique(data[,3][data[,2]==l]))
      u <- u[u != missing_data]
      data[,3][data[,2]==l] <- match(data[,3][data[,2]==l], u)
      Jl[l] <- length(u)
    }
    data[,3][is.na(data[,3])] <- 0
  }
  
  
  # ---------- create argument lists ----------
  
  # define data list
  args_data <- list(data = mat_to_rcpp(data),
                    n = proj$data$n,
                    L = proj$data$L,
                    Jl = Jl)
  
  # get parameters from active set
  args_model <- proj$parameter_sets[[s]]
  
  # add to model parameters list
  args_model <- c(args_model, list(burnin = burnin,
                                   samples = samples,
                                   rungs = rungs,
                                   auto_converge = auto_converge,
                                   coupling_on = coupling_on,
                                   scaffold_on = scaffold_on,
                                   scaffold_n = scaffold_n,
                                   split_merge_on = split_merge_on,
                                   solve_label_switching_on = solve_label_switching_on,
                                   precision = precision,
                                   GTI_pow = GTI_pow,
                                   silent = silent,
                                   output_format = output_format,
                                   cluster = cluster))
  
  # get COI_model in numeric form
  args_model$COI_model_numeric <- match(args_model$COI_model, c("uniform", "poisson" ,"nb"))
  
  # R functions to pass to Rcpp
  args_functions <- list(forty_winks = forty_winks,
                         test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # define final argument list over all K
  parallel_args <- list()
  for (i in 1:length(K)) {
    
    # create progress bars
    pb_scaf <- txtProgressBar(min = 0, max = scaffold_n, initial = NA, style = 3)
    pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
    pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
    args_pb <- list(pb_scaf = pb_scaf,
                    pb_burnin = pb_burnin,
                    pb_samples = pb_samples)
    
    # create argument list
    parallel_args[[i]] <- c(list(K = K[i]), args_data, args_model, args_functions, args_pb)
  }
  
  
  # ---------- run MCMC ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) {  # run in parallel
    if (!inherits(cluster, "cluster")) {
      stop("expected a cluster object")
    }
    clusterEvalQ(cluster, library(MALECOT))
    if (data_format == "biallelic") {
      output_raw <- clusterApplyLB(cl = cluster, parallel_args, run_mcmc_biallelic_cpp)
    }
    if (data_format == "multiallelic") {
      output_raw <- clusterApplyLB(cl = cluster, parallel_args, run_mcmc_multiallelic_cpp)
    }
  } else {  # run in serial
    if (data_format == "biallelic") {
      output_raw <- lapply(parallel_args, run_mcmc_biallelic_cpp)
    }
    if (data_format == "multiallelic") {
      output_raw <- lapply(parallel_args, run_mcmc_multiallelic_cpp)
    }
  }
  
  return(output_raw)
  #------------------------

  # begin processing results
  message("Processing results\n")

  # loop through K
  ret <- list()
  for (i in 1:length(K)) {

    # get loglikelihood
    loglike <- t(rcpp_to_mat(output_raw[[i]]$loglike))

    # split into burn-in vs sampling iterations
    loglike_burnin <- mcmc(loglike[1:burnin,])
    loglike_sampling <- mcmc(loglike[-(1:burnin),])

    # get quantiles over sampling loglikelihoods
    loglike_quantiles <- t(apply(loglike_sampling, 2, function(x){quantile(x, probs=c(0.025, 0.5, 0.975))}))

    # process mu draws
    mu <- rcpp_to_mat(output_raw[[i]]$mu)

    # process Q-matrix
    qmatrix <- rcpp_to_mat(output_raw[[i]]$qmatrix)

    # process acceptance rates
    mc_accept <- 100 * output_raw[[i]]$mc_accept/(burnin+samples)
    scaf_accept <- 100 * output_raw[[i]]$scaf_accept/(burnin+samples)
    splitmerge_accept <- 100 * output_raw[[i]]$splitmerge_accept/(burnin+samples)

    # store results of this K
    ret[[i]] <- list(loglike_burnin = loglike_burnin, loglike_sampling = loglike_sampling, loglike_quantiles = loglike_quantiles, mu = mu, qmatrix = qmatrix, mc_accept = mc_accept, scaf_accept = scaf_accept, splitmerge_accept = splitmerge_accept)
  }

  # return as list
  return(ret)
}

#------------------------------------------------
# geweke_pvalue
# return p-value of Geweke's diagnostic convergence statistic, estimated from package coda
# (not exported)
#' @noRd
geweke_pvalue <- function(x) {
  ret <- 2*pnorm(abs(geweke.diag(x)$z), lower.tail=FALSE)
  return(ret)
}

#------------------------------------------------
# test convergence
# check that geweke p-value non-significant on values x[1:n]
# (not exported)
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
# (not exported)
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

#------------------------------------------------
# call_hungarian
# calls C++ implementation of the Hungarian algorithm for binding best matching
# in a linear sum assigment problem. This is function is used in testing.
# (not exported)
#' @noRd
call_hungarian <- function(x) {
  args <- list(cost_mat = mat_to_rcpp(x))
  call_hungarian_cpp(args)
}
