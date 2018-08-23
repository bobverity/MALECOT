
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib MALECOT
#' @import parallel
#' @import coda
#' @import ggplot2
#' @import gridExtra
#' @importFrom plotly plot_ly
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices colorRampPalette
#' @import stats
#' @import utils
#' @import graphics
NULL

#------------------------------------------------
#' @title Bind bi-allelic data to project
#'
#' @description Bind data in bi-allelic format to MALECOT project. Data should
#'   be formatted as a dataframe with samples in rows and loci in columns.
#'   Genetic data should be coded as 1 (homozygote REF allele), 0 (homozygote
#'   ALT allele), or 0.5 (heterozygote). Additional meta-data columns can be
#'   specified, including a column for sample IDs and a column for sampling
#'   population.
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param df a dataframe containing genetic information and optional meta-data
#' @param ID_col which column of the input data contains the sample IDs. If NULL
#'   then IDs must be defined seperately through the \code{ID} argument
#' @param pop_col which column of the input data contains the ostensible 
#'   population of the samples. If NULL then populations must be defined 
#'   seperately through the \code{pop} argument
#' @param data_cols which columns of the input data contain genetic information.
#'   Defaults to all remaining columns of the data once special columns have 
#'   been accounted for
#' @param ID sample IDs, if not using the \code{ID_col} option
#' @param pop ostensible populations, if not using the \code{pop_col} option
#' @param missing_data what value represents missing data (defaults to -9). Must
#'   be a positive or negative integer, and cannot equal 0 or 1 as these are 
#'   reserved for genetic data.
#' @param name optional name of the data set to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data
#'
#' @export
#' @examples
#' # TODO

bind_data_biallelic <- function(project, df, ID_col = 1, pop_col = NULL, data_cols = NULL, ID = NULL, pop = NULL, missing_data = -9, name = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  assert_dataframe(df)
  n <- nrow(df)
  if (!is.null(ID_col)) {
    assert_single_pos_int(ID_col, zero_allowed = FALSE)
    assert_leq(ID_col, ncol(df))
    ID <- as.character(df[,ID_col])
  }
  if (!is.null(pop_col)) {
    assert_single_pos_int(pop_col, zero_allowed = FALSE)
    assert_leq(pop_col, ncol(df))
    pop <- df[,pop_col]
  }
  if (is.null(data_cols)) {
    data_cols <- setdiff(1:ncol(df), c(ID_col, pop_col))
  } else {
    assert_pos_int(data_cols, zero_allowed = FALSE)
    assert_leq(data_cols, ncol(df))
    assert_noduplicates(data_cols)
  }
  ID <- define_default(ID, paste0("sample", zero_pad_simple(1:n, nchar(n))))
  assert_string(ID)
  assert_length(ID, n)
  pop <- define_default(pop, rep(1,n))
  assert_pos_int(pop)
  assert_length(pop,n)
  assert_single_int(missing_data)
  assert_not_in(missing_data, c(0, 1))
  if (!is.null(name)) {
    assert_single_string(name)
  }
  assert_single_logical(check_delete_output)
  
  # check before overwriting existing output
  if (project$active_set>0 && check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- malecot_project()
  }
  
  # extract genetic data
  dat <- as.matrix(df[, data_cols, drop = FALSE])
  L <- ncol(dat)
  locus_names <- colnames(dat)
  apply(dat, 1, function(x) {
    assert_in(x, c(0, 0.5, 1, missing_data), name_x = "genetic data")
  })
  
  # process genetic data
  w <- which(dat == missing_data)
  dat[dat == 1] <- 1
  dat[dat == 0.5] <- 2
  dat[dat == 0] <- 3
  dat[w] <- 0
  
  # create dat_processed list
  dat_processed <- list(data = dat,
                        n = nrow(dat),
                        L = ncol(dat),
                        ID = ID,
                        pop = pop,
                        alleles = rep(2,L),
                        locus_names = locus_names,
                        data_format = "biallelic",
                        name = name)
  
  # add data to project
  project$data <- df
  project$data_processed <- dat_processed
  
  return(project)
}

#------------------------------------------------
#' @title Bind multi-allelic format data to project
#'
#' @description Bind data in multi-allelic format to MALECOT project. Data 
#'   should be formatted as a dataframe with three columns: "sample_ID", "locus"
#'   and "haplotype". Each row of this dataframe specifies a haplotype that was
#'   observed at that locus in that individual (i.e. long format). Haplotypes
#'   should be coded as positive integers.
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param df a dataframe containing genetic information and optional meta-data
#' @param pop ostensible populations, if not using the \code{pop_col} option
#' @param missing_data what value represents missing data (defaults to -9). Must
#'   be a positive or negative integer, and cannot equal 0 or 1 as these are 
#'   reserved for genetic data.
#' @param name optional name of the data set to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data
#'
#' @export
#' @examples
#' # TODO

bind_data_multiallelic <- function(project, df, pop = NULL, missing_data = -9, name = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  assert_dataframe(df)
  pop <- define_default(pop, rep(1,n))
  assert_pos_int(pop)
  assert_length(pop,n)
  assert_single_int(missing_data)
  if (!is.null(name)) {
    assert_single_string(name)
  }
  assert_single_logical(check_delete_output)
  
  # check before overwriting existing output
  if (project$active_set>0 && check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- malecot_project()
  }
  
  # check data format
  assert_dataframe(df)
  assert_ncol(df, 3)
  assert_eq(names(df), c("sample_ID", "locus", "haplotype"))
  
  # check sample_ID column
  assert_string(df$sample_ID)
  ID <- unique(df$sample)
  n <- length(ID)
  
  # check locus column
  locus_names <- unique(subset(df, df$sample_ID == df$sample_ID[1])$locus)
  L <- length(locus_names)
  good_loci <- mapply(function(x) {
                        return(isTRUE(all.equal(unique(x), locus_names)))
                      }, split(df$locus, f = df$sample_ID))
  if (!all(good_loci)) {
    stop(sprintf("all samples must contain loci in the range 1:%s", L))
  }
  
  # check haplotype column
  if (any(df$haplotype <= 0 && df$haplotype != missing_data)) {
    stop("for the multi-allelic format haplotypes must be coded as positive integers (or as missing data)")
  }
  
  # process genetic data. Re-factor loci as increasing integers, and haplotypes
  # as increasing integers with 0 indicating missing data
  df_processed <- df
  df_processed$locus <- match(df_processed$locus, locus_names)
  df_processed$haplotype[df_processed$haplotype == missing_data] <- NA
  alleles <- mapply(function(x) {
                      length(unique(x[!is.na(x)]))
                    }, split(df_processed$haplotype, f = df_processed$locus))
  if (all(alleles == 2)) {
    message("Note: data contains two alleles at every locus. Consider reformatting in bi-allelic format to speed up MCMC")
  }
  new_haplotype_order <- mapply(function(x) {
                                  order(unique(x[!is.na(x)]))
                                }, split(df_processed$haplotype, f = df_processed$locus), SIMPLIFY = FALSE)
  new_haplotype <- mapply(function(x,y) {
                            new_haplotype_order[[x]][y]
                          }, df_processed$locus, df_processed$haplotype)
  df_processed$haplotype <- new_haplotype
  df_processed$haplotype[is.na(df_processed$haplotype)] <- 0
  
  # add data to project
  project$data <- df
  project$data_processed <- list(data = df_processed,
                                 n = n,
                                 L = L,
                                 ID = ID,
                                 pop = pop,
                                 alleles = alleles,
                                 locus_names = locus_names,
                                 data_format = "multiallelic",
                                 name = name)
  
  return(project)
}

#------------------------------------------------
#' @title Create new MALECOT parameter set
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param name the name of the parameter set
#' @param lambda shape parameter(s) governing the prior on allele frequencies. 
#'   This prior is a Beta distribution for the bi-allelic case or a Dirichlet 
#'   distribution for the multi-allelic case. \code{lambda} can be a single
#'   scalar value, in which case the same shape paremeter is applied to all loci
#'   and all alleles (i.e. a symmetric Beta or Dirichlet prior), or
#'   \code{lambda} can be a list of length \code{L} containing vectors of length
#'   equal to the number of alleles at each locus, allowing different shape
#'   parameters to be specified for each allele individually
#' @param COI_model the type of prior on COI. Must be one of "uniform",
#'   "poisson", or "nb" (negative binomial)
#' @param COI_max TODO
#' @param COI_manual TODO
#' @param estimate_COI_mean TODO
#' @param COI_mean TODO
#' @param COI_dispersion must be > 1
#' @param estimate_error TODO
#' @param e1 TODO
#' @param e2 TODO
#' @param e1_max TODO
#' @param e2_max TODO
#'
#' @export
#' @examples
#' # TODO

new_set <- function(project, name = "(no name)", lambda = 1.0, COI_model = "poisson", COI_max = 20, COI_manual = NULL, estimate_COI_mean = TRUE, COI_mean = 3, COI_dispersion = 1, estimate_error = FALSE, e1 = 0, e2 = 0, e1_max = 0.2, e2_max = 0.2) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  n <- project$data_processed$n
  L <- project$data_processed$L
  assert_single_string(name)
  assert_in(length(lambda), c(1,L))
  assert_pos(unlist(lambda), zero_allowed = FALSE)
  assert_single_string(COI_model)
  assert_in(COI_model, c("uniform", "poisson", "nb"))
  assert_single_pos_int(COI_max, zero_allowed = FALSE)
  COI_manual <- define_default(COI_manual, rep(-1,n))
  assert_vector(COI_manual)
  assert_length(COI_manual, n)
  assert_int(COI_manual)
  assert_bounded(COI_manual[COI_manual != -1], left = 1, right = COI_max, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_logical(estimate_COI_mean)
  if (estimate_COI_mean) {
    if (COI_model != "uniform") {
      assert_single_pos(COI_mean)
      assert_bounded(COI_mean, left = 1, right = COI_max, inclusive_left = TRUE, inclusive_right = TRUE)
      if (COI_model == "nb") {
        assert_single_pos(COI_dispersion)
        assert_gr(COI_dispersion, 1)
      }
    }
  }
  assert_single_logical(estimate_error)
  assert_single_numeric(e1)
  assert_bounded(e1, left = 0.0, right = 1.0, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_numeric(e2)
  assert_bounded(e2, left = 0.0, right = 1.0, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_bounded(e1_max)
  assert_bounded(e2_max)
  if (estimate_error) {
    assert_bounded(e1, right = e1_max)
    assert_bounded(e2, right = e2_max)
  }
  
  # check that lambda specified correctly
  if (length(lambda) == L) {
    lambda_length <- mapply(length, lambda)
    assert_eq(project$data_processed$alleles, lambda_length, message = sprintf("lambda does not match number of alleles at every locus. Length of lambda: %s. Number of alleles: %s.", nice_format(lambda_length), nice_format(project$data_processed$alleles)))
  } else {
    assert_length(lambda, 1)
  }
  
  # expand lambda to list
  if (length(lambda) == 1) {
    lambda <- mapply(rep, lambda, times = project$data_processed$alleles, SIMPLIFY = FALSE)
  }
  
  # count current parameter sets and add one
  s <- length(project$parameter_sets) + 1
  
  # make new set active
  project$active_set <- s
  
  # create new parameter set
  project$parameter_sets[[s]] <- list(name = name,
                                      lambda = lambda,
                                      COI_model = COI_model,
                                      COI_max = COI_max,
                                      COI_manual = COI_manual,
                                      estimate_COI_mean = estimate_COI_mean,
                                      COI_mean = COI_mean,
                                      COI_dispersion = COI_dispersion,
                                      estimate_error = estimate_error,
                                      e1 = e1,
                                      e2 = e2,
                                      e1_max = e1_max,
                                      e2_max = e2_max)
  
  names(project$parameter_sets)[s] <- paste0("set", s)
  
  # create new output corresponding to this set
  GTI_logevidence <- data.frame(K = numeric(),
                                mean = numeric(),
                                SE = numeric())
  class(GTI_logevidence) <- "malecot_GTI_logevidence"
  
  GTI_posterior <- data.frame(K = numeric(),
                              Q2.5 = numeric(),
                              Q50 = numeric(),
                              Q97.5 = numeric())
  class(GTI_posterior) <- "malecot_GTI_posterior"
  
  project$output$single_set[[s]] <- list(single_K = list(),
                                         all_K = list(GTI_logevidence = GTI_logevidence,
                                                      GTI_posterior = GTI_posterior))
  
  names(project$output$single_set) <- paste0("set", 1:length(project$output$single_set))
  
  # expand summary output over all parameter sets
  GTI_logevidence_model <- rbind(project$output$all_sets$GTI_logevidence_model, data.frame(set = s, name = name, mean = NA, SE = NA, stringsAsFactors = FALSE))
  class(GTI_logevidence_model) <- "malecot_GTI_logevidence_model"
  project$output$all_sets$GTI_logevidence_model <- GTI_logevidence_model
  
  GTI_posterior_model <- rbind(project$output$all_sets$GTI_posterior_model, data.frame(set = s, name = name, Q2.5 = NA, Q50 = NA, Q97.5 = NA, stringsAsFactors = FALSE))
  class(GTI_posterior_model) <- "malecot_GTI_posterior_model"
  project$output$all_sets$GTI_posterior_model <- GTI_posterior_model
  
  # return invisibly
  invisible(project)
}

#------------------------------------------------
#' @title Delete parameter set
#'   
#' @description Delete a given parameter set from a MALECOT project.
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param set which set to delete. Defaults to the current active set
#' @param check_delete_output whether to prompt the user before deleting any 
#'   existing output
#'
#' @export
#' @examples
#' # TODO

delete_set <- function(project, set = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  assert_single_logical(check_delete_output)
  
  # set index to active_set by default
  set <- define_default(set, project$active_set)
  
  # further checks
  assert_single_pos_int(set, zero_allowed = FALSE)
  assert_leq(set, length(project$parameter_sets))
  
  # check before overwriting existing output
  if (project$active_set>0 & check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no(sprintf("Any existing output for set %s will be deleted. Continue? (Y/N): ", set))) {
      return(project)
    }
  }
  
  # drop chosen parameter set
  project$parameter_sets[[set]] <- NULL
  
  # drop chosen output
  project$output$single_set[[set]] <- NULL
  
  GTI_logevidence_model <- as.data.frame(unclass(project$output$all_sets$GTI_logevidence_model))[-set,]
  class(GTI_logevidence_model) <- "malecot_GTI_logevidence_model"
  project$output$all_sets$GTI_logevidence_model <- GTI_logevidence_model
  
  GTI_posterior_model <- as.data.frame(unclass(project$output$all_sets$GTI_posterior_model))[-set,]
  class(GTI_posterior_model) <- "malecot_GTI_posterior_model"
  project$output$all_sets$GTI_posterior_model <- GTI_posterior_model
  
  # make new final set active
  project$active_set <- length(project$parameter_sets)
  
  # recalculate evidence over sets if needed
  if (project$active_set > 0) {
    project <- recalculate_evidence(project)
  }
  
  # return invisibly
  invisible(project)
}

#------------------------------------------------
#' @title Run main MCMC
#'
#' @description Run the main MALECOT MCMC. Model parameters are taken from the 
#'   current active parameter set, and MCMC parameters are passed in as 
#'   arguments. All output is stored within the project.
#'   
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K the values of K that the MCMC will explore
#' @param precision the level of precision at which allele frequencies are 
#'   represented in the bi-allelic case. This allows the use of look-up tables
#'   for the likelihood, which significantly speeds up the MCMC. Set to 0 to use
#'   exact values (up to C++ "double" precision) rather than using look-up 
#'   tables
#' @param burnin the number of burn-in iterations
#' @param samples the number of sampling iterations
#' @param rungs the number of temperature rungs
#' @param GTI_pow the power used in the generalised thermodynamic integration 
#'   method. Must be greater than 1.1
#' @param auto_converge whether convergence should be assessed automatically 
#'   every \code{converge_test} iterations, leading to termination of the 
#'   burn-in phase. If \code{FALSE} then the full \code{burnin} iterations are 
#'   used
#' @param converge_test test for convergence every \code{convergence_test} 
#'   iterations if \code{auto_converge} is being used
#' @param solve_label_switching_on whether to implement the Stevens' solution to
#'   the label-switching problem. If turned off then Q-matrix output will no 
#'   longer be correct, although evidence estimates will be unaffected.
#' @param coupling_on whether to implement Metropolis-coupling over temperature 
#'   rungs
#' @param cluster option to pass in a cluster environment (see package 
#'   "parallel")
#' @param pb_markdown whether to run progress bars in markdown mode, in which 
#'   case they are updated once at the end to avoid large amounts of output.
#' @param silent whether to suppress all console output
#'
#' @export
#' @examples
#' # TODO

run_mcmc <- function(project, K = NULL, precision = 0.01, burnin = 1e3, samples = 1e3, rungs = 10, GTI_pow = 3, auto_converge = TRUE, converge_test = ceiling(burnin/10), solve_label_switching_on = TRUE, coupling_on = TRUE, cluster = NULL, pb_markdown = FALSE, silent = !is.null(cluster)) {
  
  # start timer
  t0 <- Sys.time()
  
  # ---------- check inputs ----------
  
  assert_custom_class(project, "malecot_project")
  assert_pos_int(K, zero_allowed = FALSE)
  assert_single_pos(precision, zero_allowed = TRUE)
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(rungs, zero_allowed = FALSE)
  assert_single_pos(GTI_pow)
  assert_gr(GTI_pow, 1.1)
  assert_single_logical(auto_converge)
  assert_single_pos_int(converge_test, zero_allowed = FALSE)
  assert_single_logical(solve_label_switching_on)
  assert_single_logical(coupling_on)
  if (!is.null(cluster)) {
    assert_custom_class(project, "cluster")
  }
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # get active set
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # check that precision is either zero or leads to simple sequence
  if (precision != 0) {
    if(!((1/precision)%%1) == 0) {
      stop("1/precision must be an integer")
    }
  }
  
  # get useful quantities
  n <- project$data_processed$n
  L <- project$data_processed$L
  data_format <- project$data_processed$data_format
  dat <- project$data_processed$data
  
  # ---------- create argument lists ----------
  
  # data list
  args_data <- list(data = mat_to_rcpp(dat),
                    n = n,
                    L = L,
                    alleles = project$data_processed$alleles)
  
  # input arguments list
  args_inputs <- list(precision = precision,
                      burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      GTI_pow = GTI_pow,
                      auto_converge = auto_converge,
                      converge_test = converge_test,
                      solve_label_switching_on = solve_label_switching_on,
                      coupling_on = coupling_on,
                      pb_markdown = pb_markdown,
                      silent = silent)
  
  # combine model parameters list with input arguments
  args_model <- c(project$parameter_sets[[s]], args_inputs)
  
  # get COI_model in numeric form
  args_model$COI_model_numeric <- match(args_model$COI_model, c("uniform", "poisson" ,"nb"))
  
  # R functions to pass to Rcpp
  args_functions <- list(test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # define final argument list over all K
  parallel_args <- list()
  for (i in 1:length(K)) {
    
    # create progress bars
    pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
    pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
    args_progress <- list(pb_burnin = pb_burnin,
                          pb_samples = pb_samples)
    
    # incporporate arguments unique to this K
    args_model$K <- K[i]
    
    # create argument list
    parallel_args[[i]] <- list(args_data = args_data,
                               args_model = args_model,
                               args_functions = args_functions,
                               args_progress = args_progress)
  }
  
  # ---------- run MCMC ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) { # run in parallel
    clusterEvalQ(cluster, library(malecot))
    if (data_format == "biallelic") {
      output_raw <- clusterApplyLB(cl = cluster, parallel_args, run_mcmc_biallelic_cpp)
    }
    if (data_format == "multiallelic") {
      output_raw <- clusterApplyLB(cl = cluster, parallel_args, run_mcmc_multiallelic_cpp)
    }
  } else { # run in serial
    if (data_format == "biallelic") {
      output_raw <- lapply(parallel_args, run_mcmc_biallelic_cpp)
    }
    if (data_format == "multiallelic") {
      output_raw <- lapply(parallel_args, run_mcmc_multiallelic_cpp)
    }
  }
  
  # ---------- process results ----------
  
  # begin processing results
  if (!silent) {
    cat("Processing results\n")
  }
  
  # loop through K
  ret <- list()
  for (i in 1:length(K)) {
    
    # create name lists
    sample_names <- paste0("sample", 1:n)
    locus_names <- paste0("locus", 1:L)
    deme_names <- paste0("deme", 1:K[i])
    rung_names <- paste0("rung", 1:rungs)
    
    # ---------- raw mcmc results ----------
    
    # get loglikelihood in coda::mcmc format
    loglike_burnin <- mapply(function(x){mcmc(x)}, output_raw[[i]]$loglike_burnin)
    loglike_sampling <- mcmc(t(rcpp_to_mat(output_raw[[i]]$loglike_sampling)))
    
    # get full m trace in coda::mcmc format
    full_m <- mcmc(rcpp_to_mat(output_raw[[i]]$m_store))
    colnames(full_m) <- sample_names
    
    # get full p trace in coda::mcmc format
    full_p <- list()
    for (k in 1:K[i]) {
      full_p[[k]] <- list()
      for (j in 1:L) {
        p_j <- rcpp_to_mat(output_raw[[i]]$p_store[[k]][[j]])
        colnames(p_j) <- paste0("allele", 1:ncol(p_j))
        full_p[[k]][[j]] <- mcmc(p_j)
      }
      names(full_p[[k]]) <- locus_names
    }
    names(full_p) <- deme_names
    
    # get full e1 and e2 trace in coda::mcmc format
    full_e1 <- full_e2 <- NULL
    if (args_model$estimate_error) {
      full_e1 <- mcmc(output_raw[[i]]$e1_store)
      full_e2 <- mcmc(output_raw[[i]]$e2_store)
    }
    
    # get full COI_mean trace in coda::mcmc format
    full_COI_mean <- NULL
    if (args_model$COI_model %in% c("poisson", "nb")) {
      full_COI_mean <- rcpp_to_mat(output_raw[[i]]$COI_mean_store)
      colnames(full_COI_mean) <- deme_names
      full_COI_mean <- mcmc(full_COI_mean)
    }
    
    # ---------- summary results ----------
    
    # get quantiles over sampling loglikelihoods
    loglike_quantiles <- t(apply(loglike_sampling, 2, quantile_95))
    rownames(loglike_quantiles) <- rung_names
    class(loglike_quantiles) <- "malecot_loglike_quantiles"
    
    # get quantiles over m
    m_quantiles <- t(apply(full_m, 2, quantile_95))
    class(m_quantiles) <- "malecot_m_quantiles"
    
    # get quantiles over p
    p_quantiles <- list()
    for (k in 1:K[i]) {
      p_quantiles[[k]] <- list()
      for (j in 1:L) {
        p_quantiles[[k]][[j]] <- t(apply(full_p[[k]][[j]], 2, quantile_95))
      }
      names(p_quantiles[[k]]) <- locus_names
      class(p_quantiles[[k]]) <- "malecot_p_quantiles"
    }
    names(p_quantiles) <- deme_names
    
    # get quantiles over COI_mean
    if (is.null(full_COI_mean)) {
      fake_COI_mean <- matrix((args_model$COI_max + 1)/2, ncol = K[i])
      colnames(fake_COI_mean) <- deme_names
      COI_mean_quantiles <- t(apply(fake_COI_mean, 2, quantile_95))
    } else {
      COI_mean_quantiles <- t(apply(full_COI_mean, 2, quantile_95))
    }
    
    # get quantiles over e1 and e2
    if (args_model$estimate_error) {
      e_quantiles <- rbind(quantile_95(full_e1), quantile_95(full_e2))
    } else {
      e_quantiles <- rbind(quantile_95(args_model$e1), quantile_95(args_model$e2))
    }
    rownames(e_quantiles) <- c("e1", "e2")
    class(e_quantiles) <- "malecot_e_quantiles"
    
    # process Q-matrix
    qmatrix <- rcpp_to_mat(output_raw[[i]]$qmatrix)
    colnames(qmatrix) <- deme_names
    rownames(qmatrix) <- sample_names
    class(qmatrix) <- "malecot_qmatrix"
    
    # ---------- GTI path and model evidence ----------
    
    # get ESS
    ESS <- effectiveSize(loglike_sampling)
    ESS[ESS == 0] <- samples # if no variation then assume zero autocorrelation
    ESS[ESS > samples] <- samples # ESS cannot exceed actual number of samples taken
    names(ESS) <- rung_names
    
    # weight likelihood according to GTI_pow
    loglike_weighted <- loglike_sampling
    for (j in 1:rungs) {
      beta_j <- j/rungs
      loglike_weighted[,j] <- GTI_pow*beta_j^(GTI_pow-1) * loglike_sampling[,j]
    }
    
    # calculate GTI path mean and SE
    GTI_path_mean <- colMeans(loglike_weighted)
    GTI_path_var <- apply(loglike_weighted, 2, var)
    GTI_path_SE <- sqrt(GTI_path_var/ESS)
    GTI_path <- data.frame(mean = GTI_path_mean, SE = GTI_path_SE)
    rownames(GTI_path) <- rung_names
    class(GTI_path) <- "malecot_GTI_path"
    
    # calculate GTI estimate of log-evidence
    GTI_vec <- 0.5*loglike_weighted[,1]/rungs
    if (rungs>1) {
      for (j in 2:rungs) {
        GTI_vec <- GTI_vec + 0.5*(loglike_weighted[,j]+loglike_weighted[,j-1])/rungs
      }
    }
    GTI_logevidence_mean <- mean(GTI_vec)
    
    # calculate standard error of GTI estimate
    GTI_ESS <- as.numeric(effectiveSize(GTI_vec))
    if (GTI_ESS==0) {
      GTI_ESS <- samples # if no variation then assume perfect mixing
    }
    GTI_logevidence_SE <- sqrt(var(GTI_vec)/GTI_ESS)
    
    # produce final GTI_logevidence object
    GTI_logevidence <- data.frame(estimate = GTI_logevidence_mean,
                                  SE = GTI_logevidence_SE)
    
    # ---------- acceptance rates ----------
    
    # process acceptance rates
    p_accept <- rcpp_to_mat(output_raw[[i]]$p_accept)/samples
    e1_accept <- output_raw[[i]]$e1_accept/samples
    e2_accept <- output_raw[[i]]$e2_accept/samples
    coupling_accept <- output_raw[[i]]$coupling_accept/samples
    scaf_trials <- output_raw[[i]]$scaf_trials/samples
    scaf_accept <- output_raw[[i]]$scaf_accept/samples
    split_merge_accept <- output_raw[[i]]$split_merge_accept/samples
    
    # ---------- save arguments ----------
    
    output_args <- list(precision = precision,
                        burnin = burnin,
                        samples = samples,
                        rungs = rungs,
                        GTI_pow = GTI_pow,
                        auto_converge = auto_converge,
                        converge_test = converge_test,
                        solve_label_switching_on = solve_label_switching_on,
                        pb_markdown = pb_markdown,
                        silent = silent)
    
    # ---------- save results ----------
    
    # add to project
    project$output$single_set[[s]]$single_K[[K[i]]] <- list()
    
    project$output$single_set[[s]]$single_K[[K[i]]]$summary <- list(qmatrix = qmatrix,
                                                                    loglike_quantiles = loglike_quantiles,
                                                                    m_quantiles = m_quantiles,
                                                                    p_quantiles = p_quantiles,
                                                                    e_quantiles = e_quantiles,
                                                                    COI_mean_quantiles = COI_mean_quantiles,
                                                                    ESS = ESS,
                                                                    GTI_path = GTI_path,
                                                                    GTI_logevidence = GTI_logevidence)
    
    project$output$single_set[[s]]$single_K[[K[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                                                loglike_sampling = loglike_sampling,
                                                                m = full_m,
                                                                p = full_p,
                                                                e1 = full_e1,
                                                                e2 = full_e2,
                                                                COI_mean = full_COI_mean,
                                                                p_accept = p_accept,
                                                                e1_accept = e1_accept,
                                                                e2_accept = e2_accept,
                                                                coupling_accept = coupling_accept)
    
    project$output$single_set[[s]]$single_K[[K[i]]]$function_call <- list(args = output_args,
                                                                          call = match.call())
    
  } # end loop over K
  
  # name output over K
  K_all <- length(project$output$single_set[[s]]$single_K)
  names(project$output$single_set[[s]]$single_K) <- paste0("K", 1:K_all)
  
  # ---------- tidy up and end ----------
  
  # reorder qmatrices
  project <- align_qmatrix(project)
  
  # recalculate evidence over K
  project <- recalculate_evidence(project)
  
  # end timer
  tdiff <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (tdiff<60) {
    message(sprintf("Total run-time: %s seconds", round(tdiff, 2)))
  } else {
    message(sprintf("Total run-time: %s minutes", round(tdiff/60, 2)))
  }
  
  # return invisibly
  invisible(project)
}

#------------------------------------------------
# extract GTI_logevidence from all K within a given parameter set
#' @noRd
get_GTI_logevidence_K <- function(proj, s) {
  
  # extract objects of interest
  x <- proj$output$single_set[[s]]$single_K
  if (length(x)==0) {
    return(proj)
  }
  
  # get log-evidence over all K
  GTI_logevidence_raw <- mapply(function(y) {
    GTI_logevidence <- y$summary$GTI_logevidence
    if (is.null(GTI_logevidence)) {
      return(rep(NA,2))
    } else {
      return(unlist(GTI_logevidence))
    }
  }, x)
  GTI_logevidence <- as.data.frame(t(GTI_logevidence_raw))
  names(GTI_logevidence) <- c("mean", "SE")
  rownames(GTI_logevidence) <- NULL
  GTI_logevidence <- cbind(K = 1:nrow(GTI_logevidence), GTI_logevidence)
  class(GTI_logevidence) <- "malecot_GTI_logevidence"
  
  # save result in project
  proj$output$single_set[[s]]$all_K$GTI_logevidence <- GTI_logevidence
  
  # return modified project
  return(proj)
}

#------------------------------------------------
# compute posterior over several log-evidence estimates
#' @noRd
get_GTI_posterior <- function(x) {
  
  # return NULL if all NA
  if (length(x$mean)==0 || all(is.na(x$mean))) {
    return(NULL)
  }
  
  # produce posterior estimates by simulation
  w <- which(!is.na(x$mean))
  GTI_posterior_raw <- GTI_posterior_K_sim_cpp(list(mean = x$mean[w],
                                                    SE = x$SE[w],
                                                    reps = 1e6))$ret
  # get posterior quantiles in dataframe
  GTI_posterior_quantiles <- t(mapply(quantile_95, GTI_posterior_raw))
  GTI_posterior <- data.frame(Q2.5 = rep(NA, nrow(GTI_posterior_quantiles)), Q50 = NA, Q97.5 = NA)
  GTI_posterior[w,] <- GTI_posterior_quantiles
  
  return(GTI_posterior)
}

#------------------------------------------------
# call get_GTI_posterior over values of K
#' @noRd
get_GTI_posterior_K <- function(proj, s) {
  
  # calculate posterior K
  GTI_posterior <- get_GTI_posterior(proj$output$single_set[[s]]$all_K$GTI_logevidence)
  if (is.null(GTI_posterior)) {
    return(proj)
  }
  GTI_posterior <- cbind(K = 1:nrow(GTI_posterior), GTI_posterior)
  class(GTI_posterior) <- "malecot_GTI_posterior"
  proj$output$single_set[[s]]$all_K$GTI_posterior <- GTI_posterior
  
  # return modified project
  return(proj)
}

#------------------------------------------------
# call get_GTI_posterior over models
#' @noRd
get_GTI_posterior_model <- function(proj) {
  
  # calculate posterior model
  GTI_posterior_model_raw <- get_GTI_posterior(proj$output$all_sets$GTI_logevidence_model)
  if (is.null(GTI_posterior_model_raw)) {
    return(proj)
  }
  proj$output$all_sets$GTI_posterior_model$Q2.5 <- GTI_posterior_model_raw$Q2.5
  proj$output$all_sets$GTI_posterior_model$Q50 <- GTI_posterior_model_raw$Q50
  proj$output$all_sets$GTI_posterior_model$Q97.5 <- GTI_posterior_model_raw$Q97.5
  
  # return modified project
  return(proj)
}

#------------------------------------------------
# integrate multiple log-evidence estimates by simulation
#' @noRd
integrate_GTI_logevidence <- function(x) {
  
  # return NULL if all NA
  if (length(x$mean)==0 || all(is.na(x$mean))) {
    return(NULL)
  }
  
  # produce integrated estimates by simulation
  w <- which(!is.na(x$mean))
  if (length(w)==1) {
    ret <- list(mean = x$mean[w], SE = x$SE[w])
  } else {
    ret <- GTI_integrated_K_sim_cpp(list(mean = x$mean[w], SE = x$SE[w], reps = 1e6))
  }
  
  return(ret)
}

#------------------------------------------------
# log-evidence estimates over K
#' @noRd
integrate_GTI_logevidence_K <- function(proj, s) {
  
  # integrate over K
  integrated_raw <- integrate_GTI_logevidence(proj$output$single_set[[s]]$all_K$GTI_logevidence)
  if (is.null(integrated_raw)) {
    return(proj)
  }
  proj$output$all_sets$GTI_logevidence_model$mean[s] <- integrated_raw$mean
  proj$output$all_sets$GTI_logevidence_model$SE[s] <- integrated_raw$SE
  
  # return modified project
  return(proj)
}

#------------------------------------------------
# align qmatrices over all K
#' @noRd
align_qmatrix <- function(proj) {
  
  # get active set
  s <- proj$active_set
  
  # extract objects of interest
  x <- proj$output$single_set[[s]]$single_K
  
  # find values with output
  null_output <- mapply(function(y) {is.null(y$summary$qmatrix)}, x)
  w <- which(!null_output)
  
  # set template to first qmatrix
  template_qmatrix <- x[[w[1]]]$summary$qmatrix
  n <- nrow(template_qmatrix)
  c <- ncol(template_qmatrix)
  
  # loop through output
  best_perm <- NULL
  for (i in w) {
    
    # expand template
    qmatrix <- unclass(x[[i]]$summary$qmatrix)
    template_qmatrix <- cbind(template_qmatrix, matrix(0, n, i-c))
    
    # calculate cost matrix
    cost_mat <- matrix(0,i,i)
    for (k1 in 1:i) {
      for (k2 in 1:i) {
        cost_mat[k1,k2] <- sum(qmatrix[,k1] * (log(qmatrix[,k1]+1e-100) - log(template_qmatrix[,k2]+1e-100)))
      }
    }
    
    # get lowest cost permutation
    best_perm <- call_hungarian(cost_mat)$best_matching + 1
    best_perm_order <- order(best_perm)
    
    # reorder elements
    deme_names <- paste0("deme", 1:ncol(qmatrix))
    qmatrix <- qmatrix[, best_perm_order, drop = FALSE]
    colnames(qmatrix) <- deme_names
    
    p_raw <- x[[i]]$raw$p[best_perm_order]
    names(p_raw) <- deme_names
    proj$output$single_set[[s]]$single_K[[i]]$raw$p <- p_raw
    
    p_quantiles <- x[[i]]$summary$p_quantiles[best_perm_order]
    names(p_quantiles) <- deme_names
    proj$output$single_set[[s]]$single_K[[i]]$summary$p_quantiles <- p_quantiles
    
    COI_mean_raw <- x[[i]]$raw$COI_mean
    if (!is.null(COI_mean_raw)) {
      COI_mean_raw <- COI_mean_raw[, best_perm_order, drop = FALSE]
      colnames(COI_mean_raw) <- deme_names
      proj$output$single_set[[s]]$single_K[[i]]$raw$COI_mean <- COI_mean_raw
    }
    
    # qmatrix becomes template for next level up
    template_qmatrix <- qmatrix
    
    # store result
    class(qmatrix) <- "malecot_qmatrix"
    proj$output$single_set[[s]]$single_K[[i]]$summary$qmatrix <- qmatrix
  }
  
  # return modified project
  return(proj)
}

#------------------------------------------------
#' @title Recalculate evidence and posterior estimates
#'
#' @description When a new value of K is added in to the analysis it affects all downstream evidence estimates that depend on this K - for example the overall model evidence integrated over K. This function therefore looks through all values of K in the active set and recalculates all downstream elements as needed.
#'
#' @param project a MALCOT project, as produced by the function 
#'   \code{malecot_project()}
#' 
#' @export

recalculate_evidence <- function(project) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  
  # get active set
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # get log-evidence over all K and load into project
  project <- get_GTI_logevidence_K(project, s)
  
  # produce posterior estimates of K by simulation and load into project
  project <- get_GTI_posterior_K(project, s)
  
  # get log-evidence over all parameter sets
  project <- integrate_GTI_logevidence_K(project, s)
  
  # get posterior over all parameter sets
  project <- get_GTI_posterior_model(project)
  
  # return modified project
  return(project)
}

#------------------------------------------------
#' @title Match groupings using Hungarian algorithm
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param proj TODO
#' @param target_group TODO
#'
#' @export
#' @examples
#' # TODO

fix_labels <- function(proj, target_group) {
  
  # TODO - checks on inputs
  
  # create target matrix from group
  n <- length(target_group)
  target_mat <- matrix(0, n, max(target_group))
  target_mat[cbind(1:n, target_group)] <- 1
  
  # get active set and check non-zero
  s <- proj$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # check there is output for active set
  if (is.null(proj$output[[s]])) {
    stop("currently no output corresponding to active set")
  }
  
  # fix labels in all output for this set
  for (i in 1:length(proj$output[[s]])) {
    
    # extract Qmatrix
    qmatrix <- unclass(proj$output[[s]][[i]]$summary$qmatrix)
    K <- ncol(qmatrix)
    
    # check same number of rows
    if (nrow(qmatrix)!=nrow(target_mat)) {
      stop("target_group must have same number of elements as qmatrix output")
    }
    
    # expand qmatrix and/or target_mat to get same number of cols
    if (K<ncol(target_mat)) {
      qmatrix <- cbind( qmatrix, matrix(0, nrow(qmatrix), ncol(target_mat)-K) )
    }
    if (K>ncol(target_mat)) {
      target_mat <- cbind( target_mat, matrix(0, nrow(qmatrix), K-ncol(target_mat)) )
    }
    
    # calculate cost matrix
    cost_mat <- matrix(0, ncol(qmatrix), ncol(qmatrix))
    for (k1 in 1:ncol(cost_mat)) {
      for (k2 in 1:ncol(cost_mat)) {
        cost_mat[k1,k2] <- sum(qmatrix[,k1]*(1-target_mat[,k2]))
      }
    }
    
    # run Hungarian algorithm in Rcpp
    best_perm <- fix_labels_cpp(list(cost_mat = mat_to_rcpp(cost_mat)))$best_perm
    best_perm_order <- order(best_perm)
    
    # limit best_perm_order to K
    best_perm_order <- best_perm_order[best_perm_order<=K]
    
    # re-order qmatrix
    qmatrix <- qmatrix[,best_perm_order, drop=FALSE]
    colnames(qmatrix) <- paste0("deme",1:ncol(qmatrix))
    class(qmatrix) <- "malecot_qmatrix"
    proj$output[[s]][[i]]$summary$qmatrix <- qmatrix
    
    # re-order allele frequencies
    p_quantiles <- proj$output[[s]][[i]]$summary$p_quantiles
    old_names <- names(p_quantiles)
    proj$output[[s]][[i]]$summary$p_quantiles <- p_quantiles[best_perm_order]
    names(proj$output[[s]][[i]]$summary$p_quantiles) <- old_names
    
    # TODO - re-order COI_mean
    #COI_mean_quantiles <- proj[["output"]][[s]][[i]]$summary$COI_mean_quantiles
    #oldNames <- colnames(COI_mean_quantiles)
    #COI_mean_quantiles <- COI_mean_quantiles[,bestPermOrder,drop=FALSE]
    #colnames(COI_mean_quantiles) <- oldNames
    #class(COI_mean_quantiles) <- "malProject_COI_mean_quantiles"
    #proj[["output"]][[s]][[i]]$summary$COI_mean_quantiles <- COI_mean_quantiles
    
    #TODO - reorder full output
    
  }
  
  # return invisibly
  invisible(proj)
}


