
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib MALECOT
#' @import parallel
#' @import coda
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
  assert_length(ID,n)
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
  n <- nrow(dat)
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
  locus_names <- unique(subset(df, sample_ID == sample_ID[1])$locus)
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
    lambda <- mapply(rep, lambda, times = project$data$alleles, SIMPLIFY = FALSE)
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
#' @title Generate scaffolds
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param proj TODO
#' @param scaffold_n TODO
#' @param K the number of mixture components
#' @param burnin the number of burn-in iterations
#' @param auto_converge whether burn-in should be diagnosed automatically
#' @param split_merge_on whether to implement a split-merge proposal
#' @param precision the level of precision that allele frequencies are represented at
#' @param silent supresses output messages
#' @param cluster pass in a cluster environment
#'
#' @export
#' @examples
#' # TODO

generate_scaffolds <- function(proj, scaffold_n = 10, K = NULL, burnin = 1e3, auto_converge = TRUE, split_merge_on = TRUE, precision = 0.01, silent = FALSE, cluster = NULL) {
  
  # ---------- check inputs ----------
  
  # check input arguments
  assert_custom_class(proj, "malecot_project")
  assert_pos_int(scaffold_n, zero_allowed = FALSE)
  assert_pos_int(burnin)
  assert_greq(burnin, 10)
  
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
                                   samples = 1,
                                   rungs = 1,
                                   auto_converge = auto_converge,
                                   coupling_on = FALSE,
                                   scaffold_on = FALSE,
                                   scaffold_n = scaffold_n,
                                   split_merge_on = split_merge_on,
                                   solve_label_switching_on = FALSE,
                                   precision = precision,
                                   GTI_pow = 1,
                                   silent = silent,
                                   output_format = 1,
                                   cluster = cluster))
  
  # get COI_model in numeric form
  args_model$COI_model_numeric <- match(args_model$COI_model, c("uniform", "poisson" ,"nb"))
  
  # R functions to pass to Rcpp
  args_functions <- list(test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # define final argument list over all K
  parallel_args <- list()
  for (i in 1:length(K)) {
    
    # create progress bars
    pb_scaf <- txtProgressBar(min = 0, max = scaffold_n, initial = NA, style = 3)
    args_pb <- list(pb_scaf = pb_scaf)
    
    # list of arguments unique to this K
    args_this_K <- list(K = K[i],
                        scaffold_group = list(),
                        scaffold_group_n = 0)
    
    # create argument list
    parallel_args[[i]] <- c(args_this_K, args_data, args_model, args_functions, args_pb)
  }
  
  
  # ---------- generate scaffolds ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) {  # run in parallel
    if (!inherits(cluster, "cluster")) {
      stop("expected a cluster object")
    }
    clusterEvalQ(cluster, library(MALECOT))
    if (data_format == "biallelic") {
      output_raw <- clusterApplyLB(cl = cluster, parallel_args, generate_scaffolds_biallelic_cpp)
    }
    if (data_format == "multiallelic") {
      output_raw <- clusterApplyLB(cl = cluster, parallel_args, generate_scaffolds_multiallelic_cpp)
    }
  } else {  # run in serial
    if (data_format == "biallelic") {
      output_raw <- lapply(parallel_args, generate_scaffolds_biallelic_cpp)
    }
    if (data_format == "multiallelic") {
      output_raw <- lapply(parallel_args, generate_scaffolds_multiallelic_cpp)
    }
  }
  
  
  # ---------- process results ----------
  
  # begin processing results
  if (!silent) {
    message("Processing results\n")
  }
  
  # loop through K
  ret <- list()
  for (i in 1:length(K)) {
    
    # get scaffold groups and loglikelihoods
    scaffold_group_new = rcpp_to_mat(output_raw[[i]]$scaffold_group)
    scaffold_loglike_new = output_raw[[i]]$scaffold_loglike
    
    # initialise or append to scaffold storage
    if (is.null(proj$scaffolds[[s]][[i]])) {
      scaffold_group <- scaffold_group_new
      scaffold_loglike <- scaffold_loglike_new
    } else {
      scaffold_group <- rbind(proj$scaffolds[[s]][[i]]$group, scaffold_group_new)
      scaffold_loglike <- c(proj$scaffolds[[s]][[i]]$loglike, scaffold_loglike_new)
    }
    
    # remove duplicates
    dup <- duplicated(scaffold_group)
    scaffold_loglike <- scaffold_loglike[!dup]
    scaffold_group <- scaffold_group[!dup,,drop = FALSE]
    
    # re-order to be decreasing in likelihood
    scaffold_group <- scaffold_group[order(scaffold_loglike),,drop = FALSE]
    scaffold_loglike <- scaffold_loglike[order(scaffold_loglike)]
    
    # add to project
    proj$scaffolds[[s]][[i]]$group <- scaffold_group
    proj$scaffolds[[s]][[i]]$loglike <- scaffold_loglike
  }
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Run MALECOT MCMC
#'
#' @description TODO
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

run_mcmc <- function(proj, K = NULL, burnin = 1e3, samples = 1e3, rungs = 1, auto_converge = TRUE, coupling_on = TRUE, scaffold_on = TRUE, split_merge_on = TRUE, solve_label_switching_on = TRUE, precision = 0.01, GTI_pow = 2, silent = FALSE, output_format = 1, cluster = NULL) {
  
  # ---------- check inputs ----------
  
  # check input arguments
  assert_custom_class(proj, "malecot_project")
  assert_pos_int(K)
  assert_pos_int(burnin)
  assert_greq(burnin, 0)
  assert_pos_int(samples, zero_allowed = FALSE)
  assert_pos_int(rungs, zero_allowed = FALSE)
  assert_bounded(GTI_pow, left = 2, right = 10)
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
  K_which <- match(K, proj$parameter_sets[[s]]$K_range)
  
  # ---------- process data ----------
  
  # if bi-allelic
  data <- proj$data$raw
  n <- proj$data$n
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
                                   scaffold_n = 0,
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
  args_functions <- list(test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # define final argument list over all K
  parallel_args <- list()
  for (i in 1:length(K)) {
    
    # create progress bars
    pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
    pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
    args_pb <- list(pb_burnin = pb_burnin,
                    pb_samples = pb_samples)
    
    # load scaffold groupings where possible
    scaffold_group <- proj$scaffolds[[s]][[i]]$group
    if (is.null(scaffold_group)) {
      scaffold_group <- list()
    } else {
      scaffold_group <- mat_to_rcpp(scaffold_group)
    }
    
    # list of arguments unique to this K
    args_this_K <- list(K = K[i],
                        scaffold_group = scaffold_group,
                        scaffold_group_n = length(scaffold_group))
    
    # create argument list
    parallel_args[[i]] <- c(args_this_K, args_data, args_model, args_functions, args_pb)
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
  
  # ---------- process results ----------
  
  # begin processing results
  message("Processing results\n")
  
  # loop through K
  ret <- list()
  for (i in 1:length(K)) {
    
    # create name lists
    ind_names <- paste0("ind", 1:n)
    locus_names <- paste0("locus", 1:L)
    deme_names <- paste0("deme", 1:K[i])
    rung_names <- paste0("rung", 1:rungs)
    
    # ---------- raw mcmc results ----------
    
    # get loglikelihood in coda::mcmc format
    loglike_burnin <- mcmc(t(rcpp_to_mat(output_raw[[i]]$burnin_loglike)))
    loglike_sampling <- mcmc(t(rcpp_to_mat(output_raw[[i]]$sampling_loglike)))
    
    # get full m trace in coda::mcmc format
    full_m <- mcmc(rcpp_to_mat(output_raw[[i]]$m_store))
    colnames(full_m) <- ind_names
    
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
    
    # get quantiles over e1 and e2
    e_quantiles <- NULL
    if (args_model$estimate_error) {
      e_quantiles <- rbind(quantile_95(full_e1), quantile_95(full_e2))
      rownames(e_quantiles) <- c("e1", "e2")
    }
    
    # process Q-matrix
    q_matrix <- rcpp_to_mat(output_raw[[i]]$q_matrix)
    colnames(q_matrix) <- deme_names
    rownames(q_matrix) <- ind_names
    class(q_matrix) <- "malecot_q_matrix"
    
    # ---------- GTI path and model evidence ----------
    
    # only possible if more than 1 temperature rung
    ESS <- GTI_path <- GTI_logevidence <- NULL
    if (rungs>1) {
      
      # get ESS
      ESS <- effectiveSize(loglike_sampling)
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
      for (j in 2:rungs) {
        GTI_vec <- GTI_vec + 0.5*(loglike_weighted[,j]+loglike_weighted[,j-1])/rungs
      }
      GTI_logevidence_mean <- mean(GTI_vec)
      
      # calculate standard error of GTI estimate
      GTI_ESS <- as.numeric(effectiveSize(GTI_vec))
      GTI_logevidence_SE <- sqrt(var(GTI_vec)/GTI_ESS)
      
      # produce final GTI_logevidence object
      GTI_logevidence <- data.frame(estimate = GTI_logevidence_mean,
                                    SE = GTI_logevidence_SE)
    }
    
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
    
    output_args <- list(burnin = burnin,
                        samples = samples,
                        rungs = rungs,
                        auto_converge = auto_converge,
                        coupling_on = coupling_on,
                        scaffold_on = scaffold_on,
                        split_merge_on = split_merge_on,
                        solve_label_switching_on = solve_label_switching_on,
                        precision = precision,
                        GTI_pow = GTI_pow,
                        silent = silent,
                        output_format = output_format,
                        cluster = cluster)
    
    # ---------- save results ----------
    
    # add to project
    proj$output[[s]][[K_which[i]]]$summary <- list(loglike_quantiles = loglike_quantiles,
                                          m_quantiles = m_quantiles,
                                          p_quantiles = p_quantiles,
                                          e_quantiles = e_quantiles,
                                          q_matrix = q_matrix,
                                          ESS = ESS,
                                          GTI_path = GTI_path,
                                          GTI_logevidence = GTI_logevidence)
    
    proj$output[[s]][[K_which[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                      loglike_sampling = loglike_sampling,
                                      full_m = full_m,
                                      full_p = full_p,
                                      full_e1 = full_e1,
                                      full_e2 = full_e2,
                                      full_COI_mean = full_COI_mean,
                                      p_accept = p_accept,
                                      e1_accept = e1_accept,
                                      e2_accept = e2_accept,
                                      coupling_accept = coupling_accept,
                                      scaf_trials = scaf_trials,
                                      scaf_accept = scaf_accept,
                                      split_merge_accept = split_merge_accept)
    
    proj$output[[s]][[K_which[i]]]$function_call <- list(args = output_args,
                                                call = match.call())
  }
  
  # ---------- recalculate evidence ----------
  
  # get logevidence over all K
  GTI_logevidence <- mapply(function(x){
    GTI_logevidence <- x$summary$GTI_logevidence
    if (is.null(GTI_logevidence)) {
      return(rep(FALSE,3))
    } else {
      return(c(TRUE, GTI_logevidence))
    }
  }, proj$output[[s]])
  
  # if there are logevidence estimates
  if (any(unlist(GTI_logevidence[1,]))) {
    
    # produce raw evidence draws by simulation
    w <- which(unlist(GTI_logevidence[1,]))
    GTI_evidence_raw <- GTI_evidence_sim_cpp(list(mean = unlist(GTI_logevidence[2,w]),
                                                  SE = unlist(GTI_logevidence[3,w]),
                                                  reps = 1e6))
    
    # get quantiles and load back into output
    for (i in 1:length(w)) {
      GTI_evidence <- quantile_95(GTI_evidence_raw$ret[[i]])
      proj$output[[s]][[w[i]]]$summary$GTI_evidence <- GTI_evidence
    }
  }
  
  # return invisibly
  invisible(proj)
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
    q_matrix <- unclass(proj$output[[s]][[i]]$summary$q_matrix)
    K <- ncol(q_matrix)
    
    # check same number of rows
    if (nrow(q_matrix)!=nrow(target_mat)) {
      stop("target_group must have same number of elements as q_matrix output")
    }
    
    # expand q_matrix and/or target_mat to get same number of cols
    if (K<ncol(target_mat)) {
      q_matrix <- cbind( q_matrix, matrix(0, nrow(q_matrix), ncol(target_mat)-K) )
    }
    if (K>ncol(target_mat)) {
      target_mat <- cbind( target_mat, matrix(0, nrow(q_matrix), K-ncol(target_mat)) )
    }
    
    # calculate cost matrix
    cost_mat <- matrix(0, ncol(q_matrix), ncol(q_matrix))
    for (k1 in 1:ncol(cost_mat)) {
      for (k2 in 1:ncol(cost_mat)) {
        cost_mat[k1,k2] <- sum(q_matrix[,k1]*(1-target_mat[,k2]))
      }
    }
    
    # run Hungarian algorithm in Rcpp
    best_perm <- fix_labels_cpp(list(cost_mat = mat_to_rcpp(cost_mat)))$best_perm
    best_perm_order <- order(best_perm)
    
    # limit best_perm_order to K
    best_perm_order <- best_perm_order[best_perm_order<=K]
    
    # re-order q_matrix
    q_matrix <- q_matrix[,best_perm_order, drop=FALSE]
    colnames(q_matrix) <- paste0("deme",1:ncol(q_matrix))
    class(q_matrix) <- "malecot_q_matrix"
    proj$output[[s]][[i]]$summary$q_matrix <- q_matrix
    
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


