
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib MALECOT
#' @import parallel
#' @import coda
#' @import ggplot2
#' @import gridExtra
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices colorRampPalette grey
#' @import stats
#' @import utils
#' @import graphics
NULL

#------------------------------------------------
#' @title Check that MALECOT package has loaded successfully
#'
#' @description Simple function to check that MALECOT package has loaded 
#'   successfully. Prints "MALECOT loaded successfully!" if so.
#'
#' @export

check_MALECOT_loaded <- function() {
  message("MALECOT loaded successfully!")
}

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
#'   Defaults to all remaining columns of the data once meta-data columns have 
#'   been accounted for
#' @param ID sample IDs. Ignored if using the \code{ID_col} option
#' @param pop ostensible populations. Ignored if using the \code{pop_col} option
#' @param missing_data what value represents missing data. Defaults to -9. Must 
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
  
  # check for invariant loci
  invariant <- apply(dat, 2, function(x) {
                      all(x == 1) || all(x == 3)
                    })
  if (any(invariant)) {
    stop("data cannot contain invariant loci (i.e. loci for which the same single haplotype is observed in every sample)")
  }
  
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
#'   observed at that locus in that individual. Haplotypes should be coded as
#'   positive integers.
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param df a dataframe with three columns, as decribed above
#' @param pop ostensible populations of the samples
#' @param missing_data what value represents missing data. Defaults to -9. Must
#'   be a positive or negative integer
#' @param alleles the number of alleles at each locus. If scalar then the same 
#'   number of alleles is assumed at all loci. If NULL then the number of
#'   alleles is inferred directly from data as the maximum observed value per
#'   locus
#' @param name optional name of the data set to aid in record keeping
#' @param check_delete_output whether to prompt the user before overwriting 
#'   existing data
#'
#' @export
#' @examples
#' # TODO

bind_data_multiallelic <- function(project, df, pop = NULL, missing_data = -9, alleles = NULL, name = NULL, check_delete_output = TRUE) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  assert_dataframe(df)
  n <- nrow(df)
  pop <- define_default(pop, rep(1,n))
  assert_pos_int(pop)
  assert_length(pop,n)
  assert_single_int(missing_data)
  if (!is.null(alleles)) {
    assert_pos_int(alleles, zero_allowed = FALSE)
  }
  if (!is.null(name)) {
    assert_single_string(name)
  }
  assert_single_logical(check_delete_output)
  
  # check before overwriting existing output
  if (project$active_set > 0 && check_delete_output) {
    
    # ask before overwriting. On abort, return original project
    if (!user_yes_no("All existing output and parameter sets for this project will be lost. Continue? (Y/N): ")) {
      return(project)
    }
    
    # replace old project with fresh empty version
    project <- malecot_project()
  }
  
  # check data format
  assert_ncol(df, 3)
  assert_eq(names(df), c("sample_ID", "locus", "haplotype"))
  
  # check sample_ID column
  assert_string(df$sample_ID)
  ID <- unique(df$sample)
  n <- length(ID)
  
  # check locus column
  assert_pos_int(df$locus)
  locus_names <- unique(subset(df, df$sample_ID == df$sample_ID[1])$locus)
  L <- length(locus_names)
  good_loci <- mapply(function(x) {
                        return(isTRUE(all.equal(unique(x), locus_names)))
                      }, split(df$locus, f = df$sample_ID))
  if (!all(good_loci)) {
    stop(sprintf("all samples must contain loci in the same range, i.e. 1:%s", L))
  }
  
  # expand alleles to vector if specified as scalar
  if (!is.null(alleles)) {
    if (length(alleles) == 1) {
      alleles <- rep(alleles, L)
    }
  }
  
  # check haplotype column
  if (any(df$haplotype <= 0 & df$haplotype != missing_data)) {
    stop("for the multi-allelic format haplotypes must be coded as positive integers or as missing data")
  }
  observed_n_alleles <- mapply(function(x) {
                                  length(unique(x[!is.na(x)]))
                                }, split(df$haplotype, f = df$locus))
  if (any(observed_n_alleles == 1)) {
    stop("data cannot contain invariant loci (i.e. loci for which the same single haplotype is observed in every sample)")
  }
  if (all(observed_n_alleles == 2)) {
    message("Note: data contains two alleles at every locus. Consider reformatting in bi-allelic format to speed up MCMC")
  }
  observed_max_alleles <- mapply(function(x) {
                                    max(x[!is.na(x)])
                                  }, split(df$haplotype, f = df$locus))
  if (is.null(alleles)) {
    alleles <- observed_max_alleles
  }
  if (any(observed_max_alleles > alleles)) {
    stop("observed number of alleles exceeds 'alleles' argument at some loci")
  }
  
  # process genetic data. Re-factor loci as increasing integers, and haplotypes
  # as increasing integers with 0 indicating missing data
  df_processed <- df
  df_processed$locus <- match(df_processed$locus, locus_names)
  df_processed$haplotype[df_processed$haplotype == missing_data] <- NA
  haplotype_uniques <- mapply(function(x) {
                                sort(unique(x[!is.na(x)]))
                              }, split(df_processed$haplotype, f = df_processed$locus), SIMPLIFY = FALSE)
  new_haplotype <- mapply(function(x,y) {
                            match(y, haplotype_uniques[[x]])
                          }, df_processed$locus, df_processed$haplotype)
  df_processed$haplotype <- new_haplotype
  df_processed$haplotype[is.na(df_processed$haplotype)] <- 0
  
  # reformat into list over samples and loci
  data_processed <- mapply(function(x) {
    tmp <- mapply(function(y) y$haplotype, split(x, f = x$locus), SIMPLIFY = FALSE)
    names(tmp) <- NULL
    tmp
  }, split(df_processed, f = df_processed$sample_ID), SIMPLIFY = FALSE)
  names(data_processed) <- NULL
  
  # add data to project
  project$data <- df
  project$data_processed <- list(data = data_processed,
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
#' @param lambda the shape parameter(s) of the prior on allele frequencies. This
#'   prior is Beta in the bi-allelic case, and Dirichlet in the multi-allelic 
#'   case. \code{lambda} can be:
#'   \itemize{
#'     \item{a single scalar value, in which case the same value is used for 
#'     every allele and every locus (i.e. the prior is symmetric)}
#'     \item{a vector of values, in which case the same vector is used for every
#'     locus. Only works if the same number of alleles applies at every locus}
#'     \item{a list of vectors specifying the shape parameter separately for
#'     each allele of each locus. The list must of length \code{L}, and must
#'     contain vectors of length equal to the number of alleles at that locus}
#'   }
#' @param COI_model the type of prior on COI. Must be one of "uniform",
#'   "poisson", or "nb" (negative binomial)
#' @param COI_max the maximum COI allowed for any given sample
#' @param COI_manual A vector of length n (where n is the number of samples)
#'   allowing the COI to be specified manually. Positive values indicate fixed
#'   COIs that should not be updated as part of the MCMC, while -1 values
#'   indicate that COIs should be estimated. Defaults to \code{rep(-1,n)},
#'   meaning all COIs will be esimated
#' @param estimate_COI_mean whether the mean COI should be estimated for each 
#'   subpopulation as part of the MCMC, otherwise the value \code{COI_mean} is 
#'   used for all subpopulations. Defaults to \code{TRUE}. Note that mean COI 
#'   estimation is only possible under the Poisson and negative binomial models
#'   (see \code{COI_model})
#' @param COI_mean single scalar value specifying the mean COI for all
#'   subpopulations (see \code{estimate_COI_mean} above)
#' @param COI_dispersion  the ratio of the variance to the mean of the prior on
#'   COI. Only applies under the negative binomial model. Must be >1, as a ratio
#'   of 1 can be achieved by using the Poisson distribution
#' @param estimate_error whether to estimate error probabilities \code{e1} and
#'   \code{e2}
#' @param e1 the probability of a true homozygote being incorrectly called as a
#'   heterozygote
#' @param e2 the probability of a true heterozygote being incorrectly called as a
#'   homozygote
#' @param e1_max the maximum possible value of \code{e1}
#' @param e2_max the maximum possible value of \code{e2}
#'
#' @export
#' @examples
#' # TODO

new_set <- function(project, name = "(no name)", lambda = 1.0, COI_model = "poisson", COI_max = 20, COI_manual = NULL, estimate_COI_mean = TRUE, COI_mean = 3.0, COI_dispersion = 2.0, estimate_error = FALSE, e1 = 0.0, e2 = 0.0, e1_max = 0.2, e2_max = 0.2) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  n <- project$data_processed$n
  L <- project$data_processed$L
  data_format <- project$data_processed$data_format
  assert_single_string(name)
  
  assert_pos(unlist(lambda), zero_allowed = FALSE)
  if (is.list(lambda)) {
    assert_length(lambda, L)
    lambda_nalleles <- mapply(length, lambda)
    if (!all(lambda_nalleles == project$data_processed$alleles)) {
      stop("when lambda is a list it must contain one entry per allele at every locus")
    }
  } else {
    if (length(lambda) > 1) {
      if (!all(project$data_processed$alleles == length(lambda))) {
        stop("when lambda is a vector its length must equal the number of alleles at every locus")
      }
    }
  }
  
  if (is.list(lambda)) {
    assert_length(lambda,L)
    assert_pos(unlist(lambda), zero_allowed = FALSE)
  } else {
    assert_pos(lambda, zero_allowed = FALSE)
  }
  
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
      assert_gr(COI_mean, 1)
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
  if (is.list(lambda)) {
    lambda_length <- mapply(length, lambda)
    assert_eq(project$data_processed$alleles, lambda_length, message = sprintf("when in list form, lambda must have one value for every allele at every locus."))
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
  
  # name parameter set
  names(project$parameter_sets)[s] <- paste0("set", s)
  
  # create new output at all_K level
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
  
  # name new output
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
#' @title Change the active set of a MALECOT project
#'
#' @description Change the active set of a MALECOT project
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param set the new active set
#'
#' @export
#' @examples
#' # TODO

active_set <- function(project, set) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  assert_single_pos_int(set)
  assert_leq(set, length(project$parameter_sets), message = "chosen set outside range")
  
  # change active set
  project$active_set <- set
  message(sprintf("active set = %s", set))
  
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
#'   case they are updated once at the end to avoid large amounts of output
#' @param store_acceptance whether to store acceptance rates for all parameters
#'   updated by Metropolis-Hastings. Proposal distributions are tuned adaptively
#'   with a target acceptance rate of 23\%
#' @param store_raw whether to store raw MCMC output in addition to summary 
#'   output. Setting to FALSE can considerably reduce output size in memory
#' @param silent whether to suppress all console output
#'
#' @export
#' @examples
#' # TODO

run_mcmc <- function(project, K = NULL, precision = 0.01, burnin = 1e3, samples = 1e3, rungs = 1, GTI_pow = 3, auto_converge = TRUE, converge_test = 1e2, solve_label_switching_on = TRUE, coupling_on = TRUE, cluster = NULL, pb_markdown = FALSE, store_acceptance = TRUE, store_raw = TRUE, silent = FALSE) {
  
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
  assert_single_logical(store_acceptance)
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
  if (data_format == "biallelic") {
    dat = mat_to_rcpp(dat)
  }
  
  # ---------- create argument lists ----------
  
  # data list
  args_data <- list(data = dat,
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
                      store_acceptance = store_acceptance,
                      silent = !is.null(cluster))
  
  # combine model parameters list with input arguments
  args_model <- c(project$parameter_sets[[s]], args_inputs)
  
  # get COI_model in numeric form
  args_model$COI_model_numeric <- match(args_model$COI_model, c("uniform", "poisson" ,"nb"))
  
  # record type of lambda (scalar, vector, list)
  lambda_type <- 1
  if (!is.list(args_model$lambda)) {
    if (length(args_model$lambda) != 1) {
      lambda_type <- 2
    }
  } else {
    lambda_type <- 3
  }
  args_model$lambda_type <- lambda_type
  
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
  all_converged <- TRUE
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
    
    # get full COI trace in coda::mcmc format
    full_COI <- mcmc(rcpp_to_mat(output_raw[[i]]$m_store))
    colnames(full_COI) <- sample_names
    
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
    if (data_format == "biallelic" && args_model$estimate_error) {
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
    
    # get whether rungs have converged
    converged <- output_raw[[i]]$rung_converged
    if (all_converged && any(!converged)) {
      all_converged <- FALSE
    }
    
    # ---------- summary results ----------
    
    # get 95% credible intervals over sampling loglikelihoods
    loglike_intervals <- t(apply(loglike_sampling, 2, quantile_95))
    rownames(loglike_intervals) <- rung_names
    class(loglike_intervals) <- "malecot_loglike_intervals"
    
    # get 95% credible intervals over COI
    COI_intervals <- t(apply(full_COI, 2, quantile_95))
    class(COI_intervals) <- "malecot_COI_intervals"
    
    # get COI matrix
    COI_matrix <- apply(full_COI, 2, function(x) tabulate(x, nbins = args_model$COI_max))/samples
    
    # get 95% credible intervals over p
    p_intervals <- list()
    for (k in 1:K[i]) {
      p_intervals[[k]] <- list()
      for (j in 1:L) {
        p_intervals[[k]][[j]] <- t(apply(full_p[[k]][[j]], 2, quantile_95))
      }
      names(p_intervals[[k]]) <- locus_names
      class(p_intervals[[k]]) <- "malecot_p_intervals"
    }
    names(p_intervals) <- deme_names
    
    # get 95% credible intervals over COI_mean
    if (is.null(full_COI_mean)) {
      fake_COI_mean <- matrix((args_model$COI_max + 1)/2, ncol = K[i])
      colnames(fake_COI_mean) <- deme_names
      COI_mean_intervals <- t(apply(fake_COI_mean, 2, quantile_95))
    } else {
      COI_mean_intervals <- t(apply(full_COI_mean, 2, quantile_95))
    }
    class(COI_mean_intervals) <- "malecot_COI_mean_intervals"
    
    # get 95% credible intervals over e1 and e2
    e_intervals <- NULL
    if (data_format == "biallelic") {
      if (args_model$estimate_error) {
        e_intervals <- rbind(quantile_95(full_e1), quantile_95(full_e2))
      } else {
        e_intervals <- rbind(quantile_95(args_model$e1), quantile_95(args_model$e2))
      }
      rownames(e_intervals) <- c("e1", "e2")
      class(e_intervals) <- "malecot_e_intervals"
    }
    
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
    p_accept <- NULL
    m_accept <- NULL
    e_accept <- NULL
    coupling_accept <- NULL
    if (store_acceptance) {
      p_accept <- rcpp_to_mat(output_raw[[i]]$p_accept)/samples
      m_accept <- output_raw[[i]]$m_accept/samples
      e_accept <- output_raw[[i]]$e_accept/samples
      coupling_accept <- output_raw[[i]]$coupling_accept/samples
    }
    
    # ---------- save arguments ----------
    
    output_args <- list(precision = precision,
                        burnin = burnin,
                        samples = samples,
                        rungs = rungs,
                        GTI_pow = GTI_pow,
                        auto_converge = auto_converge,
                        converge_test = converge_test,
                        solve_label_switching_on = solve_label_switching_on,
                        coupling_on = coupling_on,
                        pb_markdown = pb_markdown,
                        store_acceptance = store_acceptance,
                        store_raw = store_raw,
                        silent = silent)
    
    # ---------- save results ----------
    
    # add to project
    project$output$single_set[[s]]$single_K[[K[i]]] <- list()
    
    project$output$single_set[[s]]$single_K[[K[i]]]$summary <- list(qmatrix = qmatrix,
                                                                    loglike_intervals = loglike_intervals,
                                                                    COI_matrix = COI_matrix,
                                                                    COI_intervals = COI_intervals,
                                                                    p_intervals = p_intervals,
                                                                    e_intervals = e_intervals,
                                                                    COI_mean_intervals = COI_mean_intervals,
                                                                    ESS = ESS,
                                                                    GTI_path = GTI_path,
                                                                    GTI_logevidence = GTI_logevidence,
                                                                    converged = converged)
    if (store_acceptance) {
      project$output$single_set[[s]]$single_K[[K[i]]]$summary$p_accept <- p_accept
      project$output$single_set[[s]]$single_K[[K[i]]]$summary$m_accept <- m_accept
      project$output$single_set[[s]]$single_K[[K[i]]]$summary$e_accept <- e_accept
      project$output$single_set[[s]]$single_K[[K[i]]]$summary$coupling_accept <- coupling_accept
    }
    
    if (store_raw) {
      project$output$single_set[[s]]$single_K[[K[i]]]$raw <- list(loglike_burnin = loglike_burnin,
                                                                  loglike_sampling = loglike_sampling,
                                                                  COI = full_COI,
                                                                  p = full_p,
                                                                  e1 = full_e1,
                                                                  e2 = full_e2,
                                                                  COI_mean = full_COI_mean)
    }
    
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
  if (!silent) {
    tdiff <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (tdiff < 60) {
      message(sprintf("Total run-time: %s seconds", round(tdiff, 2)))
    } else {
      message(sprintf("Total run-time: %s minutes", round(tdiff/60, 2)))
    }
  }
  
  # warning if any rungs in any MCMCs did not converge
  if (!all_converged && !silent) {
    message("\n**WARNING** at least one MCMC run did not converge within specified burn-in\n")
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
    
    if (!is.null(x[[i]]$raw)) {
      
      p_raw <- x[[i]]$raw$p[best_perm_order]
      names(p_raw) <- deme_names
      proj$output$single_set[[s]]$single_K[[i]]$raw$p <- p_raw
      
      COI_mean_raw <- x[[i]]$raw$COI_mean
      if (!is.null(COI_mean_raw)) {
        COI_mean_raw <- COI_mean_raw[, best_perm_order, drop = FALSE]
        colnames(COI_mean_raw) <- deme_names
        proj$output$single_set[[s]]$single_K[[i]]$raw$COI_mean <- COI_mean_raw
      }
    }
    
    p_intervals <- x[[i]]$summary$p_intervals[best_perm_order]
    names(p_intervals) <- deme_names
    proj$output$single_set[[s]]$single_K[[i]]$summary$p_intervals <- p_intervals
    
    COI_mean_quantiles <- x[[i]]$summary$COI_mean_quantiles
    if (!is.null(COI_mean_quantiles)) {
      COI_mean_quantiles <- COI_mean_quantiles[best_perm_order, , drop = FALSE]
      rownames(COI_mean_quantiles) <- deme_names
      class(COI_mean_quantiles) <- "malecot_COI_mean_quantiles"
      proj$output$single_set[[s]]$single_K[[i]]$summary$COI_mean_quantiles <- COI_mean_quantiles
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
#' @title Match grouping against q-matrix
#'
#' @description Compares qmatrix output for a chosen value of K against a 
#'   \code{target_group} vector. Returns the order of \code{target_group} 
#'   groups, such that there is the best possible alignment against the qmatrix.
#'   For example, if the vector returned is \code{c(2,3,1)} then the second
#'   group in the target vector should be matched against the first group in the
#'   qmatrix, followed by the third group in the target vector against the
#'   second group in the qmatrix, followed by the first group in the target
#'   vector against the third group in the qmatrix.
#'
#' @param project a MALCOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K compare against qmatrix output for this value of K
#' @param target_group the target group to be aligned against the qmatrix
#'
#' @export
#' @examples
#' # TODO

get_group_order <- function(project, K, target_group) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  assert_single_pos_int(K, zero_allowed = FALSE)
  assert_vector(target_group)
  assert_pos_int(target_group)
  
  # create target_qmatrix from target_group
  n <- length(target_group)
  target_qmatrix <- matrix(0, n, max(target_group))
  target_qmatrix[cbind(1:n, target_group)] <- 1
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # check output exists for chosen K
  qmatrix <- project$output$single_set[[s]]$single_K[[K]]$summary$qmatrix
  if (is.null(qmatrix)) {
    stop(sprintf("no qmatrix output for K = %s of active set", K))
  }
  
  # check same number of rows
  if (nrow(qmatrix) != nrow(target_qmatrix)) {
    stop("target_group must have same number of elements as qmatrix output")
  }
  n <- nrow(qmatrix)
  
  # expand qmatrix and/or target_qmatrix to get same number of cols
  if (ncol(target_qmatrix) > ncol(qmatrix)) {
    qmatrix <- cbind( qmatrix, matrix(0, n, ncol(target_qmatrix) - ncol(qmatrix)) )
  }
  if (ncol(qmatrix) > ncol(target_qmatrix)) {
    target_qmatrix <- cbind( target_qmatrix, matrix(0, n, ncol(qmatrix) - ncol(target_qmatrix)) )
  }
  
  # calculate cost matrix
  cost_mat <- matrix(0,K,K)
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      cost_mat[k1,k2] <- sum(qmatrix[,k1]*(1-target_qmatrix[,k2]))
    }
  }
  
  # get lowest cost permutation
  best_perm <- call_hungarian(cost_mat)$best_matching + 1
  best_perm_order <- order(best_perm)
  
  return(best_perm)
}

#------------------------------------------------
#' @title Get ESS
#'
#' @description Returns effective sample size (ESS) of chosen model run.
#'
#' @param project a MALCOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K get ESS for this value of K
#'
#' @export
#' @examples
#' # TODO

get_ESS <- function(project, K = NULL) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$ESS)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no ESS output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  ESS <- project$output$single_set[[s]]$single_K[[K]]$summary$ESS
  if (is.null(ESS)) {
    stop(sprintf("no ESS output for K = %s of active set", K))
  }
  
  return(ESS)
}
