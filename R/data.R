
#------------------------------------------------
#' @title Simulate genetic data
#'
#' @description Simulate genetic data from the same model used in the MALECOT inference step.
#'
#' @details TODO
#'
#' @param n the number of samples
#' @param L the number of loci per sample
#' @param K the number of subpopulations
#' @param data_format whether to produce data in "biallelic" or "multiallelic"
#'   format. Note that if biallelic format is chosen then \code{alleles} is
#'   always set to 2
#' @param pop_col_on TODO
#' @param alleles the number of alleles at each locus. Can be a vector of length
#'   \code{L} specifying the number of alleles at each locus, or a single scalar
#'   value specifying the number of alleles at all loci
#' @param lambda the shape parameter(s) of the prior on allele frequencies. This
#'   prior is Beta in the bi-allelic case, and Dirichlet in the multi-allelic
#'   case. \code{lambda} can be a list of length \code{L} containing vectors of
#'   length equal to the number of alleles at that locus, or it can be a single
#'   scalar value. If a list then \code{lambda} specifies all shape parameters
#'   of the prior at each locus separately, if a scalar value then the same
#'   shape parameter is used over all loci and all alleles
#' @param COI_model TODO
#' @param COI_max TODO
#' @param COI_manual option to override the MCMC and set the COI of one or more
#'   samples manually, in which case they are not updated. Vector of length
#'   \code{n} specifing the integer valued COI of each sample, with -1
#'   indicating that all samples should be estimated
#' @param COI_mean TODO
#' @param COI_dispersion must be > 1. Only used under the negative binomial
#'   model
#' @param e1 TODO
#' @param e2 TODO
#' @param prop_missing TODO
#'
#' @export
#' @examples
#' # TODO

sim_data <- function(n = 100, L = 24, K = 3, data_format = "biallelic", pop_col_on = TRUE, alleles = 2, lambda = 1, COI_model = "poisson", COI_max = 20, COI_manual = rep(-1,n), COI_mean = 3, COI_dispersion = 2, e1 = 0, e2 = 0, prop_missing = 0) {
  
  ##### CHECK INPUTS #####
  
  assert_single_pos_int(n, zero_allowed = FALSE)
  assert_single_pos_int(L, zero_allowed = FALSE)
  assert_single_pos_int(K, zero_allowed = FALSE)
  assert_in(data_format, c("biallelic", "multiallelic"))
  if (data_format != "biallelic") {
    assert_vector(alleles)
    assert_pos_int(alleles, zero_allowed = FALSE)
    assert_gr(alleles, 1)
    assert_in(length(alleles), c(1, L))
  }
  assert_pos(unlist(lambda), zero_allowed = FALSE)
  assert_single_string(COI_model)
  assert_in(COI_model, c("uniform", "poisson", "nb"))
  assert_single_pos_int(COI_max, zero_allowed = FALSE)
  assert_vector(COI_manual)
  assert_length(COI_manual, n)
  assert_int(COI_manual)
  assert_bounded(COI_manual[COI_manual != -1], left = 1, right = COI_max, inclusive_left = TRUE, inclusive_right = TRUE)
  if (COI_model != "uniform") {
    assert_vector(COI_mean)
    assert_in(length(COI_mean), c(1, K))
    assert_bounded(COI_mean, left = 1, right = COI_max, inclusive_left = TRUE, inclusive_right = TRUE)
    if (COI_model == "nb") {
      assert_vector(COI_dispersion)
      assert_in(length(COI_dispersion), c(1, K))
      assert_gr(COI_dispersion, 1)
    }
  }
  assert_single_numeric(e1)
  assert_bounded(e1, left = 0.0, right = 1.0, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_numeric(e2)
  assert_bounded(e2, left = 0.0, right = 1.0, inclusive_left = TRUE, inclusive_right = TRUE)
  assert_single_numeric(prop_missing)
  assert_bounded(prop_missing, left = 0.0, right = 1.0, inclusive_left = TRUE, inclusive_right = TRUE)
  
  # force alleles = 2 if biallelic format selected
  if (data_format == "biallelic") {
    if (!all(alleles == 2)) {
      message(sprintf("Note: biallelic format specified, therefore overriding argument 'alleles = %s'", alleles))
    }
    alleles <- 2
  }
  
  # force alleles to vector over loci
  if (length(alleles) == 1) {
    alleles <- rep(alleles, L)
  }
  
  # force lambda to list
  if (!is.list(lambda)) {
    lambda <- lapply(alleles, function(x) {rep(lambda,x)})
  }
  
  # check that lambda compatible with number of alleles
  if (!all(mapply(length, lambda)==alleles)) {
    stop("lambda not compatible with number of alleles")
  }
  
  # define COI mean and dispersion under uniform and Poisson priors
  if (COI_model =="uniform") {
    COI_mean <- (COI_max+1)/2
    COI_dispersion <- 1/6*(COI_max+1)*(2*COI_max+1)/COI_mean - COI_mean
  } else if (COI_model == "poisson") {
    COI_dispersion <- 1 - 1/COI_mean
  }
  
  # force COI mean and dispersion to vector
  if (length(COI_mean) == 1) {
    COI_mean <- rep(COI_mean, K)
  }
  if (length(COI_dispersion) == 1) {
    COI_dispersion <- rep(COI_dispersion, K)
  }
  
  
  ##### SIMULATE DATA #####
  
  # create names
  samp_names <- paste0("samp", zero_pad_simple(1:n, nchar(n)))
  locus_names <- paste0("locus", 1:L)
  pop_names <- paste0("pop", 1:K)
  
  # generate true grouping and allele frequencies
  true_group <- sort(sample(K, n, replace = TRUE))
  true_p <- list()
  for (l in 1:L) {
    allele_names <- paste0("allele", 1:alleles[l])
    if (alleles[l] == 1) {
      true_p[[l]] <- matrix(1, nrow = K)
    } else {
      true_p[[l]] <- t(mapply(rdirichlet, replicate(K, lambda[[l]], simplify = FALSE)))
    }
    rownames(true_p[[l]]) <- pop_names
    colnames(true_p[[l]]) <- allele_names
  }
  names(true_group) <- samp_names
  names(true_p) <- locus_names
  
	# generate true COIs
  switch(COI_model,
         "uniform" = {
           true_m <- sample(1:COI_max, n, replace = TRUE)
         },
         "poisson" = {
           mu <- COI_mean[true_group]
           true_m <- rpois(n, lambda=mu-1) + 1
         },
         "nb" = {
           mu <- COI_mean[true_group]
           v <- mu * COI_dispersion[true_group]
           true_m <- rnbinom(n, size=(mu-1)^2/(v-mu+1), prob=(mu-1)/v) + 1
         })
  true_m[COI_manual != -1] <- COI_manual[COI_manual != -1]
  names(true_m) <- samp_names
  
  # truncate COIs at COI_max
  true_m[true_m > COI_max] <- COI_max
  
  # simulate bi-allelic or multi-allelic data
  if (data_format == "biallelic") {
    dat_raw <- sim_data_biallelic(n = n,
                                  L = L,
                                  K = K,
                                  true_group = true_group,
                                  true_p = true_p,
                                  true_m = true_m,
                                  locus_names = locus_names,
                                  samp_names = samp_names,
                                  e1 = e1,
                                  e2 = e2,
                                  prop_missing = prop_missing,
                                  pop_col_on)
  } else {
    dat_raw <- sim_data_multiallelic(n = n,
                                     L = L,
                                     K = K,
                                     true_group = true_group,
                                     true_p = true_p,
                                     true_m = true_m,
                                     locus_names = locus_names,
                                     samp_names = samp_names,
                                     e1 = e1,
                                     e2 = e2,
                                     prop_missing = prop_missing,
                                     pop_col_on)
  }
  
  # return simulated data and true parameter values
  output <- list()
  output$data <- dat_raw$df
  output$n <- n
  output$L <- L
  output$data_uncorrupted <- dat_raw$df_uncorrupted
  output$true_group <- true_group
  output$true_m <- true_m
  output$true_p <- true_p
  output$call <- match.call()
  
  return(output)
}

#------------------------------------------------
# simulate bi-allelic data
#' @noRd
sim_data_biallelic <- function(n, L, K, true_group, true_p = true_p, true_m, locus_names, samp_names, e1, e2, prop_missing, pop_col_on) {
  
  # simulate raw data
  dat <- NULL
  for (k in 1:K) {
    if (any(true_group == k)) {
      
      # draw raw numbers of REF allele for individuals in this deme
      true_m_k <- true_m[true_group == k]
      true_p_k <- mapply(function(x) {x[k,1]}, true_p)
      dat_raw_k <- t(mapply(rbinom, n = L, size = true_m_k, MoreArgs = list(prob = true_p_k)))
      
      # convert to matrix of {0.0, 0.5, 1.0}
      dat_k <- matrix(0.5, length(true_m_k), L)
      dat_k[dat_raw_k == matrix(true_m_k, length(true_m_k), L)] <- 1
      dat_k[dat_raw_k == 0] <- 0
      
      # append data
      dat <- rbind(dat, dat_k)
    }
  }
  colnames(dat) <- locus_names
  
  # add errors and missing data
  dat_uncorrupted <- NULL
  if (e1>0 || e2>0 || prop_missing>0) {
    dat_uncorrupted <- dat
    
    # error1 - homo missclassified as het
    if (e1 > 0) {
      homo1 <- dat[dat_uncorrupted == 1]
      homo1[runif(length(homo1)) < e1] <- 0.5
      dat[dat_uncorrupted == 1] <- homo1
      
      homo2 <- dat[dat_uncorrupted == 0]
      homo2[runif(length(homo2)) < e1] <- 0.5
      dat[dat_uncorrupted == 0] <- homo2
    }
    
    # error2 - het missclassified as homo
    if (e2 > 0) {
      het <- dat[dat_uncorrupted == 0.5]
      rand1 <- (runif(length(het)) < e2)
      if (any(rand1)) {
        het[rand1] <- sample(c(0,1), sum(rand1), replace = TRUE)
      }
      dat[dat_uncorrupted == 0.5] <- het
    }
    
    # missing data
    if (prop_missing > 0) {
      prop_missing_round <- round(prop_missing*n*L)
      dat[sample.int(n*L, prop_missing_round )] <- -9
    }
  }
  
  # convert dat and dat_uncorrupted to dataframe
  df <- data.frame(sample_ID = samp_names, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  if (pop_col_on) {
    df$pop <- true_group
  }
  df_uncorrupted <- NULL
  if (!is.null(dat_uncorrupted)) {
    df_uncorrupted <- cbind(df, dat_uncorrupted)
  }
  df <- cbind(df, dat)
  
  # return list
  return(list(df = df, df_uncorrupted = df_uncorrupted))
}

#------------------------------------------------
# simulate multi-allelic data
#' @noRd
sim_data_multiallelic <- function(n, L, K, true_group, true_p = true_p, true_m, locus_names, samp_names, e1, e2, prop_missing, pop_col_on) {
  
  # simulate raw data
  df <- NULL
  for (l in 1:L) {
    true_p_group <- lapply(true_group, function(i) {true_p[[l]][i,]})
    haplotypes <- mapply(function(x,y) {
      sort(unique(sample.int(length(x), y, replace = TRUE, prob = x)))
    }, true_p_group, y = true_m, SIMPLIFY = FALSE)
    df_l <- data.frame(sample_ID = rep(samp_names, times = sapply(haplotypes,length)), locus = l, haplotype = unlist(haplotypes), stringsAsFactors = FALSE)
    df <- rbind(df, df_l)
  }
  df <- df[order(df$sample),]
  row.names(df) <- NULL
  
  # add missing data
  df_uncorrupted <- NULL
  if (prop_missing > 0) {
    df_uncorrupted <- df
    
    prop_missing_round <- round(prop_missing*n*L)
    missing_index <- expand.grid(unique(df$sample_ID), 1:L, -9)[sample.int(n*L, prop_missing_round, replace = FALSE),]
    names(missing_index) <- c("sample_ID", "locus", "haplotype")
    df <- subset(df, !( paste(df$sample_ID, df$locus, sep=".") %in% paste(missing_index$sample_ID, missing_index$locus, sep=".") ))
    df <- rbind(df, missing_index)
    df <- df[order(df$locus),]
    df <- df[order(df$sample),]
    row.names(df) <- NULL
  }
  
  # return list
  return(list(df = df, df_uncorrupted = df_uncorrupted))
}

#------------------------------------------------
#' @title Simulate genetic data subject to constraints
#'
#' @description TODO - text
#'
#' @details TODO
#'
#' @param ... TODO
#' @param data_format TODO
#' @param no_invariant_loci TODO
#' @param no_missing_samples TODO
#' @param no_missing_loci TODO
#' @param max_attempts TODO
#'
#' @export
#' @examples
#' # TODO

sim_data_safe <- function(..., data_format = "biallelic", no_invariant_loci = TRUE, no_missing_samples = TRUE, no_missing_loci = TRUE, max_attempts = 1e3) {

  # attempt to simulate satisfactory data a finite number of times
  for (i in 1:max_attempts) {

    # simulate data
    sim1 <- sim_data(..., data_format = data_format)
    data <- sim1$data
    n <- sim1$n
    L <- sim1$L

    # bi-allelic data
    if (data_format=="biallelic") {

      # check no invariant loci
      if (no_invariant_loci & any(colSums(data==1 | data==0)==nrow(data)) ) {
        next
      }

      # check no missing samples
      if ( no_missing_samples & any(colSums(data==-1)==nrow(data)) ) {
        next
      }

      # check no missing loci
      if ( no_missing_loci & any(rowSums(data==-1)==ncol(data)) ) {
        next
      }
    }

    # multi-allelic data
    if (data_format=="multiallelic") {

      # check no invariant loci
      n_haplotypes <- mapply(function(i) {
        s <- data$haplotype[data$locus==i]
        length(unique(s[s>0]))
        }, 1:L)
      if (no_invariant_loci & any(n_haplotypes==1)) {
        next
      }

      # check no missing samples
      n_nonmissing_samples <- mapply(function(i){
        sum(data$sample==i & data$haplotype>0)
        }, 1:n)
      if (no_missing_samples & any(n_nonmissing_samples==0)) {
        next
      }

      # check no missing loci
      n_nonmissing_loci <- mapply(function(i){
        sum(data$locus==i & data$haplotype>0)
        }, 1:L)
      if ( no_missing_loci & any(n_nonmissing_loci==0) ) {
        next
      }
    }

    # if made it to here then data passed all checks
    return(sim1)
  }

  stop(paste("Unable to produce data set satisfying constraints within", max_attempts, "random draws"))
}
