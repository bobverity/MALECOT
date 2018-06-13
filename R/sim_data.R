
#------------------------------------------------
#' @title Simulate genetic data
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param n TODO
#' @param L TODO
#' @param K TODO
#' @param data_format TODO
#' @param alleles TODO
#' @param lambda TODO
#' @param COI_model TODO
#' @param COI_manual TODO
#' @param COI_max TODO
#' @param COI_mean TODO
#' @param COI_dispersion TODO
#' @param e1 TODO
#' @param e2 TODO
#' @param prop_missing TODO
#'
#' @export
#' @examples
#' # TODO

sim_data <- function(n = 100, L = 24, K = 3, data_format = "biallelic", alleles = 2, lambda = 1, COI_model = "poisson", COI_manual = NULL, COI_max = 20, COI_mean = 3, COI_dispersion = 1, e1 = 0, e2 = 0, prop_missing = 0) {
  
  # check inputs and force certain formats
  assert_pos_int(n, zero_allowed = FALSE)
  assert_pos_int(L, zero_allowed = FALSE)
  assert_pos_int(K, zero_allowed = FALSE)
  assert_in(data_format, c("biallelic", "multiallelic"))
  assert_pos_int(alleles, zero_allowed = FALSE)
  assert_gr(alleles, 1)
  assert_in(length(alleles), c(1, L))
  assert_bounded(unlist(lambda), left = 0, right = 100)
  assert_in(COI_model, c("uniform", "poisson", "nb"))
  if (!is.null(COI_manual)) {
    assert_that(length(COI_manual) == n)
    assert_gr(COI_manual, 1)
    COI_mean <- NA
    COI_dispersion <- NA
  } else {
    assert_pos_int(COI_max, zero_allowed = FALSE)
    assert_bounded(COI_mean, left = 1, right = COI_max)
    assert_in(length(COI_mean), c(1, K))
    assert_gr(COI_dispersion, 1 - 1/COI_mean)
    assert_in(length(COI_dispersion), c(1, K))
  }
  assert_bounded(e1)
  assert_bounded(e2)
  assert_bounded(prop_missing)
  
  # fix alleles if biallelic format
  if (data_format=="biallelic") {
    alleles <- 2
  }
  
  # define COI mean and dispersion under uniform and Poisson priors
  if (COI_model =="uniform") {
    COI_mean <- (COI_max+1)/2
    COI_dispersion <- 1/6*(COI_max+1)*(2*COI_max+1)/COI_mean - COI_mean
  } else if (COI_model =="poisson") {
    COI_dispersion <- 1 - 1/COI_mean
  }
  
  # force alleles to vector over loci
  if (length(alleles)==1) {
    alleles <- rep(alleles, L)
  }
  
  # force lambda to list
  if (is.list(lambda)) {
    if (length(lambda)!=L) {
      stop("lamba must be either a scalar value or a list of length L")
    }
  } else {
    if (length(lambda)!=1) {
      stop("lamba must be either a scalar value or a list of length L")
    }
    lambda <- lapply(alleles, function(x){rep(lambda,x)})
  }
  
  # check that lambda compatible with number of alleles
  if (!all(mapply(length, lambda)==alleles)) {
    stop("lambda not compatible with number of alleles")
  }
  
  # force COI mean and dispersion to vector
  if (length(COI_mean)==1) {
    COI_mean <- rep(COI_mean, K)
  }
  if (length(COI_dispersion)==1) {
    COI_dispersion <- rep(COI_dispersion, K)
  }

  # create names
  ind_names <- paste0("ind", 1:n)
  locus_names <- paste0("locus", 1:L)
  deme_names <- paste0("deme", 1:K)
  
  # generate true grouping and allele frequencies
  true_group <- sort(sample(K, n, replace = TRUE))
  true_p <- list()
  for (l in 1:L) {
    if (alleles[l]==1) {
      true_p[[l]] <- matrix(1,nrow=K)
    } else {
      true_p[[l]] <- t(mapply(rdirichlet, replicate(K,lambda[[l]],simplify=FALSE)))
    }
    rownames(true_p[[l]]) <- deme_names
    colnames(true_p[[l]]) <- paste0("allele", 1:alleles[l])
  }
  
	# generate true COIs
  if (!is.null(COI_manual)) {
    true_m <- COI_manual
  } else {
    if (COI_model == "uniform") {
      true_m <- COI_manual
    } else if (COI_model == "poisson") {
      mu <- COI_mean[true_group]
      true_m <- rpois(n, lambda=mu-1) + 1
    } else if (COI_model == "nb") {
      mu <- COI_mean[true_group]
      v <- mu * COI_dispersion[true_group]
      true_m <- rnbinom(n, size=(mu-1)^2/(v-mu+1), prob=(mu-1)/v) + 1
    }
  }
  
  # truncate COIs at COI_max
  true_m[true_m>COI_max] <- COI_max
  
  # initialise loglikelihood
  loglike <- 0
  
  # name outputs
  names(true_group) <- ind_names
  names(true_m) <- ind_names
  names(true_p) <- locus_names
  
  # simulate bi-allelic data
  if (data_format=="biallelic") {

    # simulate raw data
    data <- NULL
    for (k in 1:K) {
      if (any(true_group==k)) {
        true_m_k <- true_m[true_group==k]
        true_p_k <- mapply(function(x){x[k,1]}, true_p)
        data_raw_k <- t(mapply(rbinom, n=L, size=true_m_k, MoreArgs=list(prob=true_p_k)))
        data_k <- matrix(0.5, length(true_m_k), L)
        data_k[data_raw_k==matrix(true_m_k, length(true_m_k), L)] <- 1
        data_k[data_raw_k==0] <- 0
        data <- rbind(data, data_k)
      }
    }
    colnames(data) <- locus_names
    rownames(data) <- ind_names

    # add errors and missing data
    data_uncorrupted <- NULL
    if (e1>0 || e2>0 || prop_missing>0) {
      data_uncorrupted <- data

      # error1 - homo missclassified as het
      if (e1>0) {
        homo1 <- data[data_uncorrupted==1]
        homo1[runif(length(homo1))<e1] <- 0.5
        data[data_uncorrupted==1] <- homo1

        homo2 <- data[data_uncorrupted==0]
        homo2[runif(length(homo2))<e1] <- 0.5
        data[data_uncorrupted==0] <- homo2
      }

      # error2 - het missclassified as homo
      if (e2>0) {
        het <- data[data_uncorrupted==0.5]
        rand1 <- (runif(length(het))<e2)
        if (any(rand1)) {
          het[rand1] <- sample(c(0,1), sum(rand1), replace = TRUE)
        }
        data[data_uncorrupted==0.5] <- het
      }

      # missing data
      if (prop_missing>0) {
        prop_missing_round <- round(prop_missing*n*L)
        data[sample.int(n*L, prop_missing_round )] <- -1
      }
    }
    
    # calculate loglikelihood
    loglike <- 0
    for (i in 1:n) {
      this_group <- true_group[i]
      this_m <- true_m[i]
      for (j in 1:L) {
        this_p <- true_p[[j]][this_group,1]
        if (data[i,j] == 1) {
          delta <- this_m*log(this_p)
        } else if (data[i,j] == 0) {
          delta <- this_m*log(1-this_p)
        } else if (data[i,j] == 0.5) {
          delta <- log(1 - this_p^this_m - (1-this_p)^this_m)
        }
        loglike <- loglike + delta
      }
    }
    
  } # end biallelic simulation

  # simulate multi-allelic data
  if (data_format=="multiallelic") {

    # simulate raw data
    data <- NULL
    for (l in 1:L) {
      true_p_group <- lapply(true_group, function(i){true_p[[l]][i,]})
      haplotypes <- mapply(function(x,y) {
        sort(unique(sample.int(length(x), y, replace=TRUE, prob=x)))
        }, true_p_group, y=true_m)
      data_l <- data.frame(sample=rep(1:n, times=sapply(haplotypes,length)), locus=l, haplotype=unlist(haplotypes))
      data <- rbind(data, data_l)
    }
    data <- data[order(data$sample),]
    row.names(data) <- NULL

    # add missing data
    data_uncorrupted <- NULL
    if (prop_missing>0) {
      data_uncorrupted <- data

      prop_missing_round <- round(prop_missing*n*L)
      missing_index <- expand.grid(1:n,1:L,-1)[sample.int(n*L, prop_missing_round, replace = FALSE),]
      names(missing_index) <- c("sample", "locus", "haplotype")
      data <- subset(data, !( paste(data$sample, data$locus, sep=".") %in% paste(missing_index$sample, missing_index$locus, sep=".") ))
      data <- rbind(data, missing_index)
      data <- data[order(data$locus),]
      data <- data[order(data$sample),]
      row.names(data) <- NULL
    }
  } # end multiallelic simulation
  
  # return simulated data and true parameter values
  output <- list()
  output$data <- data
  output$n <- n
  output$L <- L
  output$data_uncorrupted <- data_uncorrupted
  output$true_group <- true_group
  output$true_m <- true_m
  output$true_p <- true_p
  output$loglike <- loglike
  output$call <- match.call()
  
  return(output)
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
