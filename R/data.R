
#------------------------------------------------
#' @title Simulate genetic data
#'
#' @description Generates simulated data in the format accepted by MALECOT.
#'   Simulated data can be either bi-allelic or multi-allelic. In addition to
#'   returning the simulated data set, this function returns intermediate
#'   outputs used in data generation, including the true complexity of infection
#'   and true subpopulation membership of all samples. These additional outputs
#'   will not generally be present in real-world data sets, but can be useful in
#'   assessing the accuracy of the MCMC inference step by "truthing" estimated
#'   values (this is one of the major advantages of using simulated data).
#'
#' @details The \code{sim_data()} function returns a list containing multiple
#'   elements The main output is the element \code{data}, which is the final
#'   simulated data set, along with the number of samples (\code{n}) and the
#'   number of loci (\code{L}). If the parameters \code{e1} or \code{e2} are
#'   non-zero, or if \code{prop_missing} is greater than zero, then an
#'   uncorrupted version of the data is also returned in the element
#'   \code{data_uncorrupted}. This object represents what the true underlying
#'   data looked like before errors and missing data were introduced. Other
#'   "truthing" elements include the true subpopulation grouping, true COI, and
#'   true allele frequencies that were simulated to produce the final data. The
#'   accuracy of MCMC-based inference can be assessed by comparing posterior
#'   estimates against these "true" values used in simulation. Finally, the
#'   original call to the function is also stored in the element \code{call} for
#'   the sake of producing a complete record of all simulation parameters used.
#'
#'   The prior on COI can be uniform, Poisson, or negative binomial. In the
#'   uniform case there is an equal chance of any given sample having a COI
#'   between 1 and \code{COI_max} (inclusive). In the Poisson and negative
#'   binomial cases it is important to note that the distribution is over
#'   (COI-1), rather than over COI. This is because both Poisson and negative
#'   binomial distributions allow for 0 values, which cannot be the case here
#'   because observed samples must contain at least 1 genotype.
#'
#'   The full probability mass distribution for the Poisson case with
#'   \code{COI_mean}\eqn{=\mu} and \code{COI_max}\eqn{=M} can be written \deqn{
#'   Pr(COI=n) = z (\mu-1)^(n-1) exp(-(\mu-1)) / (n-1)! } where z is a
#'   normalising constant that ensures the distribution sums to unity: \deqn{1/z
#'   = \sum_{i=1}^M (\mu-1)^(i-1) exp(-(\mu-1)) / (i-1)! }
#'
#'   The mean of this distribution will generally be very close to \eqn{\mu},
#'   and the variance will be close to \eqn{\mu-1} (strictly it will approach
#'   these values as \eqn{M} tends to infinity).
#'
#'   The full probability mass distribution for the negative binomial case with
#'   \code{COI_mean}\eqn{=\mu}, \code{COI_dispersion}\eqn{=v/\mu} and
#'   \code{COI_max}\eqn{=M} can be written \deqn{ Pr(COI=n) = z \Gamma(n-1+N)/(
#'   \Gamma(N)(n-1)! ) p^N (1-p)^(n-1) } where \eqn{N=(\mu-1)^2/(v-\mu+1)},
#'   \eqn{p=(\mu-1)/v}, and z is a normalising constant that ensures the
#'   distribution sums to unity: \deqn{1/z = \sum_{i=1}^M \Gamma(i-1+N)/(
#'   \Gamma(N)(i-1)! ) p^N (1-p)^(i-1) }
#'
#'   The mean of this distribution will generally be very close to \eqn{\mu},
#'   and the variance will be close to \eqn{v} (strictly it will approach these
#'   values as \eqn{M} tends to infinity).
#'
#' @references Chang, H.H., Worby, C.J., Yeka, A., Nankabirwa, J., Kamya, M.R.,
#' Staedke, S.G., Dorsey, G., Murphy, M., Neafsey, D.E., Jeffreys, A.E. and
#' Hubbart, C., 2017. THE REAL McCOIL: A method for the concurrent estimation of
#' the complexity of infection and SNP allele frequency for malaria parasites.
#' PLoS computational biology, 13(1), p.e1005348.
#'
#' @param n number of samples
#' @param L number of loci
#' @param K number of subpopulations. Samples are drawn with equal probability
#'   from each of these subpopulations
#' @param data_format character string indicating whether data should be
#'   \code{"biallelic"} or \code{"multiallelic"}. In the bi-allelic case data
#'   are formatted similar to The Real McCOIL (Chang et al. 2017). In the
#'   multi-allelic case data are in long format, allowing a different number of
#'   haplotypes per sample
#' @param alleles vector of length \code{L} indicating the number of possible
#'   alleles at each locus, or if a single value is used then the same value
#'   applies for all loci. In the bi-allelic case this parameter is ignored as
#'   the number of alleles is fixed at 2
#' @param lambda shape parameters governing the prior on allele frequencies,
#'   which is a Beta distribution for the bi-allelic case or a Dirichlet
#'   distribution for the multi-allelic case. In the simplest case this can be a
#'   single scalar value, in which case the same shape paremeter will be used
#'   for all loci and all alleles (i.e. a symmetric Beta or Dirichlet prior). In
#'   the more general case this can be a list of length \code{L}, containing
#'   vectors of length equal to the number of alleles at each locus, allowing
#'   different shape parameters to be defined for each allele individually
#' @param COI_model the type of prior distribution on COI in each subpopulation.
#'   Can be \code{uniform}, \code{poisson} or \code{nb} (negative binomial).
#'   Note that the poisson and negative binomial distributions use a
#'   non-standard parametrization to account for the fact that COI=0 is not
#'   possible (see Details). If \code{COI_manual} is used then this prior is
#'   ignored.
#' @param COI_manual manually define COI of all samples
#' @param COI_max the maximum possible COI in any individual, irrespective of
#'   the type of prior on COI
#' @param COI_mean the mean COI in a given subpopulation (applies to
#'   \code{poisson} and \code{nb} priors only). If vector valued then should be
#'   of length \code{K} and defines the mean separately for each subpopulation.
#'   If scalar valued then the same value applies to all subpopulations. Note
#'   that this mean value will only be correct when \code{COI_max} is
#'   sufficiently large (see Details)
#' @param COI_dispersion the index of dispersion (i.e. the variance over the
#'   mean) of the prior on COI. Applies to \code{nb} prior only - in the Poisson
#'   case the dispersion is fully defined by the \code{COI_mean}. If vector
#'   valued then should be of length \code{K} and defines the dispersion
#'   separately for each subpopulation. If scalar valued then the same value
#'   applies to all subpopulations. Note that this dispersion will only be
#'   correct when \code{COI_max} is sufficiently large (see Details)
#' @param e1 first error term: the probabilty of a true homozygous locus being
#'   mis-called as a heterozygote. Errors are only allowed under the
#'   \code{biallelic} model
#' @param e2 second error term: the probabilty of a true heterozygous locus
#'   being mis-called as a homozygote. Errors are only allowed under the
#'   \code{biallelic} model
#' @param prop_missing the proportion of missing data. Missing data is randomly
#'   distributed between samples and loci, and is coded as \code{-1}
#'
#' @export
#' @examples
#' # TODO

sim_data <- function(n = 100, L = 24, K = 3, data_format = "biallelic", alleles = 2, lambda = 2, COI_model = "poisson", COI_manual = NULL, COI_max = 20, COI_mean = 3, COI_dispersion = 1, e1 = 0, e2 = 0, prop_missing = 0) {

  # check inputs and force certain formats
  assert_pos_int(n, zero_allowed = FALSE)
  assert_pos_int(L, zero_allowed = FALSE)
  assert_pos_int(K, zero_allowed = FALSE)
  assert_in(data_format, c("biallelic", "multiallelic"))
  assert_pos_int(alleles, zero_allowed = FALSE)
  assert_gr(alleles, 1)
  assert_in(length(alleles), 1, L)
  assert_bounded(lambda, left = 0, right = 100)
  assert_in(COI_model, c("uniform", "poisson", "nb"))
  if (!is.null(COI_manual)) {
    assert_that(length(COI_manual) == n)
    assert_gr(COI_manual, 1)
    COI_mean <- NA
    COI_dispersion <- NA
  } else {
    assert_pos_int(COI_max, zero_allowed = FALSE)
    assert_bounded(COI_mean, left = 1, right = COI_max)
    assert_in(length(COI_mean), 1, K)
    assert_gr(COI_dispersion, 1 - 1/COI_mean)
    assert_in(length(COI_dispersion), 1, K)
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

  # force COI mean and dispersion to vector
  if (length(COI_mean)==1) {
    COI_mean <- rep(COI_mean, K)
  }
  if (length(COI_dispersion)==1) {
    COI_dispersion <- rep(COI_dispersion, K)
  }

  # check that lambda compatible with number of alleles
  if (!all(mapply(length, lambda)==alleles)) {
    stop("lambda not compatible with number of alleles")
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

  # name outputs
  names(true_group) <- ind_names
  names(true_m) <- ind_names
  names(true_p) <- locus_names
  names(true_p) <- locus_names

  # simulate bi-allelic data
  if (data_format=="biallelic") {

    # simulate raw data
    data <- NULL
    for (k in 1:K) {
      if (any(true_group==k)) {
        true_m_k <- true_m[true_group==k]
        data_raw_k <- t(mapply(rbinom, n=L, size=true_m_k, MoreArgs=list(prob=true_p[[k]])))
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
