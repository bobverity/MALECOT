# MALECOT
### Version 0.1.0
[![Travis-CI Build Status](https://travis-ci.org/bobverity/MALECOT.svg?branch=master)](https://travis-ci.org/bobverity/MALECOT)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/bobverity/MALECOT?branch=master&svg=true)](https://ci.appveyor.com/project/bobverity/MALECOT)
[![Coverage Status](https://img.shields.io/codecov/c/github/bobverity/MALECOT/master.svg)](https://codecov.io/github/bobverity/MALECOT?branch=master)
[![Documentation](https://github.com/bobverity/MALECOT/blob/master/R_ignore/images/documentation-click%20here!-blue.png)](https://bobverity.github.io/MALECOT/)

--------------------------------------------------------------------------------------------------------------------------------

MALECOT is an R package for inferring population structure from genetic data, specifically in situations where each sample might represent multiple genotypes. For example, in malaria it is common for a single human host to be infected by multiple parasite "strains", and so what we see when we sequence extracted parasites is a convolution of multiple parasite genotypes. MALECOT simultaneously estimates the number of parasite genotypes (i.e the complexity of infection, or COI) of each sample along with the population structure. This is important in obtaining accurate estimates of both the underlying population structure and the true COI, as failing to account for these properties can lead to biased inference.

MALECOT uses Bayesian Markov chain Monte Carlo (MCMC) to estimate the full posterior distribution of all unknown parameters. Core functions are written in C++ through the Rcpp package for increased speed. It also contains methods for estimating the number of clusters (K), and can handle biallelic and more-than-biallelic data (i.e. SNPs, microsats, haplotypes). Finally, the package contains a number of diagnostic and plotting functions for checking MCMC performance, and visualising final results.

*WARNING: The current program is live, but still very much in development!*
