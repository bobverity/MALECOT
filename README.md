# MALECOT
### Version 0.1.0
[![Travis-CI Build Status](https://travis-ci.org/bobverity/MALECOT.svg?branch=master)](https://travis-ci.org/bobverity/MALECOT)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/bobverity/MALECOT?branch=master&svg=true)](https://ci.appveyor.com/project/bobverity/MALECOT)
[![Coverage Status](https://img.shields.io/codecov/c/github/bobverity/MALECOT/master.svg)](https://codecov.io/github/bobverity/MALECOT?branch=master)
[![Documentation](https://github.com/bobverity/MALECOT/blob/master/R_ignore/images/documentation-click%20here!-blue.png)](https://bobverity.github.io/MALECOT/)

--------------------------------------------------------------------------------------------------------------------------------

MALECOT is an R package for inferring population structure from genetic data in situations where each sample might represent multiple genotypes. It does this using Bayesian Markov chain Monte Carlo, alternately estimating each of three unknown parameters:
1. Allele frequencies in all subpopulations
2. Allocation of samples to subpopulations
3. Complexity of infection (COI) of each sample


**WARNING: Although MALECOT is live, it is still very much in development!**


The current version (0.1.0) has the following limitations:
* takes biallelic data only
* patchy documentation
* limited and basic plotting functions
* some formats/pipelines likely to change slightly

Most of these limitations should be ironed out within a few weeks of this release. See the [associated vignette](https://bobverity.github.io/MALECOT/articles/basic_tutorial.html) for a worked example analysis on simulated data.
