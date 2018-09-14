# MALECOT
### Version 0.1.0
[![Travis-CI Build Status](https://travis-ci.org/bobverity/MALECOT.svg?branch=master)](https://travis-ci.org/bobverity/MALECOT)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/bobverity/MALECOT?branch=master&svg=true)](https://ci.appveyor.com/project/bobverity/MALECOT)
[![Coverage Status](https://img.shields.io/codecov/c/github/bobverity/MALECOT/master.svg)](https://codecov.io/github/bobverity/MALECOT?branch=master)
[![Documentation](https://github.com/bobverity/MALECOT/blob/master/R_ignore/images/documentation-click%20here!-blue.png)](https://bobverity.github.io/MALECOT/)

--------------------------------------------------------------------------------------------------------------------------------

MALECOT is an R package for inferring population structure from genetic data, and is tailored to problems in malaria where each sample commonly contains multiple genotypes. The basic model assumes that there are *K* subpopulations represented in the data, and that we are interested in working out which samples belong to which subpopulations. The analysis works via Bayesian Markov chain Monte Carlo, alternately estimating each of three unknown parameters:

1. Allele frequencies in all subpopulations
2. Complexity of infection (COI) of each sample, i.e. how many genotypes are represented in each sample
3. Group allocation, i.e. which samples belong to which subpopulations

To get started, take a look at the [installation instructions](https://bobverity.github.io/MALECOT), followed by a [basic tutorial](https://bobverity.github.io/MALECOT) on running the program.

