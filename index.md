# MALECOT
### Version 0.1.1
[![Travis-CI Build Status](https://travis-ci.org/bobverity/MALECOT.svg?branch=master)](https://travis-ci.org/bobverity/MALECOT)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/bobverity/MALECOT?branch=master&svg=true)](https://ci.appveyor.com/project/bobverity/MALECOT)
[![Coverage Status](https://img.shields.io/codecov/c/github/bobverity/MALECOT/master.svg)](https://codecov.io/github/bobverity/MALECOT?branch=master)
[![Documentation](https://github.com/bobverity/MALECOT/blob/master/R_ignore/images/documentation-click%20here!-blue.png)](https://bobverity.github.io/MALECOT/)

--------------------------------------------------------------------------------------------------------------------------------

MALECOT is an R package for inferring population structure from genetic data. It is tailored to problems in malaria where each sample commonly contains multiple genotypes. The basic model assumes that there are *K* subpopulations represented in the data, and that we are interested in working out which samples belong to which subpopulations, similar to programs like [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html). The analysis works via Bayesian Markov chain Monte Carlo, alternately estimating each of three unknown parameters:

1. Allele frequencies in all subpopulations
2. Group allocation, i.e. which samples belong to which subpopulations
3. Complexity of infection (COI) of each sample, i.e. how many genotypes are represented in each sample

The first two parameters are estimated by programs like STRUCTURE, but the third parameter is not, making STRUCTURE not ideally suited to malaria data where complex infections are common. There are also more complex priors available in MALECOT, and the ability to estimate to estimate *K* via a technique called generalised thermodynamic integration (GTI).

To get started, take a look at the [installation instructions](https://bobverity.github.io/MALECOT/articles/installation.html), followed by a [basic tutorial](https://bobverity.github.io/MALECOT/articles/tutorial-biallelic.html) on running the program.


