#' nmainla: Network meta-analysis using INLA
#'
#' An R package for performing network meta-analysis using INLA.
#'
#' @section Details:
#' Network meta-analysis is a generalization of pairwise meta-analysis to analyze networks
#' of trials comparing two or more treatments simultaneously (Dias et al, 2011). Bayesian
#' hierarchical models are commonly used for network meta-analysis (Dias et al, 2011). The
#' default choice for performing inference within such models are Markov Chain Monte Carlo
#' (MCMC), for example using BUGS-variants programs such as JAGS. A deterministic approach
#' to do fully Bayesian inference for latent Gaussian models (LGMs) are integrated nested
#' Laplace approximations (INLA) (Rue et al, 2009) which is a fast and accurate alternative
#' to MCMC. INLA methodology is implemented as an R package INLA (<www.r-inla.org>). Sauter
#' and Held (2015) has shown that INLA can be used for fitting many NMA models including
#' fixed effect and consistency models, node-splitting models.
#'
#' This package extends the INLA implementation of Sauter and Held (2015) to Jackson
#' model (Jackson et al, 2014) and network meta-regression and extracts the
#' features needed for NMA models from INLA R package and presents
#' in an intuitive way (Guenhan et al, in preparation). Currently, contrast-based network
#'  meta-analysis using trial-arm level data for datasets with binary, continuous, and survival outcomes are supported.
#'  Note that the installation of R package
#' 'INLA' is compulsory for successful usage. The 'INLA' package can be obtained from
#' <http://www.r-inla.org>. We recommend the testing version, which can be downloaded
#' by running: source("http://www.math.ntnu.no/inla/givemeINLA-testing.R").
#'
#' Type vignette("nmainla") to how to use this package.
#'
#' The development version of nmainla is available on GitHub <https://github.com/gunhanb/nmainla>.
#'
#'
#' @source Sauter, R. and Held, L. (2015). Network meta-analysis
#' with integrated nested Laplace approximations. Biometrical Journal 57 1038--1050.
#' @source Guenhan, B.K., Friede, T., Held, L. (in preparation)
#' A design-by-treatment interactoion model for network meta-analysis
#' using integrated nested Laplace approximations
#' @source Jackson, D., Barrett, J. K., Rice, S., White, I. R. and Higgins, J. P. (2014).
#' A design-by-treatment interaction model for network meta-analysis with random
#' inconsistency effects. Statistics in Medicine 33 3639--3654.
#' @source Rue, H., Martino, S. and Chopin, N. (2009). Approximate Bayesian inference
#' for latent Gaussian models by using integrated nested Laplace approximations.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology) 71 319--392.
#' @source Dias, S., Welton, N. J., Sutton, A. J. and Ades, A. (2011). NICE DSU
#' Technical Support Document 2: A Generalised Linear Modelling Framework for Pairwise
#' and Network Meta-analysis of Randomised Controlled Trials. Last updated September 2016.
#' @source Dias, S., Sutton, A. J., Welton, N. J. and Ades, A. E. (2013). Evidence
#' synthesis for Decision Making 3: Heterogeneity--Subgroups, Meta-Regression, Bias,
#' and Bias-Adjustment. Medical Decision Making 33 618--640.
#' @author Burak Kuersad Guenhan <burak.gunhan@med.uni-goettingen.de>
#' @docType package
#' @name nmainla
NULL
