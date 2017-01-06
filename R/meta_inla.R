#' Fitting a pairwise meta-analysis model using INLA.
#'
#' \code{meta_inla} fits a pairwise meta-analysis model using INLA
#'
#' The following model types are supported
#' \itemize{
#'   \item \code{FE}, fixed-effect model
#'   \item \code{RE}, random effects model
#'   }
#'
#' @param datINLA An object of \code{create_INLA_dat_pair}
#' @param fixed.par A numerical vector specifying the parameter of the normal prior
#' density for mean treatment effect, first value is parameter for mean, second is for variance.
#' @param tau.prior A string specifying the prior density for the heterogeneity standard deviation,
#' options are 'uniform' for uniform prior and 'half-normal' for half-normal prior.
#' @param tau.par A numerical vector specifying the parameter of the prior
#' density for heterogenety stdev.
#' \itemize{
#'   \item \code{var.par = c(u, l)}: \code{u} is lower bound and \code{l} is upper
#'   bound when \code{var.prior = 'uniform'}
#'   \item \code{var.par = c(m, v)}: \code{m} is mean and \code{v} is variance
#'   when \code{var.prior = 'uniform'}
#' }
#' @param mreg Logical indicating whether covariate(s) should be incorporated to fit a
#' meta-regression model, default \code{FALSE}
#' @param type A string indicating the type of the model, options are "FE", "RE".
#' @param approach A string indicating the approach of the model, options are "arm-based", "contrast-based"
#' @param verbose Logical indicating whether the program should run in a verbose model, default \code{FALSE}.
#' @param inla.strategy A string specfying the strategy to use for the approximations of INLA;
#' one of 'gaussian', 'simplified.laplace' (default) or 'laplace', see \code{?INLA::control.inla}.
#' @param improve.hyperpar.dz Step length in the standardized scale used in the construction of the grid, default 0.75,
#' see \code{INLA::inla.hyperpar}.
#' @param correct Logical Add correction for the Laplace approximation, default \code{FALSE},
#' see \code{INLA::inla.hyperpar}.
#' @param correct.factor Numerical Factor used in adjusting the correction factor if \code{correct=TRUE}, default 10,
#' see \code{INLA::inla.hyperpar}.
#' @return \code{meta_inla} returns a \code{meta_inla} object with components:
#'
#' @examples
#' data('TBdat')
#' ## Create the dataset suitable for INLA
#' TBdatINLA <- create_INLA_dat_pair(TBdat$TRT, TBdat$CON, TBdat$TRTTB, TBdat$CONTB)
#'
#' ## Fitting a random-effects model using arm-based approach
#' \dontrun{
#' if(requireNamespace('INLA', quietly = TRUE)){
#'  require('INLA', quietly = TRUE)
#' fit.TB.RE.INLA <- meta_inla(TBdatINLA, type = 'RE', approach = 'arm-based',
#' tau.prior = 'uniform', tau.par = c(0, 5))
#' }
#' }
#'
#' @export
meta_inla <- function(datINLA, fixed.par = c(0, 1000),
                        tau.prior = "uniform",
                        tau.par = c(0, 5),
                        type = "FE",
                        approach = "arm-based",
                        mreg = FALSE, verbose = FALSE, inla.strategy = "simplified.laplace", improve.hyperpar.dz = 0.75,
                        correct = FALSE, correct.factor = 10)


{
  if (requireNamespace("INLA", quietly = TRUE)) {
    if (!(sum(search() == "package:INLA")) == 1) {
      stop("INLA need to be loaded! \n
           Please use the following command to load INLA,\n
           library(INLA) \n")
    }
    ################ check data
    if (!is.data.frame(datINLA$data.arm)) {
      stop("Data MUST be a data frame!!!")
    }

  if(type %in% c("FE", "RE") == FALSE){
    stop("Function argument \"type\" must be equal to \"FE\" or \"RE\"!")
  }
  if(approach %in% c("arm-based", "contrast-based") == FALSE){
    stop("Function argument \"approach\" must be equal to \"arm-based\" or \"contrast-based\"!")
  }
  if(mreg == TRUE && is.null(datINLA$data.cont$cov)){
    stop("Function argument \"cov\" must not be equal to \"NULL\" !")
  }

  inla.form <- "Y ~ -1 + d"
  if(mreg == TRUE){
    inla.form <- paste(inla.form, " + cov", sep = "")
  }
  if(approach == "arm-based"){
    inla.form <- paste(inla.form, " + mu ", sep = "")
  }
  if(tau.prior == "half-normal") {
    prior.expr <- " + f(het, model=\"iid\", prior = \"logtnormal\", param = tau.par"
  }
  if(tau.prior == "uniform") {
    # Function for Uniform distribution:
    hyperunif.function <- function(x) {
      ifelse(exp(x)^-0.5 < tau.par[2] & exp(x)^-0.5 > tau.par[1], logdens <- log(1/(tau.par[2] - tau.par[1])), logdens <- log(9.98012604599318e-322))
      logdenst <- logdens + log(0.5 * exp(-x/2))
      return(logdenst)
    }
    # Set up grid to evaluate the uniform prior:
    lprec <- seq(from = -40, to = 40, len = 20000)  ## CHANGE this LINE if INLA crashes!
    # (extend grid by changing (from= , to=))
    # Create table with prior values and lprec:
    unif.prior.table <- paste(c("table:", cbind(lprec, sapply(lprec, FUN = hyperunif.function))), sep = "", collapse = " ")
    prior.expr <- " + f(het, model=\"iid\", hyper = list(theta = list(prior = unif.prior.table))"
}
  if(type == "RE"){
    inla.form <- paste(inla.form, prior.expr , ")", sep = "")
  }
  if(approach == "contrast-based"){
    prec = datINLA$data.cont$prec
    fit.inla <- INLA::inla(stats::as.formula(inla.form), data = datINLA$data.cont, family = "normal",
                           control.fixed = list(expand.factor.strategy = "inla", mean = fixed.par[1],
                                                prec = 1/fixed.par[2]), scale = prec,
                           control.family = list(hyper = list(prec = list( fixed = TRUE, initial = 0))),
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE, config = TRUE),
                           control.inla = list(strategy = inla.strategy, correct = correct, correct.factor = correct.factor))
    if (!fit.inla$ok) {
      stop("Something wrong while running model with data! Please set verbose = TRUE to check!!!!")
    }
    # improve the estimates for hyperparameters
    fit.inla <- INLA::inla.hyperpar(fit.inla, dz = improve.hyperpar.dz)
  }
  if(approach == "arm-based"){
    Ntrials = datINLA$data.arm$sampleSize
    fit.inla <- INLA::inla(stats::as.formula(inla.form), data = datINLA$data.arm, family = "binomial",
                           control.fixed = list(expand.factor.strategy = "inla", mean = fixed.par[1],
                                                prec = 1/fixed.par[2]), Ntrials = Ntrials,
                           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE, config = TRUE),
                           control.inla = list(strategy = inla.strategy, correct = correct, correct.factor = correct.factor))
    if (!fit.inla$ok) {
      stop("Something wrong while running model with data! Please set verbose = TRUE to check!!!!")
    }
    # improve the estimates for hyperparameters
    fit.inla <- INLA::inla.hyperpar(fit.inla, dz = improve.hyperpar.dz)
  }
  nu         <- as.matrix(fit.inla$summary.fixed[1, c(1, 2, 3, 4, 5)])
  if(type == "RE"){
    tau.mean <- INLA::inla.emarginal(function(x) 1/sqrt(x), fit.inla$marginals.hyperpar[["Precision for het"]])
    tau <- INLA::inla.tmarginal(function(x) 1/sqrt(x), fit.inla$marginals.hyperpar[["Precision for het"]])
    tau.stdev = sqrt(INLA::inla.emarginal(function(x) x^2, tau) - INLA::inla.emarginal(function(x) x^1, tau)^2)
    tau.quant <- as.numeric(rev(sqrt((1/summary(fit.inla)$hyperpar[1, c(3, 4, 5)]))))
    tab <- matrix(NA, 1, 5)
    tab[1, ] <- c(tau.mean, tau.stdev, tau.quant[1], tau.quant[2], tau.quant[3])
    colnames(tab) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    rownames(tab) <- "tau"
  } else tab <- NA

  if(mreg == TRUE){
    cov <- as.numeric(fit.inla$summary.fixed[2, c(1, 2, 3, 4, 5)])
  }
  fit.inla$nu <- nu
  if (mreg == TRUE) {
    fit.inla$cov <- cov
  }
  fit.inla$tau <- tab
  fit.inla$fixed.par <- fixed.par
  fit.inla$tau.prior <- tau.prior
  fit.inla$tau.par <- tau.par
  fit.inla$type <- type
  fit.inla$approach <- approach
  fit.inla$mreg <- mreg
  fit.inla$inla.strategy <- inla.strategy

  class(fit.inla) <- "meta_inla"
  return(fit.inla)

  } else {
    stop("INLA need to be installed and loaded!\n
Please use the following command to install and load INLA,\n
install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/testing\") \n
library(INLA) \n")
  }

}

#' @export
print.meta_inla <- function (x, digits = 2, ...)
{
  if (!is.element("meta_inla", class(x)))
    stop("Argument 'x' must be an object of class \"meta.inla\".")
  cat("Time used: \n")
  print(x$cpu.used)
  if(x$mreg == FALSE){
    cat("Meta analysis using INLA\n")
    cat("Treatment effect\n")
    print(round(x$nu, digits))

  } else {
    cat("Meta regression using INLA\n")
    cat("Intercept\n")
    print(round(x$nu, digits))
    cat("Covariate\n")
    print(round(x$cov, digits))
  }
    if(x$type == "RE"){
      cat("Heterogeneity stdev\n")
      print(round(x$tau, digits))
    }

}
