#' Fitting a network meta-analysis model using INLA
#'
#' \code{nma_inla} fits a network meta-analysis model using INLA.
#'
#' The following likelihood types are supported
#' \itemize{
#'   \item \code{normal}: for continuous (mean difference) data.
#'
#'   Required coloumns: \code{[mean, std.err]}
#'
#'   Result: relative mean difference
#'   \item \code{binomial}: for dichotomous data.
#'
#'   Required coloumns: \code{[responders, sampleSize]}
#'
#'   Result: log odds ratio
#'   \item \code{normal}: for event-rate (survival) data.
#'
#'   Required coloumns: \code{[responders, exposure]}
#'
#'   Result: log hazard ratio
#'   }

#' The following model types are supported
#' \itemize{
#'   \item \code{FE}, ordinary fixed effect model, assuming homogeneity between trials
#'   (Dias et al., 2013)
#'   \item \code{consistency}, ordinary consistency model, assuming consistency in the
#'   network. (Jackson et al., 2014)
#'   \item \code{jackson}, the design-by-treatment interaction model with random
#'   inconsistency parameters. (Jackson et al., 2014)
#'   }
#'
#' @param datINLA An object of \code{create_INLA_dat}
#' @param likelihood The likelihood to be used.
#' @param fixed.par A numerical vector specifying the parameter of the normal prior
#' density for basic parameters, first value is parameter for mean, second is for variance.
#' @param tau.prior A string specifying the prior density for the heterogeneity standard deviation,
#' options are 'uniform' for uniform prior and 'half-normal' for half-normal prior.
#' @param tau.par A numerical vector specifying the parameter of the prior
#' density for heterogenety stdev.
#' \itemize{
#'   \item \code{var.par = c(u, l)}: \code{u} is lower bound and \code{l} is upper
#'   bound when \code{var.prior = 'uniform'}.
#'   \item \code{var.par = c(m, v)}: \code{m} is mean and \code{v} is variance
#'   when \code{var.prior = 'half-normal'}.
#' }
#' @param kappa.prior A string specifying the prior density for the inconsistency standard deviation,
#' options are 'uniform' for uniform prior and 'half-normal' for half-normal prior.
#' @param kappa.par A numerical vector specifying the parameter of the prior.
#' density for inconsistency stdev. The definition of the priors is the same as for \code{tau.par}.
#' @param mreg Logical indicating whether covariate(s) should be incorporated to fit a
#' network meta-regression model, default \code{FALSE}.
#' @param type A string indicating the type of the model, options are "FE", "consistency" and "jackson".
#' @param verbose Logical indicating whether the program should run in a verbose model, default \code{FALSE}.
#' @param inla.strategy A string specfying the strategy to use for the approximations of INLA;
#' one of 'gaussian', 'simplified.laplace' (default) or 'laplace', see \code{?INLA::control.inla}.
#' @param improve.hyperpar.dz Step length in the standardized scale used in the construction of the grid, default 0.75,
#' see \code{INLA::inla.hyperpar}. Not used if \code{mod = 'FE'}.
#' @param correct Logical Add correction for the Laplace approximation, default \code{FALSE},
#' see \code{INLA::inla.hyperpar}. Not used if \code{mod = 'FE'}.
#' @param correct.factor Numerical Factor used in adjusting the correction factor if \code{correct=TRUE}, default 10,
#' see \code{INLA::inla.hyperpar}. Not used if \code{mod = 'FE'}.
#' @return \code{nma_inla} returns a \code{nma_inla} object.
#'
#' @examples
#' SmokdatINLA <- create_INLA_dat(dat = Smokdat, armVars = c('treatment' = 't', 'responders' = 'r'
#' ,'sampleSize' = 'n'), nArmsVar = 'na')
#' \dontrun{
#' ## Fitting a consistency model
#' if(requireNamespace('INLA', quietly = TRUE)){
#'  require('INLA', quietly = TRUE)
#' fit.Smok.cons.INLA <- nma_inla(SmokdatINLA, likelihood = 'binomial', type = 'FE',
#' tau.prior = 'uniform', tau.par = c(0, 5))
#' }
#' }
#'
#' @export
nma_inla <- function(datINLA, likelihood = NULL, fixed.par = c(0, 1000), tau.prior = "uniform",
                     tau.par = c(0, 5), kappa.prior = "uniform", kappa.par = c(0, 5), mreg = FALSE,
                     type = "consistency", verbose = FALSE, inla.strategy = "simplified.laplace", improve.hyperpar.dz = 0.75,
                     correct = FALSE, correct.factor = 10)
{
    if (requireNamespace("INLA", quietly = TRUE)) {
        if (!(sum(search() == "package:INLA")) == 1) {
            stop("INLA need to be loaded! \n
           Please use the following command to load INLA,\n
           library(INLA) \n")
        }
        ################ check data
        if (!is.data.frame(datINLA)) {
            stop("Data MUST be a data frame!!!")
        }
        ################ check likelihood type
        if (likelihood %in% c("binomial", "normal", "poisson") == FALSE) {
          stop("Function argument \"likelihood\" must be equal to \"binomial\" or \"normal\" or \"poisson\" !!!")
        }
        ################ check model type
        if (type %in% c("FE", "consistency", "jackson") == FALSE) {
            stop("Function argument \"type\" must be equal to \"FE\" or \"consistency\" or \"jackson\" !!!")
        }
        ################ check priors for hyperparameters
        if (tau.prior %in% c("uniform", "half-normal") == FALSE) {
            stop("Function argument \"tau.prior\" must be equal to \"uniform\" or \"half-normal\" !!!")
        }
        if (kappa.prior %in% c("uniform", "half-normal") == FALSE) {
            stop("Function argument \"kappa.prior\" must be equal to \"uniform\" or \"half-normal\" !!!")
        }

        cor <- 0.5  # correlation between treatment comparisons of the same multi-arm trial.
        ngroup <- max(datINLA$na) - 1
        cor.inla.init <- log((1 + cor * (ngroup - 1))/(1 - cor))
        N <- max(datINLA$study)
        d_params <- grep("d1", names(datINLA), value = TRUE)
        N_d_params <- length(d_params)
        inla.form <- paste("Y ~ -1 + mu +", paste(d_params, collapse = "+", sep = " "), sep = " ")
        if (mreg == TRUE) {
          inla.form <- paste(inla.form, " + cov", sep = "")
        }
        if (type %in% c("consistency", "jackson") == TRUE) {
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

            if (tau.prior == "uniform") {
                het.expr <- " + f(het, model=\"iid\", hyper = list(theta = list(prior = unif.prior.table)"
            }
            if (tau.prior == "half-normal") {
                # INLA uses precisions instead of variances
                tau.par[2] <- 1 / tau.par[2]
                het.expr <- " + f(het, model=\"iid\", hyper = list(theta = list(prior = \"logtnormal\", param = tau.par)"
            }

            multi.arm.expr <- ", group = g, control.group = list(model = \"exchangeable\", hyper = list(rho = list(fixed = TRUE, initial = cor.inla.init)))"
            inla.form <- paste(inla.form, het.expr, ")", sep = "")
            inla.form <- paste(inla.form, multi.arm.expr, ")", sep = "")
        }
        if (type == "jackson") {
            if (kappa.prior == "uniform") {
                # Function for Uniform distribution:
                inc.expr <- " + f(inc, model=\"iid\", hyper = list(theta = list(prior = unif.prior.table)"
            }
            if (kappa.prior == "half-normal") {
              kappa.par[2] <- 1 / kappa.par[2]
              inc.expr <- " + f(inc, model=\"iid\", hyper = list(theta = list(prior = \"logtnormal\", param = kappa.par)"
            }
            inla.form <- paste(inla.form, inc.expr, ")", sep = "")
            inla.form <- paste(inla.form, multi.arm.expr, ")", sep = "")
        }
        # inla call
        if (likelihood == "binomial") {
          datINLA$Y = datINLA$responders
          Ntrials = datINLA$sampleSize
          fit.inla <- INLA::inla(stats::as.formula(inla.form), data = datINLA, family = "binomial", verbose = verbose,
                               control.fixed = list(expand.factor.strategy = "inla", mean = fixed.par[1],
                                                    prec = 1/fixed.par[2]), Ntrials = Ntrials,
                               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE, config = TRUE),
                               control.inla = list(strategy = inla.strategy, correct = correct, correct.factor = correct.factor))
        }
        if (likelihood == "normal") {
          datINLA$Y = datINLA$mean
          prec = 1 / datINLA$std.err^2
          fit.inla <- INLA::inla(stats::as.formula(inla.form), family = "normal", verbose = verbose,
                                 data = datINLA,
                                 control.fixed = list(expand.factor.strategy = "inla", mean = fixed.par[1],
                                                      prec = 1/fixed.par[2]),
                                 control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                                 scale = prec,
                                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE, config = TRUE),
                                 control.inla = list(strategy = inla.strategy, correct = correct, correct.factor = correct.factor))
        }
        if (likelihood == "poisson") {
          datINLA$Y = datINLA$responders
          E = datINLA$exposure
          fit.inla <- INLA::inla(stats::as.formula(inla.form), data = datINLA, family = "poisson", verbose = verbose,
                                 control.fixed = list(expand.factor.strategy = "inla", mean = fixed.par[1],
                                                      prec = 1/fixed.par[2]), E = E,
                                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE, config = TRUE),
                                 control.inla = list(strategy = inla.strategy, correct = correct, correct.factor = correct.factor))
        }

        if (!fit.inla$ok) {
            stop("Something wrong while running model with data! Please set verbose = TRUE to check!!!!")
        }
        if (type %in% c("consistency", "jackson") == TRUE) {
        # improve the estimates for hyperparameters
        fit.inla <- INLA::inla.hyperpar(fit.inla, dz = improve.hyperpar.dz)
        }
        # Extracting estiamates of basic parameters from INLA object
        d_params <- as.matrix(fit.inla$summary.fixed[(N + 1):(N + N_d_params), c(1, 2, 3, 4, 5)])
        if (mreg == TRUE) {
        # Extracting estiamates of basic parameters
          cov <- as.matrix(fit.inla$summary.fixed[N + N_d_params + 1, c(1, 2, 3, 4, 5)])
        }
        # Estimates of hyperparameters need to be transformed
        if (type %in% c("consistency", "jackson") == TRUE) {
            tau.mean <- INLA::inla.emarginal(function(x) 1/sqrt(x), fit.inla$marginals.hyperpar[["Precision for het"]])
            tau <- INLA::inla.tmarginal(function(x) 1/sqrt(x), fit.inla$marginals.hyperpar[["Precision for het"]])
            tau.stdev = sqrt(INLA::inla.emarginal(function(x) x^2, tau) - INLA::inla.emarginal(function(x) x^1, tau)^2)
            tau.quant <- as.numeric(rev(sqrt((1/summary(fit.inla)$hyperpar[1, c(3, 4, 5)]))))
            if (type == "jackson") {
                kappa.mean <- INLA::inla.emarginal(function(x) 1/sqrt(x), fit.inla$marginals.hyperpar[["Precision for inc"]])
                # this should be fixed!
                try(kappa <- INLA::inla.tmarginal(function(x) 1/sqrt(x), fit.inla$marginals.hyperpar[["Precision for inc"]]), silent = TRUE)
                if (is.numeric(kappa) == TRUE) {
                  kappa.stdev = sqrt(INLA::inla.emarginal(function(x) x^2, kappa) - INLA::inla.emarginal(function(x) x^1, kappa)^2)
                }
                kappa.quant <- as.numeric(rev(sqrt((1/summary(fit.inla)$hyperpar[2, c(3, 4, 5)]))))
                tab <- matrix(NA, 2, 5)
                tab[1, ] <- c(tau.mean, tau.stdev, tau.quant[1], tau.quant[2], tau.quant[3])
                if (is.numeric(kappa) == TRUE) {
                  tab[2, ] <- c(kappa.mean, kappa.stdev, kappa.quant[1], kappa.quant[2], kappa.quant[3])
                } else tab[2, ] <- c(kappa.mean, NA, kappa.quant[1], kappa.quant[2], kappa.quant[3])
            } else {
              tab <- matrix(NA, 1, 5)
              rownames(tab) <- "tau"
              }
            tab[1, ] <- c(tau.mean, tau.stdev, tau.quant[1], tau.quant[2], tau.quant[3])
            colnames(tab) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
        } else tab <- NA

        fit.inla$d_params <- d_params
        if (mreg == TRUE) {
          fit.inla$cov <- cov
        }
        fit.inla$hyperpar <- tab
        fit.inla$fixed.par <- fixed.par
        fit.inla$tau.prior <- tau.prior
        fit.inla$tau.par <- tau.par
        fit.inla$kappa.prior <- kappa.prior
        fit.inla$kappa.par <- kappa.par
        fit.inla$type <- type
        fit.inla$mreg <- mreg
        fit.inla$inla.strategy <- inla.strategy
        fit.inla$N <- N
        fit.inla$N_d_params <- N_d_params

        class(fit.inla) <- "nma_inla"
        return(fit.inla)
    } else {
        stop("INLA need to be installed and loaded!\n
Please use the following command to install and load INLA,\n
install.packages(\"INLA\", repos=\"http://www.math.ntnu.no/inla/R/testing\") \n
library(INLA) \n")
    }
}
