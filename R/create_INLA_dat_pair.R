#' Prepare pairwise meta-analysis dataset for INLA.
#'
#' \code{create_INLA_dat_pair} creates two dataframes, one to use in
#' a contrast based and the other in an arm-based pairwise meta-analysis.
#'
#' The resulting data.frame can be used as data argument in
#' \code{meta_inla}.
#'
#'
#' @param ntrt Number of subjects in treatment arm
#' @param nctrl Number of subjects in control arm
#' @param ptrt Number of events in treatment arm
#' @param pctrl Number of events in treatment arm
#' @param cov Optional argument to include a covariate in the model
#' @return A list of two dataframe objects
#' @examples
#' data('TBdat')
#' ## Create the dataset suitable for INLA
#' TBdatINLA <- create_INLA_dat_pair(TBdat$TRT, TBdat$CON, TBdat$TRTTB, TBdat$CONTB)
#'
#' ## Check that the data are correct
#' print(TBdatINLA)
#' @export
create_INLA_dat_pair <- function(ntrt, nctrl, ptrt, pctrl, cov = NULL) {
    # Continuity correction incase of zero events
    zerocell1 <- function(y) {
        if (y["ptrt"] == 0) {
            y["ptrt"] <- 0.5
            y["ntrt"] <- y["ntrt"] + 1
        }
        return(y)
    }

    zerocell2 <- function(y) {
        if (y["pctrl"] == 0) {
            y["pctrl"] <- 0.5
            y["ntrt"] <- y["ntrt"] + 1
        }
        return(y)
    }

    data <- NULL
    data$ptrt <- ptrt
    data$ntrt <- ntrt
    data$pctrl <- pctrl
    data$nctrl <- nctrl
    data <- data.frame(data)
    N <- nrow(data)
    data.nozero <- t(apply(data, 1, zerocell1))
    data.nozero <- t(apply(data.nozero, 1, zerocell2))
    d <- rep(1, times = N)
    Y <- apply(data.nozero, 1, function(x) log((x[1] * (x[4] - x[3]))/((x[2] - x[1]) * x[3])))  # observed log odds ratios
    prec <- 1/apply(data.nozero, 1, function(x) 1/x[1] + 1/(x[2] - x[1]) + 1/x[3] + 1/(x[4] - x[3]))  # observed precisions
    het <- 1:nrow(data.nozero)  # ID for random effects
    # Dataset for summary-based meta-analysis
    data.cont <- data.frame(cbind(d, Y, prec, het, cov))
    Y <- as.vector(rbind(data$pctrl, data$ptrt))  # number of events
    sampleSize <- as.vector(rbind(data$nctrl, data$ntrt))  # number of all patients
    d <- rep(0:1, times = N)
    het <- as.vector(rbind(rep(NA, times = N), 1:N))  # ID for random effects
    if (!is.null(cov)) {
        cov <- as.vector(rbind(NA, cov))
    }
    # Dataset for arm-levelmeta-analysis
    data.arm <- data.frame(cbind(Y, sampleSize, d, het, cov))
    data.arm$mu <- as.factor(as.numeric(gl(n = N, k = 2)))
    datINLA <- list(data.cont = data.cont, data.arm = data.arm)
    return(datINLA)
}
