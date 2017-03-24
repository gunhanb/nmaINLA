#' Print nmainla object
#'
#' Takes an \code{nma_inla} object which is obtained by function \code{nma_inla} and print
#' the model and data information such as model type used in the model.
#'
#' The resulting data.frame can be used as data argument in \code{nma_inla}.
#'
#' @param x An \code{nma_inla} object.
#' @param digits An integer indicating the number of decimal places.
#' @param ... Further arguments passed to or from other methods.
#' @return The return value is invisible \code{NULL}
#' @export
print.nma_inla <- function(x, digits = 3, ...) {
    cat("Network meta-analysis using INLA\n")
    cat("Relative treatment effects\n")
    print(round(x$d_params, digits))
    if (x$mreg == TRUE) {
      cat("Covariate\n")
      print(round(x$cov, digits))
    }
    if (x$type == "consistency") {
        cat("Heterogeneity stdev\n")
        print(round(x$hyper[1, ], digits))
    } else if (x$type == "jackson") {
        cat("Heterogeneity stdev\n")
        print(round(x$hyperpar[1, ], digits))
        cat("Inconsistency stdev\n")
        print(round(x$hyperpar[2, ], digits))
    }
    return(invisible())
}
