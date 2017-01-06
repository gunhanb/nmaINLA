#' Prepare network meta-analysis dataset for INLA.
#'
#' \code{create_INLA_dat} converts datasets in the one-study-per-row format to one-arm-per-row format
#', then adds indicator (dummy) variables for the basic contrasts,
#' heterogeneity random effects and design-specific inconsistency random effects and for
#' correlated multi-arm trials.
#'
#' The resulting data.frame can be used as data argument in \code{nma_inla}.
#'
#' @param dat Data in one-study-per-row format.
#' @param armVars Vector of per-arm variables
#'  The name of each component will be the column name in the resulting dataset.
#' @param covariate Optional. Vector of study-specific covariate
#' @param design Optional. Vector of study-specific design. We refer design for
#' the set of treatments in each trial.
#' @param nArmsVar Variable holding the number of arms for each study.
#' @return A data frame with the generated coloumns.
#' @author Burak Kursad Gunhan, \email{burak.gunhan@med.uni-goettingen.de}
#' and Gert van Valkenhoef, \email{g.h.m.van.valkenhoef@rug.nl}
#' @seealso \code{gemtc::mtc.data.studyrow}
#' @examples
#' data('Smokdat')
#' ## Create the dataset suitable for INLA
#' SmokdatINLA <- create_INLA_dat(dat = Smokdat, armVars = c('treatment' = 't', 'responders' = 'r'
#' ,'sampleSize' = 'n'), nArmsVar = 'na')
#' ## Check that the data are correct
#' print(SmokdatINLA)
#' @export
create_INLA_dat <- function(dat = dat,
                            armVars = c(treatment = "t", responders = "r", sampleSize = "n"),
                            covariate = "cov",
                            design = "des",
                            nArmsVar = "na") {
    ######################################################## THIS CODE IS COPY PASTE FROM gemtc::mtc.data.studyrow
    studyNames = 1:nrow(dat)
    patterns = c("%s", "%s%d")
    studyVars = c()
    treatmentNames = NA
    colsOrNA <- function(row, cols) {
        rval <- rep(NA, length(cols))
        sel <- cols %in% colnames(row)
        rval[sel] <- row[cols[sel]]
        rval
    }

    nArmsCol <- sprintf(patterns[1], nArmsVar)
    studyCols <- sprintf(patterns[1], studyVars)

    out <- do.call(rbind, lapply(1:nrow(dat), function(i) {
        row <- dat[i, ]
        na <- row[nArmsCol]
        studyEntries <- row[studyCols]
        names(studyEntries) <- names(studyVars)
        do.call(rbind, lapply(1:unlist(na), function(j) {
            armCols <- sprintf(patterns[2], armVars, j)
            armEntries <- colsOrNA(row, armCols)
            names(armEntries) <- names(armVars)
            c(list(study = i), studyEntries, armEntries)
        }))
    }))

    colNames <- colnames(out)
    out <- lapply(colnames(out), function(col) {
        unlist(out[, col])
    })
    names(out) <- colNames

    out[["study"]] <- studyNames[out[["study"]]]
    if (all(!is.na(treatmentNames))) {
        out[["treatment"]] <- treatmentNames[out[["treatment"]]]
    }
    datINLA <- do.call(data.frame, out)
    #####################      mtc.data.studyrow is finished Adding indicator variable needed for INLA
    datINLA$design <- NULL
    datINLA$covariate <- NULL
    datINLA$na <- rep(dat[[paste(nArmsVar)]], times = dat[[paste(nArmsVar)]])
    N <- nrow(datINLA)
    datINLA$baseline <- rep(dat[[paste(armVars[[1]], 1, sep = "")]], times = dat[[paste(nArmsVar)]])
    # Study must be as factor!!!
    datINLA$mu <- as.factor(datINLA$study)
    # Create basic parameters (compared to treatment 1)
    unibase <- sort(unique(datINLA[["treatment"]]))
    for (k in 2:length(unibase)) {
        datINLA[[paste("d", 1, unibase[k], sep = "")]] <- rep(0, N)
        datINLA[[paste("d", 1, unibase[k], sep = "")]][which(datINLA[["treatment"]] == k & datINLA[["baseline"]] == 1)] <- 1
    }
    if(length(unibase) != 2) {
      for (k in 2:(length(unibase) - 1)) {
        for (j in (k + 1):length(unibase)) {
            datINLA[which(datINLA[["treatment"]] == j & datINLA[["baseline"]] == k), paste("d", 1, unibase[k], sep = "")] <- -1
            datINLA[which(datINLA[["treatment"]] == j & datINLA[["baseline"]] == k), paste("d", 1, unibase[j], sep = "")] <- 1
        }
      }
    }
    # Adding indicator variable for multi-arm trials 'g'
    g <- c()
    for (i in 1:length(unique(datINLA[["study"]]))) {
        g <- c(g, c(NA, seq(1:(table(datINLA[["study"]])[i] - 1))))
    }
    datINLA$g <- g
    # Adding indicator variable for trial-specific heterogeneity random effetcs 'het'
    het <- datINLA[["study"]]
    k <- 1
    for (i in 1:length(unique(datINLA[["study"]]))) {
        het[k] <- NA
        k <- sum(table(datINLA[["study"]])[1:i]) + 1
    }
    datINLA$het <- het
    # Adding indicator variable for design-specific inconsistency random effetcs 'inc'
    if (is.null(dat[[paste(design)]]) == FALSE) {
        inc <- c()
        inc <- rep(dat[[paste(design)]], times = dat[[paste(nArmsVar)]])
        for (i in 1:nrow(datINLA)) {
            if (datINLA$het[i] %in% NA)
                inc[i] <- NA
        }
        datINLA$inc <- inc
    }
    # Adding covariate
    if (is.null(dat[[paste(covariate)]]) == FALSE) {
      cov <- c()
      cov <- rep(dat[[paste(covariate)]], times = dat[[paste(nArmsVar)]])
      for (i in 1:nrow(datINLA)) {
        if (datINLA$het[i] %in% NA && datINLA$baseline[i] %in% 1)
          cov[i] <- NA
      }
      datINLA$cov <- cov
    }

    return(datINLA)
}

