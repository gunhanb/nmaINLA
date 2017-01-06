#' Data for the smoking cessation NMA-network discussed in Dias et al. (2010)
#'
#' The smoking cessation dataset describes a network comparing 4 treatments with 24 trials.
#' There are 22 pairwise comparisons and 2 thee-arm trials.
#'
#' @format A data frame with 24 observations on the following 10 variables.
#' \describe{
#' \item{\code{r1}}{Number of events (responses) in the baseline treatment (treatment
#' 1).}
#' \item{\code{r2}}{Number of events (responses) in the second study arm treatment
#' (treatment 2).}
#' \item{\code{r3}}{Number of events (responses) in the third study arm treatment
#' (treatment 3).}
#' \item{\code{n1}}{Total number of study participants in the baseline
#' treatment (treatment 1).}
#' \item{\code{n2}}{Total number of study participants in the second study arm
#' treatment (treatment 2).}
#' \item{\code{n3}}{Total number of study participants in the second study arm
#' treatment (treatment 3).}
#' \item{\code{t1}}{Indicator variable identifying treatment 1.}
#' \item{\code{t2}}{Indicator variable identifying treatment 2.}
#' \item{\code{t3}}{Indicator variable identifying treatment 3.}
#' \item{\code{na}}{Indicator with number of arms in trial.}
#' \item{\code{des}}{Design, the set of treatments included in each study}
#' }
#'
#' @source Dias, S., Welton, N. J., Caldwell, D. M. and Ades, A. E. (2010)
#' Checking consistency in mixed treatment comparison meta-analysis,
#' Statistics in Medicine, 29:932--944.
"Smokdat"
