#' Test k-mer Affinities
#'
#' @description
#' After estimating k-mer level affinities using probe-set aggregate, this
#' function tests for high affinity within individual conditions. 
#'
#' @param se SummarizedExperiment of k-mer results from \code{kmerFit}.
#' 
#' @return
#' SummarizedExperiment of significance results for testing k-mer affinity
#' strength within conditions.
#'
#' @importFrom broom tidy
#' @importFrom limma normexp.fit normexp.signal
#' @importFrom stats p.adjust
#' @importFrom dplyr select group_by left_join ungroup mutate
#' @importFrom tidyr nest unnest
#' @export
#' @author Patrick Kimes
kmerTestAffinity <- function(se) {

    kmers <- rowData(se)$seq

    ## gather data
    aest <- broom::tidy(se, "affinityEstimate", long = TRUE, .filter = 0L)
    aest <- dplyr::select(aest, seq, condition = cname, affinityEstimate = value)
    avar <- broom::tidy(se, "affinityVariance", long = TRUE, .filter = 0L)
    avar <- dplyr::select(avar, seq, condition = cname, affinityVariance = value)

    adat <- dplyr::left_join(aest, avar, by = c("condition", "seq"))
    
    adat <- tidyr::nest(adat, -condition)
    adat <- dplyr::mutate(adat,
                          nefit = lapply(data, function(x) {
                              limma::normexp.fit(x$affinityEstimate, method = "saddle")$par
                          }))
    adat <- dplyr::mutate(adat,
                          data = mapply(function(x, y) {
                              x$affinitySignal <- limma::normexp.signal(par = y, x = x$affinityEstimate)
                              x
                          }, x = data, y = nefit, SIMPLIFY = FALSE))
    adat <- dplyr::select(adat, -nefit)
    adat <- tidyr::unnest(adat)
    adat <- dplyr::mutate(adat,
                          affinityZ = affinitySignal / sqrt(affinityVariance),
                          affinityP = 2*pnorm(-abs(affinityZ)))
    adat <- dplyr::group_by(adat, condition)
    adat <- dplyr::mutate(adat, affinityQ = p.adjust(affinityP, method = "BH"))
    adat <- dplyr::ungroup(adat)

    ## tidy results to assays
    assaylist <- list(affinityEstimate = upbm:::.tidycol2mat(adat, "affinityEstimate", kmers, colnames(se)),
                      affinityVariance = upbm:::.tidycol2mat(adat, "affinityVariance", kmers, colnames(se)),
                      affinitySignal = upbm:::.tidycol2mat(adat, "affinitySignal", kmers, colnames(se)),
                      affinityZ = upbm:::.tidycol2mat(adat, "affinityZ", kmers, colnames(se)),
                      affinityP = upbm:::.tidycol2mat(adat, "affinityP", kmers, colnames(se)),
                      affinityQ = upbm:::.tidycol2mat(adat, "affinityQ", kmers, colnames(se)))

    rdat <- dplyr::select(adat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]
    
    SummarizedExperiment(assays = assaylist, rowData = rdat)
}
