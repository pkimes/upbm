#' @title Test for k-mer preferential affinities
#'
#' @description
#' Using estimated k-mer level affinities returned by \code{kmerFit}, this
#' function tests for preferential affinity within conditions for each k-mer.
#' Here, preferential affinity is defined as statistically significant
#' affinity above a normal background distribution fit to each
#' condition.
#'
#' For each condition, the k-mer affinities above background are determined by
#' first fitting a normal+exponential convolution to the distribution of estimated
#' affinities. Then, the expected exponential component for each k-mer is taken
#' as the estimate of preferential affinity. Under this formulation, the normal component
#' of the fitted convolution corresponds to the background distribution of binding
#' affinities for k-mers, while the exponential component corresponds to the signal
#' (preferential) affinity above this background.
#'
#' @param se a SummarizedExperiment of k-mer-level estimated affinities returned by
#'        \code{\link{kmerFit}}.
#' 
#' @return
#' SummarizedExperiment of k-mer preferential affinity results with the following
#' assays:
#' 
#' \itemize{
#' \item \code{"affinityEstimate"}: input k-mer affinities.
#' \item \code{"affinityVariance"}: input k-mer affinity variances.
#' \item \code{"affinitySignal"}: estimated affinity signals above background.
#' \item \code{"affinityZ"}: studentized signals (\code{affinitySignal / sqrt(affinityVariance)}).
#' \item \code{"affinityP"}: one-sided tail p-values for studentized signals.
#' \item \code{"affinityQ"}: FDR-controlling Benjamini-Hochberg adjusted p-values.
#' }
#'
#' @importFrom broom tidy
#' @importFrom limma normexp.fit normexp.signal
#' @importFrom stats p.adjust pnorm
#' @importFrom dplyr select group_by left_join ungroup mutate
#' @importFrom tidyr nest unnest
#' @export
#' @author Patrick Kimes
kmerTestAffinity <- function(se) {

    stopifnot(is(se, "SummarizedExperiment"))
    if (!all(c("affinityEstimate", "affinityVariance") %in% assayNames(se))) {
        stop("Input SummarizedExperiment is missing k-mer affinity estimates.\n",
             "SummarizedExperiment should be created by calling kmerFit(..).")
    }
    
    kmers <- rowData(se)$seq

    ## gather data
    adat <- broom::tidy(se, c("affinityEstimate", "affinityVariance"))
    adat <- dplyr::rename(adat, condition = cname)
    
    adat <- tidyr::nest(adat, data = c(-condition))
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
    adat <- tidyr::unnest(adat, cols = data)
    adat <- dplyr::mutate(adat,
                          affinityZ = affinitySignal / sqrt(affinityVariance),
                          affinityP = 2 * stats::pnorm(-abs(affinityZ)))
    adat <- dplyr::group_by(adat, condition)
    adat <- dplyr::mutate(adat, affinityQ = stats::p.adjust(affinityP, method = "BH"))
    adat <- dplyr::ungroup(adat)

    ## tidy results to assays
    assaylist <- list(affinityEstimate = .tidycol2mat(adat, "affinityEstimate", kmers, colnames(se)),
                      affinityVariance = .tidycol2mat(adat, "affinityVariance", kmers, colnames(se)),
                      affinitySignal = .tidycol2mat(adat, "affinitySignal", kmers, colnames(se)),
                      affinityZ = .tidycol2mat(adat, "affinityZ", kmers, colnames(se)),
                      affinityP = .tidycol2mat(adat, "affinityP", kmers, colnames(se)),
                      affinityQ = .tidycol2mat(adat, "affinityQ", kmers, colnames(se)))

    rdat <- dplyr::select(adat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]
    
    SummarizedExperiment(assays = assaylist, rowData = rdat)
}
