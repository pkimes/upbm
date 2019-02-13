#' @title Test for k-mer preferential affinities
#'
#' @description
#' Using estimated k-mer level affinities returned by \code{kmerFit}, this
#' function tests for preferential affinity across k-mers in each condition
#' separately. Here, preferential affinity is defined as statistically significant
#' affinity above a normal background affinity distribution fit to each
#' condition.
#'
#' For each condition, the k-mer affinities above background are determined by
#' fitting a normal+exponential convolution to the distribution of estimated
#' affinities and taking the expected exponential component for each
#' k-mer. Intuitively, the normal component of the fitted convolution captures
#' the background distribution of binding affinities for k-mers, while the exponential
#' component captures any signal (preferential) affinity above this
#' background.
#'
#' @param se a SummarizedExperiment of k-mer-level estimated affinities returned by
#'        \code{\link{kmerFit}}.
#' 
#' @return
#' SummarizedExperiment of k-mer affinity testing results with the following
#' assays:
#' 
#' \itemize{
#' \item \code{"affinityEstimate"}: input k-mer affinities.
#' \item \code{"affinityVariance"}: input k-mer affinity variances.
#' \item \code{"affinitySignal"}: estimated affinity signals above background.
#' \item \code{"affinityZ"}: studentized signals (signal / sqrt(variance)).
#' \item \code{"affinityP"}: one-sided tail p-values for studentized signals .
#' \item \code{"affinityQ"}: FDR-controlling Benjamini-Hochberg adjusted p-values.
#' }
#'
#' @importFrom broom tidy
#' @importFrom limma normexp.fit normexp.signal
#' @importFrom stats p.adjust
#' @importFrom dplyr select group_by left_join ungroup mutate
#' @importFrom tidyr nest unnest
#' @export
#' @author Patrick Kimes
kmerTestAffinity <- function(se) {

    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot(c("affinityEstimate", "affinityVariance") %in%
              assayNames(se))

    kmers <- rowData(se)$seq

    ## gather data
    aest <- broom::tidy(se, "affinityEstimate", long = TRUE)
    aest <- dplyr::select(aest, seq, condition = cname, affinityEstimate = value)
    avar <- broom::tidy(se, "affinityVariance", long = TRUE)
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
