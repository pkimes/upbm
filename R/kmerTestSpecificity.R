#' Test k-mer Specificities
#'
#' @description
#' After estimating k-mer level affinities using probe-set aggregate, this
#' function tests for differential specificity between conditions. 
#'
#' @param se SummarizedExperiment of k-mer results from \code{kmerFit}.
#' @param span span parameter to be passed to \code{limma::loessFit}.
#'        (default = 0.05)
#' @param ... other parameters to pass to \code{limma::loessFit}.
#' 
#' @return
#' SummarizedExperiment of significance results for testing k-mer specificity
#' differences between conditions.
#'
#' @importFrom broom tidy
#' @importFrom limma loessFit
#' @importFrom stats p.adjust
#' @importFrom dplyr select group_by left_join ungroup mutate
#' @importFrom tidyr nest unnest
#' @export
#' @author Patrick Kimes
kmerTestSpecificity <- function(se, span = 0.05, ...) {

    kmers <- rowData(se)$seq

    ## gather data
    cmean <- broom::tidy(se, "contrastAverage", long = TRUE, .filter = 0L)
    cmean <- dplyr::select(cmean, seq, condition = cname, contrastAverage = value)
    cdiff <- broom::tidy(se, "contrastDifference", long = TRUE, .filter = 0L)
    cdiff <- dplyr::select(cdiff, seq, condition = cname, contrastDifference = value)
    cvar <- broom::tidy(se, "contrastVariance", long = TRUE, .filter = 0L)
    cvar <- dplyr::select(cvar, seq, condition = cname, contrastVariance = value)
    
    cdat <- dplyr::left_join(cmean, cdiff, by = c("condition", "seq"))
    cdat <- dplyr::left_join(cdat, cvar, by = c("condition", "seq"))

    cdat <- tidyr::nest(cdat, -condition)
    cdat <- dplyr::mutate(cdat,
                          fit = lapply(data, function(x) {
                              limma::loessFit(x$contrastDifference, x$contrastAverage,
                                              span = span, ...)
                          }),
                          data = mapply(function(x, y) {
                              x$contrastFit <- y$fitted
                              x$contrastResidual <- y$residuals
                              x
                          }, x = data, y = fit, SIMPLIFY = FALSE))
    cdat <- dplyr::select(cdat, -fit)
    cdat <- tidyr::unnest(cdat)

    ## compute z-scores, p-values, and adjusted p-values
    cdat <- dplyr::mutate(cdat, specificityZ = contrastResidual / sqrt(contrastVariance))
    cdat <- dplyr::mutate(cdat, specificityP = 2 * pnorm(-abs(specificityZ)))
    cdat <- dplyr::group_by(cdat, condition)
    cdat <- dplyr::mutate(cdat, specificityQ = p.adjust(specificityP, method = "BH"))
    cdat <- dplyr::ungroup(cdat)

    ## tidy results to assays
    assaylist <- list(contrastAverage = upbm:::.tidycol2mat(cdat, "contrastAverage", kmers, colnames(se)),
                      contrastDifference = upbm:::.tidycol2mat(cdat, "contrastDifference", kmers, colnames(se)),
                      contrastFit = upbm:::.tidycol2mat(cdat, "contrastFit", kmers, colnames(se)),
                      contrastResidual = upbm:::.tidycol2mat(cdat, "contrastResidual", kmers, colnames(se)),
                      specificityZ = upbm:::.tidycol2mat(cdat, "specificityZ", kmers, colnames(se)),
                      specificityP = upbm:::.tidycol2mat(cdat, "specificityP", kmers, colnames(se)),
                      specificityQ = upbm:::.tidycol2mat(cdat, "specificityQ", kmers, colnames(se)))

    rdat <- dplyr::select(cdat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]
    
    SummarizedExperiment(assays = assaylist, rowData = rdat)
}
