#' @title Test for k-mer differential specificities
#'
#' @description
#' Using estimated k-mer level affinities returned by \code{kmerFit}, this
#' function tests for differential specificities across conditions for each k-mer
#' separately. The call to \code{kmerFit} must have been made with
#' \code{contrasts = TRUE}. Here, shared specificity is defined as a common ordering of
#' k-mer affinities across conditions. Therefore, differential specificity is tested by
#' assessing the statistical significance of deviations from a common trend
#' of k-mer affinities observed between conditions.
#'
#' For each condition, the common trend is estimated by fitting a LOESS curve
#' to the corresponding columns of the \code{"contrastAverage"} (x) and
#' \code{"contrastDifference"} (y) assays of the input SummarizedExperiment object.
#' Then, the difference between the observed differential affinity and the
#' LOESS-estimated specificity trend is taken as the estimate of differential
#' specificity.
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
#' SummarizedExperiment of k-mer differential specificity results with the following
#' assays:
#' 
#' \itemize{
#' \item \code{"contrastAverage"}: input k-mer average affinities.
#' \item \code{"contrastDifference"}: input k-mer differential affinities.
#' \item \code{"contrastVariance"}: input k-mer differential affinity variances.
#' \item \code{"contrastFit"}: estimated specificity trend.
#' \item \code{"contrastResidual"}: residual from estimated specificity trend (\code{contrastDifference - contrastFit}).
#' \item \code{"specificityZ"}: studentized differences (\code{contrastResidual / sqrt(contrastVariance)}).
#' \item \code{"specificityP"}: two-sided tail p-values for studentized differences.
#' \item \code{"specificityQ"}: FDR-controlling Benjamini-Hochberg adjusted p-values.
#' }
#'
#' @importFrom broom tidy
#' @importFrom limma loessFit
#' @importFrom stats p.adjust pnorm
#' @importFrom dplyr select group_by left_join ungroup mutate
#' @importFrom tidyr nest unnest
#' @export
#' @author Patrick Kimes
kmerTestSpecificity <- function(se, span = 0.05, ...) {

    stopifnot(is(se, "SummarizedExperiment"))
    if (!all(c("contrastAverage", "contrastDifference", "contrastVariance") %in%
             assayNames(se))) {
        stop("Input SummarizedExperiment is missing k-mer contrast estimates.\n",
             "SummarizedExperiment should be created by calling kmerFit(..) with 'contrasts=TRUE'.")
    }

    kmers <- rowData(se)$seq

    ## gather data
    cmean <- broom::tidy(se, "contrastAverage", long = TRUE)
    cmean <- dplyr::select(cmean, seq, condition = cname, contrastAverage = value)
    cdiff <- broom::tidy(se, "contrastDifference", long = TRUE)
    cdiff <- dplyr::select(cdiff, seq, condition = cname, contrastDifference = value)
    cvar <- broom::tidy(se, "contrastVariance", long = TRUE)
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
    cdat <- dplyr::mutate(cdat, specificityP = 2 * stats::pnorm(-abs(specificityZ)))
    cdat <- dplyr::group_by(cdat, condition)
    cdat <- dplyr::mutate(cdat, specificityQ = stats::p.adjust(specificityP, method = "BH"))
    cdat <- dplyr::ungroup(cdat)

    ## tidy results to assays
    assaylist <- list(contrastAverage = .tidycol2mat(cdat, "contrastAverage", kmers, colnames(se)),
                      contrastDifference = .tidycol2mat(cdat, "contrastDifference", kmers, colnames(se)),
                      contrastVariance = .tidycol2mat(cdat, "contrastVariance", kmers, colnames(se)),
                      contrastFit = .tidycol2mat(cdat, "contrastFit", kmers, colnames(se)),
                      contrastResidual = .tidycol2mat(cdat, "contrastResidual", kmers, colnames(se)),
                      specificityZ = .tidycol2mat(cdat, "specificityZ", kmers, colnames(se)),
                      specificityP = .tidycol2mat(cdat, "specificityP", kmers, colnames(se)),
                      specificityQ = .tidycol2mat(cdat, "specificityQ", kmers, colnames(se)))

    rdat <- dplyr::select(cdat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]
    
    SummarizedExperiment(assays = assaylist, rowData = rdat)
}
