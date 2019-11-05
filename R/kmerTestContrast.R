#' @title Test for k-mer differential affinities
#'
#' @description
#' Using estimated k-mer level affinities returned by \code{kmerFit}, this
#' function tests for differential affinities across conditions for each k-mer
#' separately. The call to \code{kmerFit} must have been made with
#' \code{contrasts = TRUE}.
#'
#' While a column for the baseline condition of the contrasts is included in the
#' returned SummarizedExperiment object to match the dimensions of the input
#' object, all values for this column are set to NA. 
#' 
#' @description
#' After estimating k-mer level affinities using probe-set aggregate, this
#' function tests for differential affinity between conditions. 
#'
#' @param se SummarizedExperiment of k-mer results from \code{\link{kmerFit}} with
#'        \code{contrasts = TRUE}.
#' 
#' @return
#' SummarizedExperiment of k-mer differential affinity results with the following
#' assays:
#' 
#' \itemize{
#' \item \code{"contrastAverage"}: input k-mer average affinities.
#' \item \code{"contrastDifference"}: input k-mer differential affinities.
#' \item \code{"contrastVariance"}: input k-mer differential affinity variances.
#' \item \code{"contrastZ"}: studentized differences (\code{contrastDifference / sqrt(contrastVariance)}).
#' \item \code{"contrastP"}: two-sided tail p-values for studentized differences.
#' \item \code{"contrastQ"}: FDR-controlling Benjamini-Hochberg adjusted p-values.
#' }
#'
#' @importFrom broom tidy
#' @importFrom stats p.adjust pnorm
#' @importFrom dplyr select group_by left_join ungroup mutate
#' @export
#' @author Patrick Kimes
kmerTestContrast <- function(se) {

    stopifnot(is(se, "SummarizedExperiment"))
    if (!all(c("contrastAverage", "contrastDifference", "contrastVariance") %in%
             assayNames(se))) {
        stop("Input SummarizedExperiment is missing k-mer contrast estimates.\n",
             "SummarizedExperiment should be created by calling kmerFit(..) with 'contrasts=TRUE'.")
    }
    
    kmers <- rowData(se)$seq

    ## gather data
    cdiff <- broom::tidy(se, "contrastDifference", long = TRUE)
    cdiff <- dplyr::select(cdiff, seq, condition = cname, contrastDifference = value)
    cvar <- broom::tidy(se, "contrastVariance", long = TRUE)
    cvar <- dplyr::select(cvar, seq, condition = cname, contrastVariance = value)
    
    cdat <- dplyr::left_join(cdiff, cvar, by = c("condition", "seq"))

    ## compute z-scores, p-values, and adjusted p-values
    cdat <- dplyr::mutate(cdat, contrastZ = contrastDifference / sqrt(contrastVariance))
    cdat <- dplyr::mutate(cdat, contrastP = 2 * stats::pnorm(-abs(contrastZ)))
    cdat <- dplyr::group_by(cdat, condition)
    cdat <- dplyr::mutate(cdat, contrastQ = stats::p.adjust(contrastP, method = "BH"))
    cdat <- dplyr::ungroup(cdat)

    ## tidy results to assays
    assaylist <- list(contrastAverage = assay(se, "contrastAverage"),
                      contrastDifference = .tidycol2mat(cdat, "contrastDifference", kmers, colnames(se)),
                      contrastVariance = .tidycol2mat(cdat, "contrastVariance", kmers, colnames(se)),
                      contrastZ = .tidycol2mat(cdat, "contrastZ", kmers, colnames(se)),
                      contrastP = .tidycol2mat(cdat, "contrastP", kmers, colnames(se)),
                      contrastQ = .tidycol2mat(cdat, "contrastQ", kmers, colnames(se)))

    rdat <- dplyr::select(cdat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]
    
    SummarizedExperiment(assays = assaylist, rowData = rdat)
}
