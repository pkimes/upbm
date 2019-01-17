#' Test k-mer Contrasts
#'
#' @description
#' After estimating k-mer level affinities using probe-set aggregate, this
#' function tests for differential affinity between conditions. 
#'
#' @param se SummarizedExperiment of k-mer results from \code{kmerFit}.
#' 
#' @return
#' SummarizedExperiment of significance results for testing k-mer affinity
#' differences between conditions.
#'
#' @importFrom broom tidy
#' @importFrom stats p.adjust
#' @importFrom dplyr select group_by left_join ungroup mutate
#' @importFrom tidyr nest unnest
#' @export
#' @author Patrick Kimes
kmerTestContrast <- function(se) {

    kmers <- rowData(se)$seq

    ## gather data
    cdiff <- broom::tidy(se, "contrastDifference", long = TRUE, .filter = 0L)
    cdiff <- dplyr::select(cdiff, seq, condition = cname, contrastDifference = value)
    cvar <- broom::tidy(se, "contrastVariance", long = TRUE, .filter = 0L)
    cvar <- dplyr::select(cvar, seq, condition = cname, contrastVariance = value)
    
    cdat <- dplyr::left_join(cdiff, cvar, by = c("condition", "seq"))

    ## compute z-scores, p-values, and adjusted p-values
    cdat <- dplyr::mutate(cdat, contrastZ = contrastDifference / sqrt(contrastVariance))
    cdat <- dplyr::mutate(cdat, contrastP = 2 * pnorm(-abs(contrastZ)))
    cdat <- dplyr::group_by(cdat, condition)
    cdat <- dplyr::mutate(cdat, contrastQ = p.adjust(contrastP, method = "BH"))
    cdat <- dplyr::ungroup(cdat)

    ## tidy results to assays
    assaylist <- list(contrastAverage = assay(se, "contrastAverage"),
                      contrastDifference = assay(se, "contrastDifference"),
                      contrastZ = upbm:::.tidycol2mat(cdat, "contrastZ", kmers, colnames(se)),
                      contrastP = upbm:::.tidycol2mat(cdat, "contrastP", kmers, colnames(se)),
                      contrastQ = upbm:::.tidycol2mat(cdat, "contrastQ", kmers, colnames(se)))

    rdat <- dplyr::select(cdat, seq)
    rdat <- rdat[match(kmers, rdat$seq), ]
    
    SummarizedExperiment(assays = assaylist, rowData = rdat)
}
