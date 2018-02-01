#' RMA Normalization of Samples
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' performs RMA-style normalization on each sample (currently, only quantile
#' normalization is supported) and returns thesame SummarizedExperiment object
#' with the normalized intensities.
#'
#' @param se SummarizedExperiment object containing PBM intensity data
#'
#' @return
#' SummarizedExperiment object with normalized intensities.
#'
#' @import SummarizedExperiment
#' @importFrom preprocessCore normalize.quantiles
#' @export
#' @author Patrick Kimes
rmaNormalize <- function(se) {
    ## perform quantile normalization
    new_assay <- preprocessCore::normalize.quantiles(as.matrix(assay(se, "gpr")))
    new_assay <- DataFrame(new_assay)
    names(new_assay) <- rownames(colData(se))
    
    ## modify input SummarizedExperiment
    assay(se, "gpr") <- new_assay

    ## add step to list
    if (! "steps" %in% names(metadata(se))) {
        metadata(se)$steps <- list()
    }
    metadata(se)$steps <- c(metadata(se)$steps, "rma normalization")

    return(se)
}
