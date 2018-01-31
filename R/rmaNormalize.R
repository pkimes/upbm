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
#' @import preprocessCore
#' @export
#' @author Patrick Kimes
rmaNormalize <- function(se) {
    ## perform quantile normalization
    new_assay <- preprocessCore::normalize.quantiles(as.matrix(assay(se, "gpr")))
    
    ## construct new SummarizedExperiment from input SummarizedExperiment
    new_se <- se

    ## replace assays
    assays(new_se) <- list(gpr = DataFrame(new_assay))
    
    ## add step to list
    if (! "steps" %in% names(metadata(new_se))) {
        metadata(new_se)$steps <- list()
    }
    metadata(new_se)$steps <- c(metadata(new_se)$steps, "rma normalization")

    return(new_se)
}
