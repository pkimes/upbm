#' RMA Normalization of Samples
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' performs RMA-style normalization on each sample (currently, only quantile
#' normalization is supported) and returns thesame SummarizedExperiment object
#' with the normalized intensities.
#'
#' @param se SummarizedExperiment object containing PBM intensity data
#' @param assay_name string name of the assay to adjust. (default = "fore")
#' @param .force logical whether to run normalization even if data
#'        has already been normalized within arrays. (default = FALSE)
#'
#' @return
#' SummarizedExperiment object with normalized intensities.
#'
#' @import SummarizedExperiment
#' @importFrom preprocessCore normalize.quantiles
#' @export
#' @author Patrick Kimes
rmaNormalize <- function(se, assay_name = "fore", .force = FALSE) {

    ## check if already normalized
    if (!.force) {
        stopifnot(is.null(metadata(se)$betweenArrayNormalization))
    }

    ## perform quantile normalization
    new_assay <- preprocessCore::normalize.quantiles(as.matrix(assay(se, assay_name)))
    new_assay <- DataFrame(new_assay)
    names(new_assay) <- rownames(colData(se))
    
    ## modify input SummarizedExperiment
    assay(se, assay_name) <- new_assay

    ## add step to metadata
    method_str <- paste("preprocessCore::normalize.quantiles ->", assay_name)
    metadata(se)$steps <- c(metadata(se)$steps, method_str)
    if (.force) {
        metadata(se)$betweenArrayNormalization <-
                       c(metadata(se)$betweenArrayNormalization, method_str)
    } else {
        metadata(se)$betweenArrayNormalization <- method_str
    }
    
    return(se)
}
