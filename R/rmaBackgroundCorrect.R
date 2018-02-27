#' RMA Background Correction of Samples
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' performs RMA-style background correction on each sample, and returns the
#' same SummarizedExperiment object with the background corrected intensities.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data
#' @param assay_name string name of the assay to adjust. (default = "fore")
#' @param .force logical whether to run correction even if data
#'        has already been corrected for background. (default = FALSE)
#'
#' @return
#' SummarizedExperiment object with background corrected intensities.
#' 
#' @import SummarizedExperiment
#' @importFrom preprocessCore rma.background.correct
#' @export 
#' @author Patrick Kimes
rmaBackgroundCorrect <- function(se, assay_name = "fore", .force = FALSE) {

    new_assay <- as.matrix(assay(se, assay_name))
    
    ## check if already corrected
    if (!.force) {
        stopifnot(is.null(metadata(se)$backgroundCorrection))
    }

    ## perform RMA background correction (normal, exponential mixture)
    new_assay <- preprocessCore::rma.background.correct(new_assay)
    new_assay <- DataFrame(new_assay)
    names(new_assay) <- rownames(colData(se))
    
    ## modify input SummarizedExperiment
    assay(se, assay_name) <- new_assay
    
    ## add step to metadata
    method_str <- paste("preprocessCore::rma.background.correct ->", assay_name)
    metadata(se)$steps <- c(metadata(se)$steps, method_str)
    if (.force) {
        metadata(se)$backgroundCorrection <-
                       c(metadata(se)$backgroundCorrection, method_str)
    } else {
        metadata(se)$backgroundCorrection <- method_str
    }
    
    return(se)
}
