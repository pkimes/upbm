#' RMA Background Correction of Samples
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' performs RMA-style background correction on each sample, and returns the
#' same SummarizedExperiment object with the background corrected intensities.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data
#'
#' @return
#' SummarizedExperiment object with background corrected intensities.
#' 
#' @import SummarizedExperiment
#' @importFrom preprocessCore rma.background.correct
#' @export 
#' @author Patrick Kimes
rmaBackgroundCorrect <- function(se) {
    ## perform RMA background correction (normal, exponential mixture)
    new_assay <- preprocessCore::rma.background.correct(as.matrix(assay(se, "gpr")))
    new_assay <- DataFrame(new_assay)
    names(new_assay) <- rownames(colData(se))
    
    ## modify input SummarizedExperiment
    assay(se, "gpr") <- new_assay
    
    ## add step to list
    if (! "steps" %in% names(metadata(se))) {
        metadata(se)$steps <- list()
    }
    metadata(se)$steps <- c(metadata(se)$steps, "background correction")

    return(se)
}
