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
#' @import preprocessCore
#' @export 
#' @author Patrick Kimes
rmaBackgroundCorrect <- function(se) {
    ## perform RMA background correction (normal, exponential mixture)
    new_assay <- rma.background.correct(as.matrix(assay(se, "gpr")))
    
    ## construct new SummarizedExperiment from input SummarizedExperiment
    new_se <- se

    ## replace assays
    assays(new_se) <- list(gpr = DataFrame(new_assay))
    
    ## add step to list
    if (! "steps" %in% names(metadata(new_se))) {
        metadata(new_se)$steps <- list()
    }
    metadata(new_se)$steps <- c(metadata(new_se)$steps, "background correction")

    return(new_se)
}
