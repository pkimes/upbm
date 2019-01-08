#' SummarizedExperiment Wrapper for limma backgroundCorrect
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' performs background correction on each sample using limma's
#' \code{backgroundCorrect} function, and returns the same
#' SummarizedExperiment object with the background corrected intensities.
#'
#' @param se SummarizedExperiment object containing PBM intensity data.
#' @param assay string name of the assay to adjust. (default = \code{SummarizedExperiment::assayNames(se)[1]})
#' @param assayb string name of the background assay to use for
#'        adjustment. (default = NULL)
#' @param ... parameters to pass to \code{limma::backgroundCorrect.matrix}.
#'        See details below for more information on main parameters.
#' @param .force logical whether to run correction even if data
#'        has already been corrected for background. (default = FALSE)
#' @param .nonnegative logical whether to restrict intensities to non-negative
#'        values by setting negative values to NA. (default = TRUE)
#' 
#' @return
#' SummarizedExperiment object with background corrected intensities.
#'
#' @details
#' The most relevant parameters specified to \code{limma::backgroundCorrect.matrix}
#' are \code{method =} and \code{normexp.method}.
#' Options for \code{method} include: "none" and "normexp". The default for
#' single-channel arrays is "normexp". Options for \code{norm.exp} include:
#' "saddle", "mle", "rma" and "rma75". The default is "saddle".
#' For other parameter, see the \code{limma::backgroundCorrect.matrix}
#' documentation.
#' 
#' @import SummarizedExperiment
#' @importFrom limma backgroundCorrect.matrix
#' @export 
#' @author Patrick Kimes
limmaBackgroundCorrect <- function(se, assay = SummarizedExperiment::assayNames(se)[1],
                                   assayb = NULL, ..., .force = FALSE,
                                   .nonnegative = TRUE) {
    
    ## check if already corrected
    if (!.force) {
        stopifnot(is.null(metadata(se)$backgroundCorrection))
    }
    stopifnot(assay %in% SummarizedExperiment::assayNames(se))
    new_assay <- as.matrix(SummarizedExperiment::assay(se, assay))

    bkgd_assay <- NULL
    if (!is.null(assayb)) {
        stopifnot(assayb %in% SummarizedExperiment::assayNames(se))
        bkgd_assay <- as.matrix(SummarizedExperiment::assay(se, assayb))
    }
    
    ## perform RMA background correction (normal, exponential mixture)
    new_assay <- limma::backgroundCorrect.matrix(E = new_assay, Eb = bkgd_assay, ...)
    if (.nonnegative) {
        new_assay[new_assay < 0] <- NA
    }
    new_assay <- DataFrame(new_assay)
    names(new_assay) <- rownames(colData(se))
    
    ## modify input SummarizedExperiment
    SummarizedExperiment::assay(se, assay) <- new_assay
    
    ## add step to metadata
    method_str <- paste("limma::backgroundCorrect ->", assay)
    metadata(se)$steps <- c(metadata(se)$steps, method_str)
    if (.force) {
        metadata(se)$backgroundCorrection <-
                       c(metadata(se)$backgroundCorrection, method_str)
    } else {
        metadata(se)$backgroundCorrection <- method_str
    }
    
    return(se)
}
