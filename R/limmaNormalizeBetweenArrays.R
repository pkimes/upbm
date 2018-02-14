#' SummarizedExperiment Wrapper for limma NormalizeBetweenArrays
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' performs normalization across samples using limma's
#' \code{normalizeBetweenArrays} function, and returns the same
#' SummarizedExperiment object with the normalized intensities.
#'
#' @param se SummarizedExperiment object containing PBM intensity data
#' @param ... parameters to pass to \code{limma::normalizeBetweenArrays}.
#'        See details below for more information on main parameters.
#' @param .force logical whether to run normalization even if data
#'        has already been normalized within arrays. (default = FALSE)
#'
#' @return
#' SummarizedExperiment object with normalized intensities.
#'
#' @details
#' The most relevant parameter specified to \code{limma::normalizeBetweenArrays}
#' is \code{method =}. Options include: "none", "scale", "quantile", and
#' "cyclic loess." The default for single-channel arrays is "quantile".
#' For other parameter, see the \code{limma::normalizeBetweenArrays}
#' documentation.
#' 
#' @import SummarizedExperiment
#' @importFrom limma normalizeBetweenArrays
#' @export
#' @author Patrick Kimes
limmaNormalizeBetweenArrays <- function(se, ..., .force = FALSE) {

    ## check if already normalized
    if (!.force) {
        stopifnot(is.null(metadata(se)$betweenArrayNormalization))
    }
    
    new_assay <- as.matrix(assay(se, "gpr"))

    ## perform quantile normalization
    new_assay <- limma::normalizeBetweenArrays(object = new_assay, ...)
    new_assay <- DataFrame(new_assay)
    names(new_assay) <- rownames(colData(se))
    
    ## modify input SummarizedExperiment
    assay(se, "gpr") <- new_assay


    ## add step to metadata
    method_str <- "limma::normalizeBetweenArrays"
    metadata(se)$steps <- c(metadata(se)$steps, method_str)
    if (.force) {
        metadata(se)$betweenArrayNormalization <-
                       c(metadata(se)$betweenArrayNormalization, method_str)
    } else {
        metadata(se)$betweenArrayNormalization <- method_str
    }
    
    return(se)
}
