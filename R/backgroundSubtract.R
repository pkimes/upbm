#' @title Background subtract intensities
#'
#' @description
#' Given a PBMExperiment with foreground and background intensities,
#' this function performs simple background subtraction and returns the same
#' PBMExperiment object with the background intensities subtracted from the
#' foreground intensities.
#'
#' @param pe a PBMExperiment object containing PBM intensity data.
#' @param assay a numeric index or string specifying the foreground assay.
#'        (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param assayb a numeric index or string specifying the background assay.
#'        (default = \code{SummarizedExperiment::assayNames(pe)[2]})
#' @param keepb a logical value whether to keep the background assay
#'        after subtraction. (default = TRUE)
#' @param nonnegative a logical value whether to restrict intensities to
#'        non-negative values by setting negative values to NA. (default = TRUE)
#' @param verbose a logical value whether to print verbose output
#'        during analysis. (default = FALSE)
#' 
#' @return
#' PBMExperiment object with background subtracted intensities replacing the original
#' \code{assay} assay values.
#' 
#' @importFrom SummarizedExperiment assayNames assay
#' @export 
#' @author Patrick Kimes
backgroundSubtract <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                               assayb = SummarizedExperiment::assayNames(pe)[2],
                               keepb = TRUE, nonnegative = TRUE, verbose = FALSE) {
    stopifnot(is(pe, "PBMExperiment"))
    stopifnot(assay %in% SummarizedExperiment::assayNames(pe))
    stopifnot(assayb %in% SummarizedExperiment::assayNames(pe))

    if (verbose) {
        cat("|| upbm::backgroundSubstract \n")
        cat("|| - Starting background subtraction for", ncol(pe), "PBM scans.\n")
        cat("|| - Using foreground assay:", assay, "\n")
        cat("|| - Using background assay:", assayb, "\n")
    }
    
    ## perform quick subtraction
    new_assay <- as.matrix(SummarizedExperiment::assay(pe, assay))
    bkgd_assay <- as.matrix(SummarizedExperiment::assay(pe, assayb))
    new_assay <- new_assay - bkgd_assay

    if (nonnegative) {
        if (verbose) {
            cat("|| - Replacing all negative values after subtraction with NA (nonnegative = TRUE).\n")
        }
        new_assay[new_assay < 0] <- NA
    } else if (verbose) {
        cat("|| - Keeping all negative values after subtraction (nonnegative = FALSE).\n")
    }
    new_assay <- DataFrame(new_assay)
    names(new_assay) <- rownames(colData(pe))
    
    if (verbose) {
        cat("|| - Replacing values in", assay, "with subtracted intensities.\n")
    }
    ## modify input PBMExperiment
    SummarizedExperiment::assay(pe, assay) <- new_assay

    ## drop background assay unless dropping is specified
    if (!keepb) {
        if (verbose) {
            cat("|| - Dropping background assay (keepb = FALSE).\n")
        }
        SummarizedExperiment::assay(pe, assayb) <- NULL
    } else if (verbose) {
        cat("|| - Keeping background assay (keepb = TRUE).\n")
    }
    
    return(pe)
}
