#' @title PBM Preprocess
#'
#' @description 
#' This is the main wrapper function for the core PBM data pre-processing steps
#' implemented in the upbm package.
#'
#' @param pe a PBMExperiment object containing Alexa488 intensity data.
#' @param cy3pe an optional PBMExperiment object containing Cy3 intensity data. If
#'        corresponding Cy3 data is not available for the Alexa488 scans, this must be
#'        set to NULL.
#' @param cy3refpe an optional PBMExperiment object containing Cy3 reference data to
#'        use for Cy3 normalization with \code{\link{cy3FitEmpirical}} (recommended).
#'        To perform model-based normalization with \code{\link{cy3FitModel}},
#'        set to NULL. This option is ignored if \code{cy3pe = NULL}.
#' @param assay a numeric index or string specifying the foreground assay in
#'        both the Alexa488 and (optional) Cy3 PBMExperiment objects.
#'        Note that this assay is only used for the initial background subtraction and
#'        Cy3 normalization steps. All subsequent steps will use the default processed
#'        assays from the previous steps.
#'        (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param params a named list specifying any non-default parameters that should be used 
#'        for preprocessing the data. Only entries in the list exactly matching the name
#'        of a preprocessing step are considered valid (\code{"backgroundSubtract"},
#'        \code{"cy3Normalize"}, \code{"spatiallyAdjust"},
#'        \code{"normalizeWithinReplicates"}, or \code{"normalizeAcrossReplicates"}).
#'        For each step, non-default parameters must be specified as a list of named
#'        parameter values. (default = \code{list()}). See below for examples.
#' @param skip a string vector specifying any steps to skip during preprocessing. Strings
#'        must exactly match step names as with the \code{params} parameter.
#'        (default = \code{vector("character", 0)})
#' @param verbose an integer value specifying the level of verbosity. Must be one of
#'        0, 1, or 2. If 0, no progress messages are printed. If 1, progress messages are
#'        only provided to indicate when each preprocessing step is started and finished.
#'        If 2, in addition to the basic progress messages, \code{verbose = TRUE} is also specified
#'        for all preprocessing functions called internally. (default = \code{1L})
#' 
#' @examples
#' data(hoxc9alexa, package = "upbmData")
#' data(hoxc9cy3, package = "upbmData")
#' data(refcy3_8x60k_v1, package = "upbmAux")
#'
#' \dontrun{
#' upbmPreprocess(hoxc9alexa, hoxc9cy3, refcy3_8x60k_v1)
#'
#' upbmPreprocess(hoxc9alexa, NULL, NULL, params = list(backgroundSubtract = FALSE))
#'
#' upbmPreprocess(hoxc9alexa, hoxc9cy3, refcy3_8x60k_v1,
#'                params = list(backgroundSubtract = list(assayb = "back")),
#'                skip = "cy3Normalize")
#' }
#' 
#' @seealso \code{\link{backgroundSubtract}}, \code{\link{cy3FitEmpirical}}, \code{\link{cy3FitModel}}, \code{\link{cy3Normalize}}, \code{\link{spatiallyAdjust}}, \code{\link{normalizeWithinReplicates}}, \code{\link{normalizeAcrossReplicates}}
#' @return
#' PBMExperiment containing filtered and normalized Alexa488 probe intensities.
#'
#' @export
#' @author Patrick Kimes
upbmPreprocess <- function(pe, cy3pe, cy3refpe, assay = SummarizedExperiment::assayNames(pe)[1],
                           params = list(), skip = vector("character", 0), verbose = 1L) {

    stopifnot(is(pe, "PBMExperiment"))
    stopifnot(is.null(cy3pe) || is(cy3pe, "PBMExperiment"))
    stopifnot(is.null(cy3refpe) || is(cy3refpe, "PBMExperiment"))
    stopifnot(is(params, "list"))
    stopifnot(verbose %in% c(0L, 1L, 2L))
    stopifnot(assay %in% SummarizedExperiment::assayNames(pe))
    if (is(cy3pe, "PBMExperiment")) {
        stopifnot(assay %in% SummarizedExperiment::assayNames(cy3pe))
    }

    steps <- c("backgroundSubtract", "cy3Normalize", "spatiallyAdjust",
               "normalizeWithinReplicates", "normalizeAcrossReplicates")

    ## check validity of 'params='
    stopifnot(is(params, "list"))
    if (any(! names(params) %in% steps)) {
        stop("If the 'params=' parameter is specified, it must be a named list with names matching\n",
             "the following preprocessing steps:\n",
             paste(steps, collapse = ", "))
    }
    if (any(! vapply(params, is.list, logical(1L)))) {
        stop("If the 'params=' parameter is specified, each entry must be a list of non-default parameters.")
    }

    ## check validity of 'skip='
    stopifnot(is(skip, "vector") && is(skip, "character"))
    if (any(! names(params) %in% steps)) {
        stop("If the 'skip=' parameter is specified, it may only contain a subset of the following\n",
             " preprocessing steps:\n",
             paste(steps, collapse = ", "))
    }
    
    if (verbose > 0L) {
        cat("|| upbm::upbmPreprocess \n")
        cat("|| - Starting PBM data preprocessing.\n")
        cat("|| -", ncol(pe), "Alexa488 PBM scans.\n")
        if (is.null(cy3pe)) {
            cat("|| - no Cy3 PBM scans.\n")
        } else {
            cat("|| -", ncol(cy3pe), "Cy3 PBM scans.\n")
        }
    }

    npe <- pe
    
    ## background subtract
    if ("backgroundSubtract" %in% skip) {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Skipping background intensity subtraction. \n")
        }
    } else {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Background intensity subtraction ... \n")
        }

        p <- list(pe = npe, assay = assay, verbose = (verbose > 1L))
        if (!is.null(params[["backgroundSubtract"]])) {
            p <- replace(p, names(params[["backgroundSubtract"]]), params[["backgroundSubtract"]])
        }
        npe <- do.call(backgroundSubtract, p)
        if (!is.null(cy3pe)) {
            cy3pe <- do.call(backgroundSubtract, replace(p, "pe", list(cy3pe)))
        }

        if (verbose > 0L) {
            cat("|| ... finished!\n")
        }
    }
    
    ## cy3 normalize
    if ("cy3Normalize" %in% skip || is.null(cy3pe)) {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Skipping Cy3 normalization. \n")
        }
    } else {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Cy3 normalization using [",
                ifelse(is.null(cy3refpe), "cy3FitModel", "cy3FitEmpirical"), "] ... \n")
        }

        p <- list(pe = cy3pe, assay = assay, verbose = (verbose > 1L))
        
        if (is.null(cy3refpe)) {
            if (!is.null(params[["cy3FitModel"]])) {
                p <- replace(p, names(params[["cy3FitModel"]]), params[["cy3FitModel"]])
            }
            cy3fit <- do.call(cy3FitModel, p)
        } else {
            p$refpe <- cy3refpe
            if (!is.null(params[["cy3FitEmpirical"]])) {
                p <- replace(p, names(params[["cy3FitEmpirical"]]), params[["cy3FitEmpirical"]])
            }
            cy3fit <- do.call(cy3FitEmpirical, p)
        }

        p <- list(pe = npe, cy3pe = cy3fit, assay = assay, verbose = (verbose > 1L))
        if (!is.null(params[["cy3Normalize"]])) {
            p <- replace(p, names(params[["cy3Normalize"]]), params[["cy3Normalize"]])
        }
        cy3fit <- do.call(cy3Normalize, p)
        
        if (verbose > 0L) {
            cat("|| ... finished!\n")
        }
    }

    ## spatial adjustment
    if ("spatiallyAdjust" %in% skip) {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Skipping spatial adjustment. \n")
        }
    } else {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Spatial adjustment ... \n")
        }
        
        p <- list(pe = npe, verbose = (verbose > 1L))
        if (!is.null(params[["spatiallyAdjust"]])) {
            p <- replace(p, names(params[["spatiallyAdjust"]]), params[["spatiallyAdjust"]])
        }
        npe <- do.call(spatiallyAdjust, p)

        if (verbose > 0L) {
            cat("|| ... finished!\n")
        }
    }
    
    ## normalize within replicates
    if ("normalizeWithinReplicates" %in% skip) {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Skipping within-replicate normalization. \n")
        }
    } else {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Within-replicate normalization ... \n")
        }

        p <- list(pe = npe, verbose = (verbose > 1L))
        if (!is.null(params[["normalizeWithinReplicates"]])) {
            p <- replace(p, names(params[["normalizeWithinReplicates"]]),
                         params[["normalizeWithinReplicates"]])
        }
        npe <- do.call(normalizeWithinReplicates, p)

        if (verbose > 0L) {
            cat("|| ... finished!\n") 
        }
    }
    
    ## normalize across replicates
    if ("normalizeAcrossReplicates" %in% skip) {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Skipping cross-replicate normalization. \n")
        }
    } else {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Cross-replicate normalization ... \n")
        }
        p <- list(pe = npe, verbose = (verbose > 1L))
        if (!is.null(params[["normalizeAcrossReplicates"]])) {
            p <- replace(p, names(params[["normalizeAcrossReplicates"]]),
                         params[["normalizeAcrossReplicates"]])
        }
        npe <- do.call(normalizeAcrossReplicates, p)
        if (verbose > 0L) {
            cat("|| ... finished!\n") 
        }
    }
    
    ## return results
    if (verbose > 0L) {
        cat("||\n") 
        cat("|| Finished PBM data preprocessing.\n")
        cat("|| Returning PBMExperiment with", nrow(npe), "rows and", ncol(npe), "columns.\n")
    }
    return(npe)
}
