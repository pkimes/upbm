#' @title PBM Preprocess
#'
#' @description 
#' This is the main wrapper function for the core PBM data pre-processing steps
#' implemented in the upbm package.
#' 
#' @return
#' PBMExperiment containing filtered and normalized Alexa488 probe intensities.
#'
#' @author Patrick Kimes
upbmPreprocess <- function(pe, assayfore = SummarizedExperiment::assayNames(pe)[1],
                           
                           ## background subtraction
                           assayback = SummarizedExperiment::assayNames(pe)[2],
                           keepb = TRUE,
                           nonnegative = TRUE,
                           
                           ## cy3 fitting
                           cy3pe = NULL,
                           refpe = NULL,
                           assaycy3 = assayfore,
                           cy3method = "empirical", 
                           useMean = TRUE,
                           standardize = TRUE,
                           threshold = 1/2,
                           refit = TRUE,
                           
                           ## cy3 normalization
                           match_by,
                           filter = TRUE,
                           scale = TRUE,

                           ## spatial adjustment
                           k = 15L,
                           returnBias = TRUE,

                           ## within-replicate normalization
                           method = c("tmm", "quantile"),
                           q = 0.6,
                           qlower = 0,
                           qdiff = 0.2,
                           group = "id",
                           stratify = "condition",
                           baseline = NULL,

                           ## cross-replicate normalization
                           pairwise = FALSE,
                           onlybaseline = FALSE, 
                           
                           verbose = 1L) {
    
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
    
    ## background subtract
    if (verbose > 0L) {
        cat("||\n") 
        cat("|| Background intensity subtraction ... \n")
    }
    npe <- backgroundSubtract(pe, assay = assayfore, assayb = assayback,
                              keepb = keepb, verbose = (verbose > 1L),
                              nonnegative = nonnegative)
    if (verbose > 0L) {
        cat("|| ... finished!\n")
    }

    ## cy3 normalize
    if (is.null(cy3pe)) {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Skipping Cy3 normalization. \n")
        }
    } else {
        if (verbose > 0L) {
            cat("||\n") 
            cat("|| Cy3 normalization using [cy3method =", paste0("\"", cy3method, "\"]"), "... \n")
        }
        if (cy3method == "empirical") {
            cy3fit <- cy3FitEmpirical(cy3pe, refpe, assay = assaycy3,
                                      useMean = useMean, standardize = standardize,
                                      threshold = threshold, verbose = (verbose > 1L))
        } else {
            cy3fit <- cy3FitModel(cy3pe, assay = assaycy3, refit = refit,
                                  threshold = threshold, verbose = (verbose > 1L))
        }
        npe <- cy3Normalize(npe, cy3fit, assay = assayfore, match_by = match_by,
                            filter = filter, scale = scale, verbose = (verbose > 1L))
    if (verbose > 0L) {
        cat("|| ... finished!\n")
    }
    }
    
    ## spatial adjustment
    if (verbose > 0L) {
        cat("||\n") 
        cat("|| Spatial adjustment ... \n")
    }
    npe <- spatiallyAdjust(npe, k = k, returnBias = returnBias, verbose = (verbose > 1L))

    ## normalize within replicates
    if (verbose > 0L) {
        cat("|| ... finished!\n") 
        cat("||\n") 
        cat("|| Within-replicate normalization ... \n")
    }
    npe <- normalizeWithinReplicates(npe,
                                     method = method, q = q, qlower = qdiff,
                                     group = group, stratify = stratify, baseline = baseline,
                                     verbose = (verbose > 1L))
    
    ## normalize across replicates
    if (verbose > 0L) {
        cat("|| ... finished!\n") 
        cat("||\n") 
        cat("|| Cross-replicate normalization ... \n")
    }
    npe <- normalizeAcrossReplicates(npe, group = group, stratify = stratify, baseline = baseline,
                                     pairwise = pairwise, onlybaseline = onlybaseline,
                                     verbose = (verbose > 1L))

    ## return results
    if (verbose > 0L) {
        cat("|| ... finished!\n") 
        cat("||\n") 
        cat("|| Finished PBM data preprocessing.\n")
        cat("|| Returning PBMExperiment with", nrow(npe), "rows and", ncol(npe), "columns.\n")
    }
    return(npe)
}
