#' @title Normalize within replicates
#'
#' @description
#' Universal PBM experiments are often performed with several conditions of interest,
#' e.g. various allelic variants of the same transcription factor, assayed on arrays of
#' the same plate, with few replicate plates (commonly 2 or 3). Within individual replicates (plates,
#' observed probe intensities can vary greatly across conditions for biologically
#' uninteresting reasons, such as concentration differences. To explicitly correct for these
#' differences, normalization is performed in two steps.
#'
#' First, normalization is performed within replicates (plates) with the assumption that
#' biologically uninteresting differences only affect probe intensities multiplicatively.
#' Normalization factors are estimated for each sample relative to a baseline condition on each
#' plate. The baseline should ideally be a replicate wild type or other natural reference condition
#' included on each plate. This function includes approaches for performing this step of
#' normalization.
#'
#' Second, normalization is performed across replicates (plates). More detail on this procedure 
#' can be found in the \code{\link{normalizeAcrossReplicates}} documentation.
#'
#' The approaches to normalization implemented in this function make a fundamental assumption
#' that lower-tail probe intensities are distributed similarly across the conditions being
#' normalized. This assumption is generally satisfied for allelic variants of the same transcription
#' factor or transcription factors with similar binding affinities. However, this assumption may
#' not always hold, e.g. if comparing proteins of completely different families. \emph{In these cases,
#' normalization should be performed with caution, and analyses and plots comparing the distributions
#' of lower-tail probe intensities should be explored.}
#'  
#' @param pe a SummarizedExperiment object containing GPR intensity information.
#' @param assay a string name of the assay to normalize.
#'        (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param method a string specifying the method to use for normalization. Must be one of
#'        \code{"tmm"} or \code{"quantile"}. Details on the methods are provided below.
#'        (default = \code{"tmm"})
#' @param q a percentile between 0 and 1 specifying either the upper quantile of probes to include 
#'        for normalization when \code{method = "tmm"} or the quantile to use for aligning samples
#'        when \code{method = "quantile"}. (default = 0.6)
#' @param qlower a percentile between 0 and 1-\code{q} specifying the lower quantile of
#'        probes to include for normalization when \code{method = "tmm"}. This parameter is
#'        ignored when \code{method = "quantile"}. (defalut = 0)
#' @param qdiff a percentile between 0 and 0.5 specifying the additional fraction of lower-tail
#'        probes to filter based on the deviation from the baseline condition when \code{method = "tmm"}.
#'        Probes with the \code{qdiff} smallest (most negative) and the \code{qdiff} largest (most positive)
#'        deviations from baseline condition will be filtered from normalization.
#'        This parameter is ignored when \code{method = "quantile"}. (default = 0.2)
#' @param group a character string specifying a column in \code{colData(pe)} to use for grouping replicates.
#'        (default = \code{"id"})
#' @param stratify a character string specifying a column in \code{colData(pe)} to use for determining
#'        the unique baseline scan within each \code{group}. (default = \code{"condition"})
#' @param baseline a character string specifying the baseline condition in the \code{stratify} column to normalize
#'        other conditions against within each \code{group}. If not specified and set to NULL, the baseline
#'        value is guessed by looking for ``ref" in any value of the \code{stratify} column. If multiple
#'        matching values are found, an error is thrown. If the baseline condition is missing from any
#'        \code{group}, an error is thrown. (default = NULL)
#' @param verbose a logical value whether to print verbose output during analysis. (default = FALSE)
#'
#' @details
#' The trimmed mean of M-values (\code{"tmm"}) method implemented in this function for cross-sample normalization
#' within replicates is based on the popular TMM method for RNA-seq data included
#' in the \code{edgeR} package. Very simply, a normalization factor is estimated as the trimmed mean
#' of probe-level log-scale differences between the baseline condition and sample using the lower
#' \code{[qlower, q]} percentile probes. Probes are ordered by the log-scale average intensity across
#' the baseline condition and sample. The trimmed mean is calculated excluding the top and bottom \code{qdiff}
#' probes.
#'
#' Unlike RNA-seq expression estimates, PBM data show near-constant variance in log-scale differences
#' as a function of the log-scale mean intensities. Therefore, a simplified variant of the original
#' TMM method is used, where precision weights are not introduced.
#'
#' The quantile-based (\code{"quantile"}) method \emph{should not be confused with what is commonly referred to
#' as ``quantile normalization."} Here, quantile-based normalization computes scaling factors across 
#' 
#' @return
#' Original PBMExperiment object with assay containing within-replicate normalized intensities
#' (\code{"normalized"}) and a new column added to the colData, \code{"withinRepScale"},
#' containing the inverse of the scaling factors used to normalize intensities.
#' If an assay with the same name is already included in the object, it will be overwritten.
#'
#' @seealso \code{\link{normalizeAcrossReplicates}}
#' @importFrom stats quantile
#' @importFrom dplyr as_tibble tibble mutate group_by filter do select ungroup left_join
#' @importFrom tidyr gather spread nest
#' @importFrom rlang enquo quo_name
#' @export
#' @author Dongyuan Song, Patrick Kimes
normalizeWithinReplicates <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                                      method = c("tmm", "quantile"), q = 0.6, qlower = 0, qdiff = 0.2,
                                      group = "id", stratify = "condition", baseline = NULL,
                                      verbose = FALSE) {
    stopifnot(is(pe, "PBMExperiment"))
    method <- match.arg(method)
    
    stopifnot(q >= 0, q < 1)
    stopifnot(qlower < 1 - q)
    stopifnot(qdiff >= 0, qdiff <= 0.5)
    stopifnot(assay %in% SummarizedExperiment::assayNames(pe))

    ## check normalization groups
    stopifnot(group %in% names(colData(pe)))
    
    if (verbose) {
        cat("|| upbm::normalizeWithinReplicates \n")
        cat("|| - Starting within-replicate normalization for", ncol(pe), "PBM scans.\n")
    }

    if (verbose) {
        cat("|| - Filtering probes according to", length(pe@probeFilter),
            "probeFilter rule(s).\n")
        ntotal <- nrow(pe)
    }

    ## filter probes - only for computing shift/scale factors (return original pe)
    fpe <- pbmFilterProbes(pe)
    rdat <- dplyr::as_tibble(as.data.frame(rowData(fpe)[, fpe@probeCols, drop = FALSE], optional = TRUE))

    if (verbose) {
        cat("|| - Data filtered from", ntotal, "probes to", nrow(fpe), "probes.\n")
    }

    ## check normalization stratification settings
    strats <- .pbmCheckStratify(fpe, stratify, baseline, group)
    coldat <- strats$coldat
    baseline <- strats$baseline

    if (verbose) {
        cat("|| - Performing within-replicate normalization with:\n")
        cat("||     - replicate group column:", group, "\n")
        cat("||     - condition column:", stratify, "\n")
        cat("||     - baseline condition:", baseline, "\n")
    }

    ## tidy up data for computing factors
    scale_assay <- as.data.frame(assay(fpe, assay), optional = TRUE)
    scale_assay <- dplyr::as_tibble(scale_assay)
    scale_assay <- dplyr::bind_cols(scale_assay, rdat)

    ## need to merge with columns from row data 
    scale_assay <- tidyr::gather(scale_assay, sample, value, -(!! fpe@probeCols))
    scale_assay <- dplyr::left_join(scale_assay, coldat, by = "sample")
    
    assay_fits <- dplyr::filter(scale_assay, !is.na(value), value > 0)

    if (method == "tmm") {
        if (verbose) {
            cat("|| - Performing TMM normalization with:\n")
            cat("||     - q:", q, "\n")
            cat("||     - qlower:", qlower, "\n")
            cat("||     - qdiff:", qdiff, "\n")
        }
        bl_assay <- dplyr::filter(assay_fits, Stratify == baseline)
        bl_assay <- dplyr::select_at(bl_assay, c("value", "Group", fpe@probeCols))
        assay_fits <- dplyr::left_join(assay_fits, bl_assay, by = c("Group", fpe@probeCols),
                                       suffix = c("", ".bl"))
        assay_fits <- dplyr::filter(assay_fits, !is.na(value.bl))
        assay_fits <- dplyr::mutate(assay_fits,
                                    M.value = (log2(value) - log2(value.bl)),
                                    A.value = (log2(value) + log2(value.bl))/2)  
        assay_fits <- tidyr::nest(assay_fits, -sample, -Stratify, -Group)

        assay_fits <- dplyr::mutate(assay_fits,
                                    withinRepScale = vapply(data, function(x) {
                                        z <- dplyr::filter(x, A.value < quantile(A.value, q))
                                        z <- dplyr::filter(z, A.value > quantile(A.value, qlower))
                                        z <- dplyr::summarise(z, est_scale = 2^mean(M.value, trim = qdiff, na.rm = TRUE))
                                        as.numeric(z$est_scale)
                                    }, numeric(1L)))
        assay_fits <- dplyr::select(assay_fits, -data)
        
    } else if (method == "quantile") {
        if (verbose) {
            cat("|| - Performing quantile-based normalization with:\n")
            cat("||     - q:", q, "\n")
        }
        assay_fits <- dplyr::group_by(scale_assay, sample, Stratify, Group)
        assay_fits <- dplyr::summarize(assay_fits, ul = quantile(value, probs = q, na.rm = TRUE))
        assay_fits <- dplyr::ungroup(assay_fits)

        assay_fits <- dplyr::group_by(assay_fits, Group)
        assay_fits <- dplyr::mutate(assay_fits, withinRepScale = ul / ul[Stratify == baseline])
        assay_fits <- dplyr::ungroup(assay_fits)
        assay_fits <- dplyr::select(assay_fits, -ul)

    } else {
        stop("Specified 'method=' parameter is invalid.")
    } 
    
    ## tidy up original data
    new_assay <- as.data.frame(assay(pe, assay), optional = TRUE)
    new_assay <- dplyr::as_tibble(new_assay)
    new_assay <- dplyr::mutate(new_assay,
                               Row = rowData(pe)[, "Row"],
                               Column = rowData(pe)[, "Column"])
    new_assay <- tidyr::gather(new_assay, sample, value, -Row, -Column)
    new_assay <- dplyr::left_join(new_assay, coldat, by = "sample")
    
    ## adjust to reference
    new_assay <- dplyr::left_join(new_assay, assay_fits, by = c("sample", "Stratify"))
    new_assay <- dplyr::mutate(new_assay, value = value / withinRepScale)
    
    ## return to square assay shape
    new_assay <- dplyr::select(new_assay, sample, value, Row, Column)
    new_assay <- tidyr::spread(new_assay, sample, value)
    
    ## match row order to rowData
    c_order <- paste(rowData(pe)$Row, rowData(pe)$Column, sep = "-")
    new_order <- match(c_order, paste(new_assay$Row, new_assay$Column, sep = "-"))
    stopifnot(!duplicated(new_order), length(new_order) == nrow(pe))
    new_assay <- new_assay[new_order, ]
    new_assay <- dplyr::select(new_assay, -Row, -Column)
    
    ## match column order to colData
    stopifnot(colnames(new_assay) %in% colnames(pe))
    new_assay <- new_assay[, colnames(pe)]
    new_assay <- DataFrame(new_assay, check.names = FALSE)
    
    ## add to input PBMExperiment
    if ("normalized" %in% assayNames(pe)) {
        SummarizedExperiment::assay(pe, "normalized") <- NULL
        if (verbose) {
            cat("|| - Original PBMExperiment object contains \"normalized\" assay.\n")
            cat("|| - Existing \"normalized\" assay will be overwritten.\n")
        }
    }
    SummarizedExperiment::assays(pe) <- c(S4Vectors::SimpleList(normalized = new_assay),
                                          SummarizedExperiment::assays(pe))
    
    if (verbose) {
        cat("|| - Adding within-replicate normalization factors to colData as \"withinRepScale\".\n")
    }

    ## add scaling factors to colData
    assay_fits <- dplyr::select(assay_fits, -Stratify, -Group)
    if ("fits" %in% names(assay_fits)) {
        assay_fits <- dplyr::select(assay_fits, -fits)
    }

    ## merge by row names
    coldat <- merge(colData(pe), data.frame(assay_fits, row.names = "sample"),
                    by = 0, all = TRUE)
    rownames(coldat) <- coldat$Row.names
    coldat$Row.names <- NULL
    
    ## match colData row order with SE col order
    coldat <- coldat[match(colnames(pe), rownames(coldat)), , drop = FALSE]
    
    stopifnot(all(rownames(coldat) == colnames(pe)))
    colData(pe) <- coldat
    
    if (verbose) {
        cat("|| - Finished within-replicate normalization.\n")
        cat("|| - Returning PBMExperiment with", nrow(pe), "rows and", ncol(pe), "columns.\n")
    }

    return(pe)
}
