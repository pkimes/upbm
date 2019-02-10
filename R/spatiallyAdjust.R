#' @title Perform spatial adjustment
#'
#' @description
#' Given a PBMExperiment containing probe intensity data with array coordinates
#' specified in rowData as \code{"Column"} and \code{"Row"}, this function
#' computes the local intensity bias for each probe, and returns the
#' same PBMExperiment object with an additional array containing spatially
#' corrected (scaled) probe intensities. Optionally, the per-probe bias is also
#' returned as an array in the bias corrected PBMExperiment object.
#'
#' The spatial bias at each probe is computed as the ratio between the median
#' intensity in a \code{k} by \code{k}
#' region surrounding the probe and the median intensity of all probes
#' across the array. This approach is taken directly from the Universal PBM
#' Analysis Suite described in Berger and Bulyk (Nature Protocols, 2008).
#' 
#' @param pe a PBMExperiment object containing PBM intensity data.
#' @param assay a string name of the assay to adjust. (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param k an integer specifying the size of the region to use to for computing
#'        local bias. Must be odd. (default = 15L)
#' @param returnBias a logical whether to include the spatial bias as an
#'        additional 'assay' (called 'spatialbias') in the returned
#'        PBMExperiment object. (default = TRUE)
#' @param verbose a logical value whether to print verbose output
#'        during analysis. (default = FALSE)
#' 
#' @return
#' Original PBMExperiment object with assay containing spatially adjusted intensities
#' (\code{"normalized"}) and a new column added to the colData,
#' \code{"spatialMedian"}, containing the global median intensity of the original
#' probe-level data used to compute spatial bias for each sample.
#' If specified, the estimated spatial bias will also be included in an additional
#' assay (\code{"spatialbias"}). If assays with the same name are already included
#' in the object, they will be overwritten.
#' 
#' @references
#' \itemize{
#' \item Berger, M. F., & Bulyk, M. L. (2009). Universal protein-binding microarrays for the comprehensive characterization of the DNA-binding specificities of transcription factors. Nature Protocols, 4(3), 393-411.
#' }
#'
#' @import SummarizedExperiment
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather spread unnest
#' @importFrom dplyr select mutate select_ do group_by left_join row_number
#' @importFrom S4Vectors SimpleList
#' @export
#' @author Patrick Kimes
spatiallyAdjust <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                            k = 15L, returnBias = TRUE, verbose = FALSE) {
    stopifnot(is(pe, "PBMExperiment")) 

    ## don't show progress when running `do` here
    base_showprogress <- getOption("dplyr.show_progress")
    options(dplyr.show_progress = FALSE)
    on.exit(options(dplyr.show_progress = base_showprogress))
    
    if (k %% 2 == 0) {
        stop("Local window size, k, must be an odd value.")
    }
    if (! all(c("Row", "Column") %in% names(rowData(pe)))) {
        stop("Specified SummarizedExperiment does not have necessary ",
             "'Row' and 'Column' information in rowData to perform ",
             "spatial adjustment")
    }
    stopifnot(assay %in% SummarizedExperiment::assayNames(pe))

    if (verbose) {
        cat("|| upbm::spatiallyAdjust \n")
        cat("|| - Starting spatial adjustment for", ncol(pe), "PBM scans.\n")
    }

    if (verbose) {
        cat("|| - Filtering probes according to", length(pe@probeFilter),
            "probeFilter rule(s).\n")
        ntotal <- nrow(pe)
    }

    ## filter using rules
    pe <- pbmFilterProbes(pe)

    if (verbose) {
        cat("|| - Data filtered from", ntotal, "probes to", nrow(pe), "probes.\n")
    }
    
    ## extract intensities for easier manipulation
    raw_intensity <- SummarizedExperiment::assay(pe, assay)
    raw_intensity <- as.data.frame(raw_intensity, optional = TRUE)
    raw_intensity <- tibble::as_tibble(raw_intensity)

    ## add row/column indices
    raw_intensity <- dplyr::mutate(raw_intensity,
                                   Row = rowData(pe)[, "Row"],
                                   Column = rowData(pe)[, "Column"])

    ## compute spatial adjustment 
    raw_intensity <- tidyr::gather(raw_intensity, sample, value, -Row, -Column)
    med_intensity <- dplyr::group_by(raw_intensity, sample)
    med_intensity <- dplyr::do(med_intensity, spatialmedian = .wrapSA(., k))
    med_intensity <- dplyr::ungroup(med_intensity)
    
    ## store global medians for future use
    raw_medians <- dplyr::mutate(med_intensity,
                                 spatialMedian = vapply(spatialmedian, `[[`, 2, FUN.VALUE = 1.1))
    raw_medians <- dplyr::select(raw_medians, sample, spatialMedian)
    
    ## finish unwrapping spatial adjustments
    med_intensity <- dplyr::mutate(med_intensity, spatialmedian = lapply(spatialmedian, `[[`, 1))
    med_intensity <- tidyr::unnest(med_intensity)
    
    ## join median deviations w/ original raw intensities, subtract
    sub_intensity <- dplyr::left_join(raw_intensity, med_intensity,
                                      by = c("Row", "Column", "sample"),
                                      suffix = c(".raw", ".med"))

    ## only replace NA w/ 1 in median intensities
    ## - if NA in raw want NA in output
    sub_intensity <- dplyr::mutate(sub_intensity,
                                   value.med = ifelse(is.na(value.med), 1, value.med))
    sub_intensity <- dplyr::mutate(sub_intensity, value = value.raw / value.med)
    sub_intensity <- dplyr::select(sub_intensity, Row, Column, sample, value)

    ## spread back so samples are in separate columns
    med_intensity <- tidyr::spread(med_intensity, sample, value)
    sub_intensity <- tidyr::spread(sub_intensity, sample, value)

    ## join to original rowData to get proper row orders 
    med_intensity <- dplyr::left_join(as.data.frame(rowData(pe)), med_intensity,
                                      by = c("Row", "Column"))
    sub_intensity <- dplyr::left_join(as.data.frame(rowData(pe)), sub_intensity,
                                      by = c("Row", "Column"))

    ## only keep data columns
    med_intensity <- dplyr::select_(med_intensity, .dots = paste0("-", names(rowData(pe))))
    sub_intensity <- dplyr::select_(sub_intensity, .dots = paste0("-", names(rowData(pe))))

    ## reorder columns to match original SummarizedExperiment
    med_intensity <- DataFrame(med_intensity, check.names = FALSE)
    med_intensity <- med_intensity[, rownames(colData(pe)), drop = FALSE]
    sub_intensity <- DataFrame(sub_intensity, check.names = FALSE)
    sub_intensity <- sub_intensity[, rownames(colData(pe)), drop = FALSE]

    ## add to input PBMExperiment
    if ("normalized" %in% assayNames(pe)) {
        SummarizedExperiment::assay(pe, "normalized") <- NULL
        if (verbose) {
            cat("|| - Original PBMExperiment object contains \"normalized\" assay.\n")
            cat("|| - Existing \"normalized\" assay will be overwritten.\n")
        }
    }
    SummarizedExperiment::assays(pe) <- c(S4Vectors::SimpleList(normalized = sub_intensity),
                                          SummarizedExperiment::assays(pe))
    if (returnBias) {
        SummarizedExperiment::assay(pe, "spatialbias") <- med_intensity
    }

    if (verbose) {
        cat("|| - Finished spatial adjustment.\n")
    }
    return(pe) 
}

## Helper function which takes a data.frame with 'value', 'Column', 'Row'
## columns and returns the spatial median within a 'k' by 'k' region
## surrounding each value, and returns the values as a similar data.frame
## with columns, 'value', 'Column', 'Row'.
##
## @param x data.frame with columns 'value', 'Column', 'Row'.
## @param k size of local region for computing medians.
## 
## @return
## a list containing a data.frame with columns 'value', 'Column', 'Row',
## and a numeric value corresponding to the global median
## 
## @author Patrick Kimes
.wrapSA <- function(x, k) {
    y <- dplyr::select(x, value, Column, Row)
    y <- tidyr::spread(y, Column, value)
    y <- dplyr::select(y, -Row)
    y <- blockmedian(as.matrix(y), k)
    yg <- y$global
    y <- y$local / y$global
    y <- tibble::as_tibble(y)
    y <- dplyr::mutate(y, Row = dplyr::row_number())
    y <- tidyr::gather(y, Column, value, -Row)
    y <- dplyr::mutate(y, Column = as.integer(gsub("V", "", Column)))
    return(list(y, yg))
}
