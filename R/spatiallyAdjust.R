#' Spatial Adjustment of Samples
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' computes the local intensity bias for each probe, and returns the
#' same SummarizedExperiment object with the bias subtracted from each probe.
#' Optionally, the per-probe bias can also be returned as an additional array
#' in the bias corrected SummarizedExperiment object. The spatial bias at each
#' probe is defined as the difference between the median intensity in a `k` by `k`
#' square region surrounding the probe, and the median intensity of all probes
#' across the array. This approach is taken directly from the original PBM
#' analysis pipeline described in Berger and Bulyk (Nature Protocols, 2008).
#' The size of the local region can be specified by the user, with a default
#' size of 15 x 15 (as used in the aforementioned publication). 
#' 
#' @param se SummarizedExperiment object containing PBM intensity data
#' @param k odd integer specifying size of region to use to for computing
#'        local bias (default = 15)
#' @param returnBias logical whether to include the spatial bias as an
#'        additional 'assay' (called 'spatialbias') in the returned
#'        SummarizedExperiment object.
#' 
#' @return
#' SummarizedExperiment object with spatially adjusted intensities.
#' Row and column metadata are copied from the original SummarizedExperiment
#' object.
#'
#' @import SummarizedExperiment
#' @importFrom tibble as_tibble
#' @importFrom tidyr gather spread unnest
#' @importFrom dplyr select mutate select_ do group_by left_join
#' @export
#' @author Patrick Kimes
spatiallyAdjust <- function(se, k = 15, returnBias = TRUE) {
    if (k %% 2 == 0) {
        stop("Local window size, k, must be an odd value.")
    }
    stopifnot(is(se, "SummarizedExperiment")) 
    if (! all(c("Row", "Column") %in% names(rowData(se)))) {
        stop("Specified SummarizedExperiment does not have necessary ",
             "'Row' and 'Column' information in rowData to perform ",
             "spatial adjustment")
    }
    stopifnot("gpr" %in% assayNames(se))

    ## extract intensities for easier manipulation
    intensity <- assay(se, "gpr")
   
    ## if ID column present, only use de Bruijn probes (mask others)
    if ("ID" %in% names(rowData(se))) {
        intensity[grepl("^dBr_", rowData(se)$ID), ] <- NA_real_
    }

    ## add row/column indicies
    intensity <- cbind(rowData(se)[, c("Row", "Column")], intensity)

    ## compute spatial adjustment 
    intensity <- tibble::as_tibble(intensity)
    intensity <- tidyr::gather(intensity, sample, value, -Row, -Column)
    intensity <- dplyr::group_by(intensity, sample)
    intensity <- dplyr::do(intensity, spatialmedian = .wrapSA(., k))
    intensity <- tidyr::unnest(intensity)

    ## spread back so samples are in separate columns
    intensity <- tidyr::spread(intensity, -Row, -Column)

    ## join to original rowData to get proper row orders 
    intensity <- dplyr::left_join(rowData(se), intensity, by = c("Row", "Column"))
    intensity <- dplyr::select_(intensity, paste0("-", names(rowData(se))))

    ## construct new SummarizedExperiment from input SummarizedExperiment
    new_se <- se
    
    ## add new assays
    if (returnBias) {
        new_assays <- list(gpr = , spatialbias = )
    } else {
        new_assays <- list(gpr = )
    }
    assays(new_se) <- new_assays
    
    ## add step to list
    if (! "steps" %in% names(metadata(new_se))) {
        metadata(new_se)$steps <- list()
    }
    metadata(new_se)$steps <- c(metadata(new_se)$steps, "spatial adjustment")
    
    return(new_se) 
}

## Helper function which takes a data.frame with 'value', 'Column', 'Row'
## columns and returns the spatial median within a 'k' by 'k' region
## surrounding each value, and returns the values as a similar data.frame
## with columns, 'value', 'Column', 'Row'.
##
## @param x data.frame with columns 'value', 'Column', 'Row'
## @param k size of local region for computing medians
##
## @return
## a data.frame with columns 'value', 'Column', 'Row'
## 
## @author Patrick Kimes
.wrapSA <- function(x, k) {
    y <- dplyr::select(x, value, Column, Row)
    y <- tidyr::spread(y, Column, value)
    y <- dplyr::select(y, -Row)
    y <- blockmedian(as.matrix(y), k, center = TRUE)
    y <- tibble::as_tibble(y) 
    y <- dplyr::mutate(y, Row = row_number())
    y <- tidyr::gather(y, Column, value, -Row)
    y <- dplyr::mutate(y, Column = as.integer(gsub("V", "", Column)))
    y
}
