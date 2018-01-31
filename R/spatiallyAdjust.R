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
    raw_intensity <- assay(se, "gpr")
   
    ## if ID column present, only use de Bruijn probes (mask others)
    if ("ID" %in% names(rowData(se))) {
        raw_intensity[grepl("^dBr_", rowData(se)$ID), ] <- NA_real_
    }

    ## add row/column indicies
    raw_intensity <- cbind(rowData(se)[, c("Row", "Column")], raw_intensity)

    ## compute spatial adjustment 
    raw_intensity <- tibble::as_tibble(as.data.frame(raw_intensity))
    raw_intensity <- tidyr::gather(raw_intensity, sample, value, -Row, -Column)
    med_intensity <- dplyr::group_by(raw_intensity, sample)
    med_intensity <- dplyr::do(med_intensity, spatialmedian = .wrapSA(., k))
    med_intensity <- tidyr::unnest(med_intensity)

    ## join median deviations w/ original raw intensities, subtract
    sub_intensity <- dplyr::left_join(raw_intensity, med_intensity,
                                      by = c("Row", "Column", "sample"),
                                      suffix = c(".raw", ".med"))
    ## only replace NA w/ 0 in median intensities - if NA in raw (i.e. non-dBr probe) want NA in output
    sub_intensity <- dplyr::mutate(sub_intensity, value.med = ifelse(is.na(value.med), 0, value.med))
    sub_intensity <- dplyr::mutate(sub_intensity, value = value.raw - value.med)
    sub_intensity <- dplyr::select(sub_intensity, Row, Column, sample, value)
    
    ## spread back so samples are in separate columns
    med_intensity <- tidyr::spread(med_intensity, sample, value)
    sub_intensity <- tidyr::spread(sub_intensity, sample, value)

    ## join to original rowData to get proper row orders 
    med_intensity <- dplyr::left_join(as.data.frame(rowData(se)), med_intensity, by = c("Row", "Column"))
    sub_intensity <- dplyr::left_join(as.data.frame(rowData(se)), sub_intensity, by = c("Row", "Column"))

    ## only keep data columns
    med_intensity <- dplyr::select_(med_intensity, .dots = paste0("-", names(rowData(se))))
    sub_intensity <- dplyr::select_(sub_intensity, .dots = paste0("-", names(rowData(se))))

    ## reorder columns to match original SummarizedExperiment
    med_intensity <- DataFrame(med_intensity)
    med_intensity <- med_intensity[, rownames(colData(se))]
    sub_intensity <- DataFrame(sub_intensity)
    sub_intensity <- sub_intensity[, rownames(colData(se))]

    ## modify input SummarizedExperiment
    assay(se, "gpr") <- sub_intensity
    print("here")
    if (returnBias) {
        assay(se, "spatialbias") <- med_intensity
    }
    print("here")
    
    ## add step to list
    if (! "steps" %in% names(metadata(se))) {
        metadata(se)$steps <- list()
    }
    metadata(se)$steps <- c(metadata(se)$steps, "spatial adjustment")
    
    return(se) 
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
