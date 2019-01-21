#' Generate Cy3 Empirical Reference
#' 
#' Cy3 normalization is performed using an empirically-derived
#' reference distribution computed using a large collection of
#' samples. By default, samples are first scaled to have a common
#' median intensity, and the reference probe-level intensities are
#' computed on the log2 scale. The reference is calculated on the
#' log2 scale rather than the raw scale because the Cy3 scans will be
#' used for filtering and scaling, actions which require examining
#' fold changes and ratios rather than raw intensity differences.
#'
#' Scaling is performed by multiplying the raw intensities of each
#' sample by a scaling factor such that the median intensity of each
#' sample is equal to the median sample-median intensity of the original
#' intensities. This multiplicative scaling is equivalent to an additive
#' shift on the log2 scale.
#'
#' 
#' @param cy3se SummarizedExperiment object containing Cy3 scan data.
#' @param assay string name of assay to use. (default = \code{SummarizedExperiment::assayNames(cy3se)[1]})
#' @param offset numeric offset to add to intensities before log2 transforming
#'        to prevent errors with zero intensities. (default = 1L)
#' @param register logical whether to scale intensities across samples.
#'        (default = TRUE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to estimating affinities. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' 
#' @return
#' SummarizedExperiment object with Cy3 probe-level reference metrics. 
#' 
#' @export
#' @importFrom stats mad median sd
#' @importFrom dplyr group_by mutate ungroup summarize left_join select
#' @author Patrick Kimes
cy3GenerateRef <- function(cy3se, assay = SummarizedExperiment::assayNames(cy3se)[1],
                           offset = 1L, register = TRUE, .filter = 1L) {

    ## verify validity of cy3se
    stopifnot(is(cy3se, "SummarizedExperiment"))
    stopifnot("scan" %in% names(colData(cy3se)))
    if (any(colData(cy3se)$scan != "Cy3")) {
        stop("Not all scans in 'cy3se' are labelled as 'Cy3' in colData.\n",
             "Only Cy3 scans should be passed to this function.")
    }
    sum_cols <- paste0("probe_", c("median", "mean", "mad", "sd"))
    if (any(names(rowData(cy3se)) %in% sum_cols)) {
        stop("rowData of 'cy3se' has column names matching the following special names:\n",
             paste0(sum_cols, collapse = ", "), "\n",
             "Please rename these columns.")
    }
    
    ## tidy scan data for parsing
    cy3vals <- tidy.SummarizedExperiment(cy3se, assay, long = TRUE, .filter = .filter)
    
    ## scale scans
    if (register) {
        ## common multiplication factor to scale samples
        col_medians <- colMedians(as.matrix(SummarizedExperiment::assay(cy3se, assay)), na.rm = TRUE)
        sfactor <- median(col_medians, na.rm = TRUE)
        
        cy3vals <- dplyr::group_by(cy3vals, id, idx)
        cy3vals <- dplyr::mutate(cy3vals, value = value / median(value, na.rm = TRUE))
        cy3vals <- dplyr::ungroup(cy3vals)
        cy3vals <- dplyr::mutate(cy3vals, value = value * sfactor)
    } else {
        sfactor <- NA
    }
    
    ## calculate metrics on log2 scale
    cy3vals <- dplyr::mutate(cy3vals, value = log2(value + offset))

    ## compute probe-level summary metrics
    cy3vals <- dplyr::group_by(cy3vals, Column, Row)
    cy3vals <- dplyr::summarize(cy3vals,
                                probe_median = median(value, na.rm = TRUE),
                                probe_mean = mean(value, na.rm = TRUE),
                                probe_mad = mad(value, na.rm = TRUE),
                                probe_sd = sd(value, na.rm = TRUE))
    cy3vals <- dplyr::ungroup(cy3vals)

    ## create SummarizedExperiment object
    rowdat <- as.data.frame(rowData(cy3se), optional = TRUE)
    cy3vals <- dplyr::left_join(dplyr::as_tibble(rowdat), cy3vals, by = c("Column", "Row"))
    rowdat <- DataFrame(dplyr::select(cy3vals, -one_of(sum_cols)))
    refassay <- dplyr::select(cy3vals, one_of(sum_cols))
    refassay <- list(ref = DataFrame(refassay))

    ## capture call parameters, combine with formal args
    params <- match.call(expand.dots = FALSE)[-1]
    params <- replace(formals(), names(params), params)

    ## return results as a SummarizedExperiment object
    SummarizedExperiment(assays = refassay, rowData = rowdat,
                         metadata = list(sfactor = sfactor, params = params))
}
