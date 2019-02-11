#' @title Generate Cy3 empirical reference
#'
#' @description
#' Given a collection of Cy3 scans, this function estimates an empirical
#' reference distribution of expected Cy3 intensities per probe. Scans
#' must all be for the same array design. For each probe, the mean,
#' median, MAD and standard deviation across scans is computed across all
#' Cy3 scans.
#'
#' The estimated reference distributions can be passed to \code{cy3FitEmpirical}
#' with Cy3 scan data to determine probe-level outliers and scaling factors
#' for the individual Cy3 scans relative to the reference. 
#'
#' @param cy3pe a PBMExperiment object containing Cy3 intensity data.
#' @param assay a numeric index or string name specifying the assay to use.
#'        (default = \code{SummarizedExperiment::assayNames(cy3pe)[1]})
#' @param offset a numeric offset to add to intensities before log2 transforming
#'        to prevent NAs with zero intensities. (default = 1L)
#' @param register logical whether to scale intensities across samples.
#'        (default = TRUE)
#' 
#' @return
#' SummarizedExperiment object with Cy3 probe-level reference metrics. 
#'
#' @details
#' By default, samples are first scaled to have a common
#' median intensity, and the reference probe-level intensities are
#' computed on the log2 scale. The reference is calculated on the
#' log2 scale rather than the raw scale because the Cy3 scans will be
#' used for filtering and scaling, actions which require examining
#' fold changes and ratios rather than raw intensity differences.
#'
#' Scaling is performed such that the median intensity of each
#' sample is equal to the median sample-median intensity of the original
#' intensities. This multiplicative scaling is equivalent to an additive
#' shift on the log2 scale. This behavior can be turned off by specifying
#' \code{register = FALSE}.
#'  
#' @export
#' @importFrom stats mad median sd
#' @importFrom dplyr group_by group_by_at mutate ungroup summarize left_join select
#' @author Patrick Kimes
cy3GenerateRef <- function(cy3pe, assay = SummarizedExperiment::assayNames(cy3pe)[1],
                           offset = 1L, register = TRUE) {

    ## verify validity of cy3pe
    stopifnot(is(cy3pe, "PBMExperiment"))
    sum_cols <- paste0("probe_", c("median", "mean", "mad", "sd"))
    if (any(names(rowData(cy3pe)) %in% sum_cols)) {
        stop("rowData of 'cy3pe' has column names matching the following special names:\n",
             paste0(sum_cols, collapse = ", "), "\n",
             "Please rename these columns.")
    }
    
    ## tidy scan data for parsing
    cy3vals <- broom::tidy(cy3pe, assay, long = TRUE)
    
    ## scale scans
    if (register) {
        ## common multiplication factor to scale samples
        col_medians <- colMedians(as.matrix(SummarizedExperiment::assay(cy3pe, assay)), na.rm = TRUE)
        sfactor <- median(col_medians, na.rm = TRUE)
        
        cy3vals <- dplyr::group_by(cy3vals, cname)
        cy3vals <- dplyr::mutate(cy3vals, value = value / median(value, na.rm = TRUE))
        cy3vals <- dplyr::ungroup(cy3vals)
        cy3vals <- dplyr::mutate(cy3vals, value = value * sfactor)
    } else {
        sfactor <- NA
    }
    
    ## calculate metrics on log2 scale
    cy3vals <- dplyr::mutate(cy3vals, value = log2(value + offset))

    ## compute probe-level summary metrics
    cy3vals <- dplyr::group_by_at(cy3vals, .vars = cy3pe@probeCols)
    cy3vals <- dplyr::summarize(cy3vals,
                                probe_median = median(value, na.rm = TRUE),
                                probe_mean = mean(value, na.rm = TRUE),
                                probe_mad = mad(value, na.rm = TRUE),
                                probe_sd = sd(value, na.rm = TRUE))
    cy3vals <- dplyr::ungroup(cy3vals)

    ## create SummarizedExperiment object
    rowdat <- as.data.frame(rowData(cy3pe), optional = TRUE)
    cy3vals <- dplyr::left_join(dplyr::as_tibble(rowdat), cy3vals, by = cy3pe@probeCols)
    rowdat <- DataFrame(dplyr::select(cy3vals, -one_of(sum_cols)))
    refassay <- dplyr::select(cy3vals, one_of(sum_cols))
    refassay <- list(ref = DataFrame(refassay))

    ## capture call parameters, combine with formal args
    params <- match.call(expand.dots = FALSE)[-1]
    params <- replace(formals(), names(params), params)

    ## return results as a SummarizedExperiment object
    PBMExperiment(assays = refassay, rowData = rowdat,
                  metadata = list(sfactor = sfactor,
                                  params = params),
                  probeFilter = cy3pe@probeFilter,
                  probeTrim = cy3pe@probeTrim,
                  probeCols = cy3pe@probeCols)
}
