#' Read PBMExperiment Data
#'
#' Read PBM data from a table containing paths to GPR files. 
#' 
#' @param tab table of samples with at least `gpr` and `vers` columns
#'        corresponding to GPR file paths and the corresponding PBM
#'        design version
#' @param useMean logical whether to use mean fluorescent intensity
#'        for each probe rather than median fluorescent intensity
#'        (default = FALSE)
#' @param useBackground logical whether to use background subtracted
#'        intensity rather than non-subtracted intensity
#'        (default = FALSE)
#' @param probes data.frame containing probe sequences in a column, 'Sequence',
#'        or a character vector specifying the probe sequences. If specified,
#'        these values will be added to the rowData of the returned
#'        SummarizedExperiment object (default = NULL)
#'
#' @return
#' SummarizedExperiment object with a single assay, 'gpr', containing the probe
#' intensities for the samples specified in the input 'tab'.
#'
#' @import SummarizedExperiment
#' @importFrom purrr reduce
#' @importFrom dplyr select
#' @export
#' @author Patrick Kimes
makePBMExperiment <- function(tab, useMean = FALSE, useBackground = FALSE, probes = NULL) {
    ## check validity of inputs
    stopifnot(is.data.frame(tab))
    stopifnot(c("vers", "gpr") %in% names(tab))
    stopifnot(is.logical(useMean))
    stopifnot(is.logical(useBackground))
    if (!is.null(probes)) {
        if (is.vector(probes, mode = "character")) {
            probes <- DataFrame(Sequence = probes)
        }
        if (!is(probes, "DataFrame") & !is(probes, "data.frame")) {
            warning("Specified 'probes' must be a DataFrame, data.frame, or character vector ",
                    "of probe sequences equal to the size of each GPR array result.\n",
                    "Ignoring specified 'probes' input.")
            probes <- NULL
        } else if (nrow(probes) != nrow(rowdat)) {
            warning("Dimension of specified 'probes' does not match GPR files.\n",
                    "Ignoring specified 'probes' input.")
            probes <- NULL
        } else if (! "Sequence" %in% names(probes)) {
            warning("Specified 'probes' must have a column named 'Sequence' with the probe sequences.\n",
                    "Ignoring specified 'probes' input.")
            probes <- NULL
        }
    }
    
    ## currently only support all scans with same design
    if (length(unique(tab$vers)) > 1) {
        stop("All samples/scans must have the same assay design version")
    }

    ## read in all GPR scans
    assay_table <- lapply(tab$gpr, readGPR, useMean = useMean, useBackground = useBackground)
    assay_table <- purrr::reduce(assay_table, left_join,
                                 by = c("Column", "Row"))

    ## probe intensities from GPR files
    assaydat <- DataFrame(dplyr::select(assay_table, -Column, -Row))
    names(assaydat) <- paste0("s", 1:ncol(assaydat))

    ## row/probe-level metadata from GPR files
    rowdat <- DataFrame(dplyr::select(assay_table, Column, Row))

    ## add probes to metadata if provided
    if (!is.null(probes)) {
        ovnames <- intersect(names(rowdat), names(probes))
        if (length(ovnames) == 0) {
            warning("Specified 'probes' does not have columns matching GPR row metadata.\n",
                    "The probes will be added to the table in the order provided.\n",
                    "WARNING: This may lead to incorrect joining.")
            rowdat <- cbind(rowdat, probes)
        } else {
            if (!all(c("Column", "Row") %in% ovnames)) {
                warning("Specified 'probes' does not have 'Column' and 'Row' columns.\n",
                        "The probes will be matched using other columns.\n",
                        "WARNING: This may lead to incorrect joining.")
            }
            rowdat <- merge(rowdat, probes, by = ovnames, all.x = TRUE)
        }
    }
    
    ## column/condition-level metadata
    coldat <- DataFrame(dplyr::select(tab, -gpr))
    rownames(coldat) <- names(assaydat)
    
    ## SummarizedExperiment
    SummarizedExperiment(assays = list(gpr = assaydat),
                         rowData = rowdat, colData = coldat)
}


#' Read GPR File as Assay
#' 
#' Helper function for reading in a single GPR file.
#'
#' @param x path to GPR fle
#' @param useMean logical whether to use mean fluorescent intensity
#'        for each probe rather than median fluorescent intensity
#'        (default = FALSE)
#' @param useBackground logical whether to use background subtracted
#'        intensity rather than non-subtracted intensity
#'        (default = FALSE)
#'
#' @return
#' tibble (data.frame-like) object of a single GPR file with three
#' columns: 'Column', 'Row', 'intensity'. ('ID' and 'Name' columns are
#' ignored as these may be incorrect in the GPR file.) 
#' 
#' @details
#' Since the name of the value column can vary across samples,
#' we pull it out using its index. Therefore, if the order or number of
#' columns in the GPR file is altered, reading may fail. 
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr select
#' @author Patrick Kimes
readGPR <- function(x, useMean = FALSE, useBackground = FALSE) {
    colt <- rep("-", 45)
    colt[c(2, 3)] <- 'i'
    colt[c(38, 40)] <- 'c'

    ## column corresponding to intensities
    if (useMean) {
        value_idx <- 10
    } else {
        value_idx <- 9
    }
    if (useBackground) {
        value_idx <- value_idx + 25
    }
    colt[c(value_idx)] <- 'd'
    colt <- paste(colt, collapse = "")

    vals <- readr::read_tsv(x, skip = 35, col_types = colt)
    names(vals)[3] <- 'intensity'
    vals
}
