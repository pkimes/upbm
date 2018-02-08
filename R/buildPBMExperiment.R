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
#' @param filterFlags logical whether to replace intensity values at probes
#'        flagged manually or automatically as being low quality
#'         ('Bad': -100, 'Absent': -75, 'Not Found': -50) with NA. (default = TRUE)
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
#' @importFrom dplyr select left_join
#' @importFrom tibble as_tibble
#' @export
#' @author Patrick Kimes
buildPBMExperiment <- function(tab, useMean = FALSE, useBackground = FALSE, filterFlags = TRUE,
                               probes = NULL) {
    ## check validity of inputs
    stopifnot(is.data.frame(tab))
    stopifnot(c("vers", "gpr") %in% names(tab))
    stopifnot(is.logical(useMean))
    stopifnot(is.logical(useBackground))
    stopifnot(is.logical(filterFlags))
    
    ## currently only support all scans with same design
    if (length(unique(tab$vers)) > 1) {
        stop("All samples/scans must have the same assay design version")
    }

    ## read in all GPR scans
    assay_table <- mapply(readGPR, gpr_path = tab$gpr, gpr_type = tab$scan,
                          useMean = useMean, useBackground = useBackground,
                          filterFlags = filterFlags,
                          SIMPLIFY = FALSE)
    assay_table <- purrr::reduce(assay_table, left_join,
                                 by = c("Column", "Row"))

    ## probe intensities from GPR files
    assaydat <- DataFrame(dplyr::select(assay_table, -Column, -Row))
    names(assaydat) <- paste0("s", 1:ncol(assaydat))

    ## row/probe-level metadata from GPR files
    rowdat <- dplyr::select(assay_table, Column, Row)

    ## check probes
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
            rowdat <- dplyr::left_join(rowdat, tibble::as_tibble(as.data.frame(probes)),
                                       by = ovnames)
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
#' @param gpr_path path to GPR fle.
#' @param gpr_type scan type; one of "Alexa488", "Cy3", "Masliner".
#' @param useMean logical whether to use mean fluorescent intensity
#'        for each probe rather than median fluorescent intensity.
#'        (default = FALSE)
#' @param useBackground logical whether to use background subtracted
#'        intensity rather than non-subtracted intensity.
#'        (default = FALSE)
#' @param filterFlags logical whether to replace intensity values at probes
#'        flagged manually or automatically as being low quality.
#'         ('Bad': -100, 'Absent': -75, 'Not Found': -50) with NA. (default = TRUE)
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
readGPR <- function(gpr_path, gpr_type, useMean = FALSE,
                    useBackground = FALSE, filterFlags = TRUE) {
    ## determine number of columns
    p <- 45
    if (gpr_type == "Masliner") {
        p <- 49
    }

    ## specify which columns to read in
    colt <- rep("-", p)
    colt[c(2, 3)] <- 'i'
    if (filterFlags) {
        colt[38] <- 'i'
    }
    
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

    vals <- readr::read_tsv(gpr_path, skip = 35, col_types = colt)
    names(vals)[3] <- 'intensity'

    ## remove negative (low quality) flagged probes
    if (filterFlags) {
        vals$intensity[vals$Flags < 0] <- NA_real_
        vals <- dplyr::select(vals, -Flags)
    }
    vals
}
