#' @title Perform Cy3 normalization
#'
#' @description
#' Given a PBMExperiment containing Alexa488 probe intensity data and a second
#' PBMExperiment containing Cy3 observed-to-expected ratios, this function
#' performs probe filtering and/or scaling of the Alexa488 data according to the
#' over/under abundance of dsDNA at each probe as quantified from the Cy3 data.
#'
#' Both the Alexa488 and Cy3 PBMExperiment objects must have a common column in the
#' colData (specified by \code{match_by=}) that can be used to match scans across the
#' two objects. While the values of the \code{match_by} column may be repeated across
#' Alexa488 scans (e.g. if an array was reused or scanned multiple times), the
#' values must be unique across the Cy3 PBMExperiment object.
#' 
#' @param pe a PBMExperiment object containing Alexa488 intensity data.
#' @param cy3pe a PBMExperiment object containing Cy3 deviation values from
#'        using \code{cy3FitEmpirical} or \code{cy3FitModel}.
#' @param assay a numeric index or string name specifying the assay to use.
#'        (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param match_by a string column name in colData of \code{pe} and \code{cy3pe}
#'        to use for matching scans across the two objects; values of
#'        column must be unique for each scan in each object. (default = \code{"id_idx"})
#' @param filter a logical value whether to filter "low quality" flagged probes by setting to NA.
#'        (default = TRUE)
#' @param scale a logical value whether to perform scaling by ratio of observed vs. expected
#'        Cy3 intensity at each probe. (default = TRUE)
#' @param verbose a logical value whether to print verbose output during
#'        analysis. (default = FALSE)
#'
#' @return
#' Original PBMExperiment object with assay containing Cy3 filtered and/or scaled intensities
#' (\code{"normalized"}). If an assay with the same name is already included in the object, it will
#' be overwritten.
#'
#' @seealso \code{\link{cy3FitEmpirical}}
#' @importFrom dplyr as_tibble mutate select left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @export
#' @author Patrick Kimes
cy3Normalize <- function(pe, cy3pe, assay = SummarizedExperiment::assayNames(pe)[1],
                         match_by = "id_idx", filter = TRUE, scale = TRUE, verbose = FALSE) {

    stopifnot(is(pe, "PBMExperiment")) 
    stopifnot(is(cy3pe, "PBMExperiment")) 
    stopifnot("ratio" %in% SummarizedExperiment::assayNames(cy3pe))
    
    ## check validity of matching colData column
    stopifnot(match_by %in% names(colData(pe)))
    stopifnot(match_by %in% names(colData(cy3pe)))
    match_vals1 <- colData(pe)[[match_by]]
    match_vals2 <- colData(cy3pe)[[match_by]]
    if (any(duplicated(match_vals2))) {
        stop("Matching variable '", match_by, "' is not unique across samples in 'cy3pe'.\n",
             "Specify a different column in colData.")
    }

    match_overlap <- intersect(match_vals1, match_vals2)
    match_diff1 <- setdiff(match_vals1, match_vals2)
    match_diff2 <- setdiff(match_vals2, match_vals1)
    
    if (length(match_overlap) == 0) {
        stop("No samples matched across the specified two experiments.")
    }
    
    if (verbose) {
        cat("|| upbm::cy3Normalize\n")
        cat("|| - Starting probe-level Cy3 normalization and filtering",
            "for", ncol(pe), "Alexa488 PBM scans.\n")
    }

    if (!filter & !scale) {
        if (verbose) {
            cat("|| - No filtering or normalization performed.\n")
            cat("|| - Finished probe-level Cy3 normalization and filtering.\n")
        }
        return(pe)
    }

    if (length(match_diff1) > 0) {
        warning("Dropping samples from pe:\n  ", paste(match_diff1, collapse = ", "))
    }
    if (length(match_diff2) > 0) {
        warning("Dropping samples from cy3pe:\n  ", paste(match_diff2, collapse = ", "))
    }

    ## filter samples
    pe <- pe[, colData(pe)[[match_by]] %in% match_overlap, drop = FALSE]
    cy3pe <- cy3pe[, colData(cy3pe)[[match_by]] %in% match_overlap, drop = FALSE]

    ## condition must be a unique column for faceting plot
    coldat1 <- data.frame(colData(pe), check.names = FALSE,
                          check.rows = FALSE, stringsAsFactors = FALSE)
    coldat1 <- tibble::rownames_to_column(coldat1, "sample")
    coldat1 <- dplyr::rename(coldat1, Match = !!match_by)
    coldat1 <- dplyr::select(coldat1, sample, Match)

    coldat2 <- data.frame(colData(cy3pe), check.names = FALSE,
                          check.rows = FALSE, stringsAsFactors = FALSE)
    coldat2 <- tibble::rownames_to_column(coldat2, "sample")
    coldat2 <- dplyr::rename(coldat2, Match = !!match_by)
    coldat2 <- dplyr::select(coldat2, sample, Match)

    new_assay <- as.data.frame(SummarizedExperiment::assay(pe, assay), optional = TRUE)
    new_assay <- dplyr::as_tibble(new_assay)
    new_assay <- dplyr::mutate(new_assay,
                               Row = rowData(pe)[, "Row"],
                               Column = rowData(pe)[, "Column"])
    new_assay <- tidyr::pivot_longer(new_assay, names_to = "sample",
                                     values_to = "value", c(-Row, -Column))
    new_assay <- dplyr::left_join(new_assay, coldat1, by = "sample")
    ##new_assay <- dplyr::select(new_assay, -sample)

    if (scale) {
        if (verbose) {
            cat("|| - Performing Cy3 filtering ('filter = TRUE' specified).\n")
        }
        scale_assay <- as.data.frame(SummarizedExperiment::assay(cy3pe, "ratio"), optional = TRUE)
        scale_assay <- dplyr::as_tibble(scale_assay)
        scale_assay <- dplyr::mutate(scale_assay,
                                   Row = rowData(cy3pe)[, "Row"],
                                   Column = rowData(cy3pe)[, "Column"])
        scale_assay <- tidyr::pivot_longer(scale_assay, names_to = "sample",
                                           values_to = "scalar", c(-Row, -Column))
        scale_assay <- dplyr::left_join(scale_assay, coldat2, by = "sample")
        scale_assay <- dplyr::select(scale_assay, -sample)
        
        new_assay <- dplyr::left_join(new_assay, scale_assay, by = c("Row", "Column", "Match"))
        new_assay <- dplyr::mutate(new_assay, value = value / scalar)
        new_assay <- dplyr::select(new_assay, -scalar)
    } else if (verbose) {
        cat("|| - No filtering performed ('filter = FALSE' specified).\n")
    }
        
    if (filter) {
        if (verbose) {
            cat("|| - Performing Cy3 normalization ('scale = TRUE' specified).\n")
        }
        filter_assay <- as.data.frame(SummarizedExperiment::assay(cy3pe, "lowq"), optional = TRUE)
        filter_assay <- dplyr::as_tibble(filter_assay)

        filter_assay <- dplyr::mutate(filter_assay,
                                      Row = rowData(cy3pe)[, "Row"],
                                      Column = rowData(cy3pe)[, "Column"])
        filter_assay <- tidyr::pivot_longer(filter_assay, names_to = "sample",
                                            values_to = "flag", c(-Row, -Column))
        filter_assay <- dplyr::left_join(filter_assay, coldat2, by = "sample")
        filter_assay <- dplyr::select(filter_assay, -sample)

        ## also filter anything NA in Cy3
        filter_assay$flag[is.na(filter_assay$flag)] <- FALSE
        new_assay <- dplyr::left_join(new_assay, filter_assay, by = c("Row", "Column", "Match"))
        new_assay <- dplyr::mutate(new_assay, value = ifelse(flag, NA, value))
        new_assay <- dplyr::select(new_assay, -flag)
    } else if (verbose) {
        cat("|| - No normalization performed ('scale = FALSE' specified).\n")
    }

    ## return to square assay shape
    new_assay <- dplyr::select(new_assay, sample, value, Row, Column)
    new_assay <- tidyr::pivot_wider(new_assay, names_from = sample, values_from = value)

    ## match row order to rowData
    c_order <- paste(rowData(pe)$Row, rowData(pe)$Column, sep = "-")
    new_order <- match(c_order, paste(new_assay$Row, new_assay$Column, sep = "-"))
    stopifnot(!duplicated(new_order), length(new_order) == nrow(pe))
    new_assay <- new_assay[new_order, ]
    new_assay <- dplyr::select(new_assay, -Row, -Column)

    ## match column order to colData
    stopifnot(colnames(new_assay) %in% colnames(pe))
    new_assay <- new_assay[, colnames(pe)]

    ## clear any existing "normalized" assay and overwite with new normalized data
    if ("normalized" %in% assayNames(pe)) {
        SummarizedExperiment::assay(pe, "normalized") <- NULL
        if (verbose) {
            cat("|| - Original PBMExperiment object contains \"normalized\" assay.\n")
            cat("|| - Existing \"normalized\" assay will be overwritten.\n")
        }
    }
    SummarizedExperiment::assays(pe) <- c(S4Vectors::SimpleList(normalized = DataFrame(new_assay, check.names = FALSE)),
                                          SummarizedExperiment::assays(pe))
    
    if (verbose) {
        cat("|| - Finished probe-level Cy3 normalization and filtering.\n")
        cat("|| - Returning PBMExperiment with", nrow(pe), "rows and", ncol(pe), "columns.\n")
    }
    return(pe)
}
