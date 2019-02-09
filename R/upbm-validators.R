#' Check K-mer Sequences
#' 
#' Simple helper function to check if specified k-mer set is valid,
#' namely a vector of equal length character strings.
#' If \code{NULL}, a default set of 8-mers is returned. If invalid, an error is thrown.
#'
#' @param kmers character vector of k-mer sequences or \code{NULL}.
#' @param verb logical whether to print extra messages.
#'
#' @return
#' Original \code{se} object with default \code{pbm_8x60k_v1} probe design added
#' if probe 'Sequence' column is missing and the 8x60k design dimensions match
#' to original \code{se}.
#'
#' 
#' @importFrom dplyr left_join
#' @export
#' @author Patrick Kimes
checkKmers <- function(kmers, verb) {
    if (is.null(kmers)) {
        if (verb) {
            cat("!! Using the default 8-mer set.\n")
        }
        return(uniqueKmers(8L))
    }
    if (!is.vector(kmers, mode = "character")) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences as character strings.")
    }
    if (length(unique(nchar(kmers))) != 1L) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences of equal length.")
    }
    return(kmers)
}


#' Check Probe Sequences
#' 
#' Simple helper function to check if probe sequences in the rowData of
#' a SummarizedExperiment object are valid. The sequences should be in
#' a column named 'Sequence'.
#'
#' @param se SummarizedExperiment object containing PBM intensity data and
#'        probe sequence information in rowData, or simply the DataFrame/data.frame
#'        corresponding to the rowData.
#' @param verb logical whether to print extra messages.
#'
#' @return
#' Original \code{se} object with default \code{pbm_8x60k_v1} probe design added
#' if probe 'Sequence' column is missing and the 8x60k design dimensions match
#' to original \code{se}.
#'
#' @importFrom dplyr left_join
#' @export
#' @author Patrick Kimes
checkProbeSequences <- function(se, verb) {
    if (is(se, "SummarizedExperiment")) {
        pd <- rowData(se)
    } else if (is(se, "DataFrame") | is(se, "data.frame")) {
        pd <- se
    } else {
        stop("se must be a SummarizedExperiment or DataFrame/data.frame")
    }

    if ("Sequence" %in% names(pd)) {
        return(se)
    }

    ## otherwise use default probe sequences
    if (verb) {
        cat("!! Using the default probe sequence set.\n")
    }
    if (! all(c("Row", "Column") %in% names(pd))) {
        stop("The default set of probe sequences could not be used because the specified SummarizedExperiment ",
             "object does not contain probe 'Row' and 'Column' information in the rowData.\n",
             "The unique 'Row' and 'Column' values (and 'Sequence') information must be added to the rowData.")
    }
    data(pbm_8x60k_v1, package = "upbmAux")
    if (nrow(se) != nrow(pbm_8x60k_v1)) {
        stop("The default set of probe sequences could not be used because the specified SummarizedExperiment ",
             "contains ", nrow(se), " rows, and the default probe set includes ", nrow(pbm_8x60k_v1), " rows.\n",
             "The unique 'Sequence' information must be added to the rowData to use this command.")
    }
    pd <- dplyr::left_join(as.data.frame(pd, optional = TRUE),
                           dplyr::distinct(pbm_8x60k_v1),
                           by = c("Column", "Row"))
    if (is(se, "SummarizedExperiment")) {
        rowData(se) <- DataFrame(rdat)
    }
    
    return(se)
}


#' Trim probe sequences 
#'
#' Simple helper function to trim probe sequences in the rowData of
#' a PBMExperiment or PBMDesign object.
#' 
#' @param pe a PBMExperiment object or PBMDesign object.
#'
#' @return
#' Original PBMExperiment or PBMDesign object with trimmed sequences.
#'
#' @importFrom stringr str_sub
#' @export
#' @author Patrick Kimes
trimProbeSequences <- function(pe) {
    stopifnot(is(pe, "PBMExperiment") || is(pe, "PBMDesign"))
    
    if (length(pe@probeTrim) == 0L) {
        return(pe)
    }

    if (is(pe, "PBMExperiment")) {
        pd <- rowData(pe)
    } else {
        pd <- pe@design
    }
    pd$Sequence <- stringr::str_sub(pd$Sequence, pe@probeTrim[1], pe@probeTrim[2])

    return(pe)
}
