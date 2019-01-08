#' Test k-mers Differences
#'
#' @description
#' After estimating k-mer affinities from probe sets using \code{kmerFit},
#' this function computes affinity differences across conditions, relative to
#' a baseline condition (typically the reference condition). 
#'
#' @param se SummarizedExperiment of K-mer affinities.
#' @param baseline string name of baseline condition to compare other conditions against.
#'        If NULL, baseline is guessed by looking for 'ref' in any of the conditions.
#'        (default = NULL)
#'
#' @return
#' SummarizedExperiment of estimated k-mer affinity differences.
#'
#' @importFrom stats p.adjust
#' @importFrom dplyr select_ group_by left_join ungroup do mutate
#' @importFrom tidyr unnest
#' @importFrom rlang enquo
#' @export
#' @author Patrick Kimes
kmerTest <- function(se, baseline = NULL) {

    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot(is.null(baseline) || is(baseline, "character"))
    
    bdat <- assay2tidy(se, "betaME", .filter = 0L)
    bdat <- tidyr::gather(bdat, condition, value, -seq, -(!!baseline))
    bdat <- dplyr::rename(bdat, "baseline" = !!baseline)
    bdat <- dplyr::mutate(bdat, M = value - baseline,
                          A = (value + baseline) / 2)
    bdat <- dplyr::select(bdat, seq, condition, M, A)

    vdat <- assay2tidy(se, "varME", .filter = 0L)
    vdat <- tidyr::gather(vdat, condition, value, -seq, -(!!baseline))
    vdat <- dplyr::rename(vdat, "baseline" = !!baseline)
    vdat <- dplyr::mutate(vdat, deltaVar =  value + baseline)
    vdat <- dplyr::select(vdat, seq, condition, deltaVar)

    dat <- dplyr::left_join(bdat, vdat, by = c("seq", "condition"))

    kmers <- rowData(se)$seq

    alist <- list(M = .tidycol2mat(dat, kmers, "M"),
                  A = .tidycol2mat(dat, kmers, "A"),
                  v = .tidycol2mat(dat, kmers, "deltaVar"))

    SummarizedExperiment(assays = alist, rowData = rowData(se))
}

## helper to turn tidy table column into matrix for SE
.tidycol2mat <- function(x, km, cn) {
    x <- dplyr::select(x, seq, condition, !!cn)
    x <- tidyr::spread(x, condition, !!cn)
    x <- x[match(km, x$seq), sort(names(x))]
    as.matrix(dplyr::select(x, -seq))
}


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
        data(pbm_8mers)
        return(pbm_8mers)
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
    data(pbm_8x60k_v1)
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


#' Trim Probe Sequences 
#'
#' Simple helper function to trim probe sequences in the rowData of
#' a SummarizedExperiment object, or alternatively, a DataFrame/data.frame
#' object. The sequences should be in a column named 'Sequence'.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data and
#'        probe sequence information in rowData, or simply the DataFrame/data.frame
#'        corresponding to the rowData.
#' @param .trim interger vector of length two specifying start and end
#'        of probe sequence \strong{to be used}. Default is based on the universal
#'        PBM probe design where only leading 36nt should be used. 
#'        Ignored if \code{NULL}. (default = \code{c(1, 36)})
#'
#' @return
#' Original \code{se} object with trimmed Sequence column.
#'
#' @importFrom Biostrings subseq
#' @export
#' @author Patrick Kimes
trimProbeSequences <- function(se, .trim = c(1, 36)) {
    if (is.null(.trim)) {
        return(se)
    }

    if (is(se, "SummarizedExperiment")) {
        pd <- rowData(se)
    } else if (is(se, "DataFrame") | is(se, "data.frame")) {
        pd <- se
    } else {
        stop("se must be a SummarizedExperiment or DataFrame/data.frame")
    }

    if (! "Sequence" %in% names(pd)) {
        stop(paste0("Must have 'Sequence' column ",
                    ifelse(is(se, "SummarizedExperiment"), "in rowData ", ""),
                    "to perform trimming."))
    }
    
    if (length(.trim) != 2 || !is.numeric(.trim)) {
        stop("If probe sequences should be trimmed, set '.trim' to a vector of length ",
             "2 corresponding to start and end of regions to be kept.")
    }
    if (.trim[1] > .trim[2] || .trim[1] < 1) {
        stop("Invalid choice of '.trim'.")
    }
    plens <- nchar(pd$Sequence)
    if (any(plens < .trim[2])) {
        stop("Choice of '.trim' is longer than at least one sequence.")
    }

    ## if made it here, trim
    if (is(se, "SummarizedExperiment")) {
        rowData(se)$Sequence <- Biostrings::subseq(pd$Sequence, .trim[1], .trim[2])
    } else {
        se$Sequence <- Biostrings::subseq(pd$Sequence, .trim[1], .trim[2])
    }

    return(se)
}


#' Simple wrapper function to process PBM probe sequences
#'
#' @param se SummarizedExperiment object containing PBM intensity data and
#'        probe sequence information in rowData, or simply the DataFrame/data.frame
#'        corresponding to the rowData.
#' @param verbose logical whether to print extra messages. (default = FALSE)
#' @param .filter integer specifying level of probe filtering to perform. See
#'        \code{pbmFilterProbes} for more details about levels of probe filtering.
#'        (default = 1L)
#' @param .trim interger vector of length two specifying start and end
#'        of probe sequence \strong{to be used}. Default is based on the universal
#'        PBM probe design where only leading 36nt should be used. 
#'        Ignored if \code{NULL}. (default = \code{c(1, 36)})
#'
#' @return 
#' Original \code{se} object with only probes passing filter, with trimmed Sequence columns.
#' If the probe sequences are not valid, an error is thrown. 
#'
#' @seealso trimProbeSequences pbmFilterProbes checkProbeSequences
#' @export
#' @author Patrick Kimes
pbmProcessProbes <- function(se, verbose = FALSE, .filter = 1L, .trim = c(1, 36)) {
    ## check Sequence info in rowData
    se <- checkProbeSequences(se, verbose)
    
    ## filter probes
    se <- pbmFilterProbes(se, .filter) 

    ## trim probe sequences
    trimProbeSequences(se, .trim)
}
