#' Test Probe Intensity Differences
#'
#' @description
#' This function uses limma to perform probe-level tests of intensity differences
#' across conditions. Contrasts should be specified as a design formula.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data.
#' @param design formula specifying design of test, where variables are
#'        taken from the \code{colData} of \code{se}.
#' @param assay_name string name of the assay to use. (default = "fore")
#' @param offset integer offset to add to intensities before log2 scaling to
#'        prevent errors with zero intensities. If set to 0, probes with
#'        zero intensities are dropped/ignored in estimation. (default = 1)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to estimating affinities. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @param baselines named list of baseline conditions for any categorical
#'        variables included in the model. Ignored if NULL. (default = NULL)
#' 
#' @return
#' SummarizedExperiment object with probe-level testing results. Each assay
#' contains the results for a single coefficient in the model specified
#' in \code{design}.
#'
#' @examples
#' \dontrun{
#' res <- probeTest(mygpr, ~ 1 + condition)
#' }
#' 
#' @md
#' @import SummarizedExperiment
#' @importFrom dplyr left_join distinct
#' @importFrom limma lmFit eBayes
#' @export
#' @author Patrick Kimes
probeTest <- function(se, design, assay_name = "fore", offset = 1L, .filter = 1L,
                       baselines = NULL) {

    stopifnot(inherits(design, "formula"))
    stopifnot(is.null(baselines) || is.list(baselines))
    
    ## filter probes
    se <- pbmFilterProbes(se, .filter) 

    ## extract probe data matrix
    datp <- assay(se, assay_name)
    datp <- log2(as.matrix(datp) + offset)
    datp[is.infinite(datp)] <- NA

    ## extact sample metadata
    metadat <- colData(se)
    
    ## processdesign matrix
    if (!is.null(baselines)) {
        for (ic in intersect(names(baselines), colnames(metadat))) {
            iclevs <- unique(metadat[[ic]])
            stopifnot(baselines[[ic]] %in% iclevs)
            iclevs <- c(baselines[[ic]], setdiff(iclevs, baselines[[ic]]))
            metadat[[ic]] <- factor(metadat[[ic]], levels = iclevs)
        }
    }
    mmat <- model.matrix(design, metadat)

    ## fit limma model with variance trend and robust estimation
    fit <- lmFit(datp, mmat)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)

    ## return SummarizedExperiment with test results - one coef per assay
    alist <- list(Amean = replicate(ncol(fit$coefficients), fit$Amean),
                  t = fit$t,
                  coefs = fit$coefficients,
                  stdev = sweep(fit$stdev.unscaled, 1, fit$s2.post^.5, `*`),
                  df = replicate(ncol(fit$coefficients), fit$df.total))
    alist <- simplify2array(alist)
    alist <- lapply(seq_len(dim(alist)[2]), function(i) alist[, i, ])
    names(alist) <- colnames(fit$coefficients)
    
    SummarizedExperiment(assays = alist, rowData = rowData(se))
}


#' Test K-mers Probe Differences
#'
#' @description
#' After performing probe-level testing, this function applies probe set
#' aggregation to obtain K-mer level inference. 
#'
#' @param se probe-level testing results
#' @param kmers character vector of k-mers to predict.
#' @param assay_name string name of the assay to aggregate. If NULL, all assays
#'        are aggregated. (default = NULL)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to estimating affinities. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @param .trim interger vector of length two specifying start and end
#'        of probe sequence to be used. Default is based on the universal
#'        PBM probe design where only leading 36nt should be used. 
#'        Ignored if NULL. (default = c(1, 36))
#'
#' @return
#' SummarizedExperiment of estimated K-mer testing results.
#'
#' @importFrom stats p.adjust
#' @importFrom dplyr select_ group_by left_join ungroup do mutate
#' @importFrom tidyr unnest
#' @importFrom rlang enquo
#' @export
#' @author Patrick Kimes
kmerTest <- function(se, kmers, assay_name = NULL, .filter = 1L,
                     .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {

    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot(is.null(assay_name) || all(assay_name %in% assayNames(se)))
    
    if (is.null(assay_name)) {
        assay_name <- assayNames(se)
    }

    ## check kmers specified
    kmers <- checkKmers(kmers, verb = FALSE)
    
    ## check Sequence info in rowData
    se <- checkProbeSequences(se, verb = FALSE)
    
    ## filter probes
    se <- pbmFilterProbes(se, .filter) 

    ## trim probe sequences
    se <- trimProbeSequences(se, .trim)

    ## find mapping between kmers and probes
    ovnames <- intersect(names(rowData(se)), c("Row", "Column", "ID", "Sequence"))
    kmermap <- mapKmers(rowData(se)[, ovnames, drop = FALSE], kmers)

    ## use ordering from input 'kmers'
    kmermap$seq <- factor(kmermap$seq, levels = kmers)
    
    adat <- lapply(assay_name, function(x, ...) { assay2tidy(assay_name = x, ...) },
                   se = se, long = FALSE, .filter = 0L)
    names(adat) <- assay_name
    adat <- dplyr::bind_rows(adat, .id = "aname")

    adat <- dplyr::select_(adat, .dots = c(setdiff(ovnames, "Sequence"), "Amean", "t", "aname"))
    adat <- dplyr::left_join(kmermap, adat, by = setdiff(ovnames, "Sequence"))

    adat <- dplyr::group_by(adat, aname, seq)
    adat <- dplyr::do(adat, zt = t.test(.$t),
                      nprobes = nrow(.),
                      Amean = mean(.$Amean, na.rm = TRUE))
    adat <- dplyr::ungroup(adat)
    adat <- tidyr::unnest(adat, nprobes, Amean)

    return(adat)
    
    adat <- dplyr::mutate(adat,
                          pval = sapply(zt, `[[`, "p.value"),
                          tstat = sapply(zt, `[[`, "statistic"),
                          effect = sapply(zt, `[[`, "estimate"),
                          df = sapply(zt, function(x) { x$parameter[["df"]] }),
                          zt = NULL)

    alist <- list(Amean = .tidymat(adat, kmers, Amean),
                  pval = .tidymat(adat, kmers, pval),
                  tstat = .tidymat(adat, kmers, tstat),
                  effect = .tidymat(adat, kmers, effect),
                  df = .tidymat(adat, kmers, df))

    alist <- simplify2array(alist)
    alist <- lapply(seq_len(dim(alist)[2]), function(i) alist[, i, ])
    names(alist) <- sort(assay_name)

    rdat <- dplyr::select(adat, seq, nprobes)
    rdat <- rdat[match(kmers, rdat$seq), ]
    
    SummarizedExperiment(assays = alist, rowData = rdat)
}

## helper to turn tidy table column into matrix for SE
.tidymat <- function(x, km, cn) {
    cn <- rlang::enquo(cn)
    x <- dplyr::select(x, seq, aname, !!cn)
    x <- tidyr::spread(x, aname, !!cn)
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
