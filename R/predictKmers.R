#' Estimate K-mer Affinities
#'
#' This function uses elastic net regularized least squares regression
#' to estimate k-mer affinities from PBM data. The base model
#' estimates affinities across multiple samples with a shared probe-level
#' effect. Estimates are reported on the scale of log2 relative affinity.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data.
#' @param assay_name string name of the assay to use. (default = "fore")
#' @param kmers character vector of k-mers to predict.
#' @param offset integer offset to add to intensities before log2 scaling to
#'        prevent errors with zero intensities. If set to 0, probes with
#'        zero intensities are dropped/ignored in estimation. (default = 1)
#' @param verbose logical whether to print extra messages during model fitting
#'        procedure. (default = FALSE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to estimating affinities. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @param .trim interger vector of length two specifying start and end
#'        of probe sequence to be used. Default is based on the universal
#'        PBM probe design where only leading 36nt should be used. 
#'        Ignored if NULL. (default = c(1, 36))
#' @param ... parameters to be passed to glmnet call - will overwrite
#'        default paramaters.
#' 
#' @return
#' SummarizedExperiment object with predicted affinities.
#'
#' @details
#' The default signature for the glmnet call is:
#' * family = 'gaussian'
#' * alpha = 0.5
#' * intercept = FALSE
#' * lambda = rev(0:50)
#' * lower.limits = 0
#' * thresh = 1e-10
#'
#' Features to be added:
#' * option to select different regression models
#' * option to return glmnet results for all alpha
#' * frozen RMA-like pooled estimate?
#'
#' @md
#' @import Matrix glmnet SummarizedExperiment
#' @importFrom Biostrings subseq
#' @importFrom utils data
#' @importFrom stats coef
#' @importFrom methods as is
#' @importFrom dplyr left_join distinct 
#' @export
#' @author Patrick Kimes
predictKmers <- function(se, assay_name = "fore", kmers = NULL, offset = 1,
                         verbose = FALSE, .filter = 1L,
                         .trim = if (.filter > 0L) { c(1, 36) } else { NULL },
                         ...) {
    
    ## check kmers specified
    kmers <- checkKmers(kmers, verbose)
    
    ## check Sequence info in rowData
    se <- checkProbeSequences(se, verbose)
    
    ## filter probes
    se <- pbmFilterProbes(se, .filter) 

    ## trim probe sequences
    se <- trimProbeSequences(se, .trim)

    ## find mapping between kmers and probes
    ovnames <- intersect(names(rowData(se)), c("Row", "Column", "ID", "Sequence"))
    kmermap <- mapKmers(rowData(se)[, ovnames, drop = FALSE], kmers)

    ## use ordering from input 'kmers'
    kmermap$seq <- factor(kmermap$seq, levels = kmers)

    ## convert kmer map to design matrix
    designmat <- sparseMatrix(i = as.numeric(kmermap$probe_idx),
                              j = as.numeric(kmermap$seq))
    nprobes <- dim(designmat)[1]
    nkmers <- dim(designmat)[2]

    ## check if design matrix dimensions match data dimensions
    if (nprobes != nrow(se)) {
        stop("Design matrix dimensions do not match probe data dimensions")
    }

    ## add variables for probe effect
    diagmat <- .sparseDiagonal(nprobes)
    
    ## create full design matrix w/ per-probe effects, per-sample-8mer effects
    designmat <- bdiag(replicate(ncol(se), designmat, simplify = FALSE))
    diagmat <- do.call(rbind, replicate(ncol(se), diagmat, simplify = FALSE))
    designmat <- cbind(designmat, diagmat)

    ## vectorize and log-transform intensities
    vec_intensities <- as.matrix(assay(se, assay_name))
    vec_intensities <- matrix(vec_intensities, ncol = 1)
    if (offset <= 0) {
        vec_intensities[vec_intensities <= 0] <- NA
        vec_intensities <- log2(vec_intensities)
    } else {
        vec_intensities <- log2(vec_intensities + offset)
    }
    
    ## remove NAs
    valid_idx <- which(!is.na(vec_intensities))
    vec_intensities <- vec_intensities[valid_idx, ]
    designmat <- designmat[valid_idx, ]
    
    ## merge input parameters with default parameter for glmnet call
    glmnet_args <- list(x = designmat, y = vec_intensities,
                        family = 'gaussian', alpha = 0.5,
                        intercept = FALSE, lambda = rev(0:50),
                        lower.limits = 0, thresh = 1e-10)
    args_input <- list(...)
    glmnet_args <- replace(glmnet_args, names(args_input), args_input)
    
    ## fit elastic net model, constraint to non-negative
    if (verbose) { print("starting glmnet fit ...") }
    efit <- do.call(glmnet, glmnet_args)
    if (verbose) { print("... finished!") }
    
    ## extracting estimated affinities
    nlambda <- length(glmnet_args$lambda)
    coefs_all <- coef(efit)[-1, nlambda]
    coefs_est <- matrix(coefs_all[ 1:(ncol(se) * nkmers) ],
                        ncol = ncol(se), nrow = nkmers)

    ## need to be careful about dropped K-mers if too many NAs
    
    ## turn coefficient estimates into table w/ rows = 8mers, cols = conditions
    coefs_est <- DataFrame(coefs_est)
    colnames(coefs_est) <- colnames(se)

    ## determine row data
    rowdat <- DataFrame(kmer = kmers)

    ## create new SummarizedExperiment
    SummarizedExperiment(assays = list(pred = coefs_est),
                         rowData = rowdat, colData = colData(se))
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
#'        probe sequence information in rowData.
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
    if ("Sequence" %in% names(rowData(se))) {
        return(se)
    }

    ## otherwise use default probe sequences
    if (verb) {
        cat("!! Using the default probe sequence set.\n")
    }
    if (! all(c("Row", "Column") %in% names(rowData(se)))) {
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
    rdat <- dplyr::left_join(as.data.frame(rowData(se), optional = TRUE),
                             dplyr::distinct(pbm_8x60k_v1),
                             by = c("Column", "Row"))
    rowData(se) <- DataFrame(rdat)

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
        pd <- rowData(pd)
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
        se <- Biostrings::subseq(pd$Sequence, .trim[1], .trim[2])
    }

    return(se)
}


#'
#'
#' Simple wrapper function to process PBM probe sequences
#'
#' @param se SummarizedExperiment object containing PBM intensity data and
#'        probe sequence information in rowData, or simply the DataFrame/data.frame
#'        corresponding to the rowData.
#' @param verb logical whether to print extra messages. (default = FALSE)
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
    se <- trimProbeSequences(se, .trim)
}
