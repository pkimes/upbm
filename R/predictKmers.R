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


## helper function to check if specified k-mer set is valid
checkKmers <- function(kmers, verb) {
    if (is.null(kmers)) {
        if (verb) {
            cat("!! Using the default 8-mer set.\n")
        }
        data(pbm_8mers)
        return(pbm_8mers)
    }
    if (!is.vector(kmers, mode = "character")) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences to estimate.")
    }
    return(kmers)
}


## helper function to check if probe sequences in SE are valid
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


## helper function to trim probe sequences in SE if specified
trimProbeSequences <- function(se, .trim) {
    if (is.null(.trim)) {
        return(se)
    }
    
    ## otherwise perform trimming
    if (length(.trim) != 2 || !is.numeric(.trim)) {
        stop("If probe sequences should be trimmed, set '.trim' to a vector of length ",
             "2 corresponding to start and end of regions to be kept.")
    }
    if (.trim[1] > .trim[2] || .trim[1] < 1) {
        stop("Invalid choice of '.trim'.")
    }
    plens <- nchar(rowData(se)$Sequence)
    if (any(plens < .trim[2])) {
        stop("Choice of '.trim' is longer than at least one sequence.")
    }
    
    rowData(se)$Sequence <- Biostrings::subseq(rowData(se)$Sequence,
                                               .trim[1], .trim[2])
    return(se)
}
