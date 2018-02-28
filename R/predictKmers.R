#' Estimation of K-Mer Intensities
#'
#' This function uses stuff to do stuff.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data.
#' @param kmers character vector of k-mers to predict.
#' @param assay_name string name of the assay to use. (default = "fore")
#' @param std_array logical whether the array is a standard PBM array.
#'        See Details for more information. (default = TRUE)
#'        (default = TRUE)
#' @param verbose logical whether to print extra messages during model fitting
#'        procedure. (default = FALSE) 
#' @param ... parameters to be passed to glmnet call - will overwrite
#'        default paramaters.
#' 
#' @return
#' SummarizedExperiment object with predicted affinities.
#'
#' @details
#' If \code{stdAssay} is TRUE, prior to estimating k-mer affinities,
#' probes are subset to de Bruijn sequences as determined by the ID column
#' of the rowData, and probe sequences are trimmed to only use the unique
#' 36nt region.
#' 
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
#' @export
#' @author Patrick Kimes
predictKmers <- function(se, kmers = NULL, assay_name = "fore",
                         std_array = TRUE, verbose = FALSE, ...) {
    if (is.null(kmers)) {
        data(pbm_8mers)
        kmers <- pbm_8mers
        if (verbose) {
            cat("!! Using the default 8mer set.\n")
        }
    } else if (!is.vector(kmers, mode = "character")) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences to estimate.")
    }

    if (! "Sequence" %in% names(rowData(se))) {
        if (verbose) {
            cat("!! Specified SummarizedExperiment object does not contain probe sequence information in the rowData.\n",
                "!! The included default data(pbm_8x60k_v1) will be used.\n", sep = "")
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
        ##ovnames <- intersect(names(pbm_8x60k_v1), names(rowData(se)))
        ovnames <- c("Column", "Row")
        rdat <- dplyr::left_join(as.data.frame(rowData(se)), dplyr::distinct(pbm_8x60k_v1),
                                 by = ovnames)
        rowData(se) <- DataFrame(rdat)
    }
    
    if (std_array) {
        if (! "ID" %in% names(rowData(se))) {
            stop("Specified SummarizedExperiment object does not contain probe ID ",
                 "information in the rowData. This information must be provided as ",
                 "the 'ID' column in the rowData to filter on de Bruijn sequences.\n",
                 "If this is not a standard PBM array and no probe filtering or sequence",
                 "trimming is necessary, re-run the analysis with 'std_array = FALSE'.")
        }
        ## subset on de Bruijn sequence probes
        se <- se[grepl("^dBr_", rowData(se)$ID), ]

        plens <- unique(nchar(rowData(se)$Sequence))
        if (length(plens) > 1 || plens != 60) {
            stop("Not all probe sequences in specified SummarizedExperiment object are 60nt.\n",
                 "If this is not a standard PBM array and no probe filtering or sequence",
                 "trimming is necessary, re-run the analysis with 'std_array = FALSE'.")
        }
        tails <- Biostrings::subseq(rowData(se)$Sequence, 37, 60)
        tails <- unique(tails)
        if (length(tails) > 1) {
            stop("Not all probe sequences in specified SummarizedExperiment object have the ",
                 "same 24nt primer sequence.\n",
                 "If this is not a standard PBM array and no probe filtering or sequence",
                 "trimming is necessary, re-run the analysis with 'std_array = FALSE'.")
        }
        rowData(se)$Sequence <- Biostrings::subseq(rowData(se)$Sequence, 1, 36)
    }

    ## find mapping between kmers and probes
    ovnames <- intersect(names(rowData(se)), c("Row", "Column", "ID", "Sequence"))
    kmermap <- mapKmers(rowData(se)[, ovnames], kmers)

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
    vec_intensities <- log2(vec_intensities)

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


