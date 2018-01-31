#' Read PBMExperiment Data
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
                                 by = c("Column", "Row", "Name", "ID"))

    ## probe intensities from GPR files
    assaydat <- DataFrame(dplyr::select(assay_table, -Column, -Row, -Name, -ID))
    names(assaydat) <- paste0("s", 1:ncol(assaydat))

    ## row/probe-level metadata from GPR files
    rowdat <- DataFrame(dplyr::select(assay_table, Column, Row, Name, ID))

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
#' Helper function for reading in each GPR file.
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
#' SummarizedExperiment object with one assay with one column,
#' and corresponding rowData
#' 
#' @details
#' Since the name of the value column can vary across samples,
#' we pull it out using its index (5th column)
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr select
#' @author Patrick Kimes
readGPR <- function(x, useMean = FALSE, useBackground = FALSE) {
    colt <- rep("-", 45)
    colt[c(2, 3)] <- 'i'
    colt[c(4, 5, 38, 40)] <- 'c'

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

    vals <- readr::read_tsv(x, skip = 35, col_types = paste(colt, collapse = ""))
    names(vals)[5] <- 'intensity'
    vals %>% dplyr::select(Column, Row, Name, ID, intensity)
}


#' Spatial Adjustment of Samples
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' computes the local intensity bias for each probe, and returns the
#' same SummarizedExperiment object with the bias subtracted from each probe.
#' Optionally, the per-probe bias can also be returned as an additional array
#' in the bias corrected SummarizedExperiment object. The spatial bias at each
#' probe is defined as the difference between the median intensity in a `k` by `k`
#' square region surrounding the probe, and the median intensity of all probes
#' across the array. This approach is taken directly from the original PBM
#' analysis pipeline described in Berger and Bulyk (Nature Protocols, 2008).
#' The size of the local region can be specified by the user, with a default
#' size of 15 x 15 (as used in the aforementioned publication). 
#' 
#' @param se SummarizedExperiment object containing PBM intensity data
#' @param k odd integer specifying size of region to use to for computing
#'        local bias (default = 15)
#' @param returnBias logical whether to include the spatial bias as an
#'        additional 'assay' (called 'spatialbias') in the returned
#'        SummarizedExperiment object.
#' 
#' @return
#' SummarizedExperiment object with spatially adjusted intensities.
#' Row and column metadata are copied from the original SummarizedExperiment
#' object.
#'
#' @export
#' @author Patrick Kimes
spatiallyAdjust <- function(se, k = 15, returnBias = TRUE) {
    if (k %% 2 == 0) {
        stop("Local window size, k, must be an odd value.")
    }
    stopifnot(is(se, "SummarizedExperiment")) 
    if (! all(c("Row", "Column") %in% names(rowData(se)))) {
        stop("Specified SummarizedExperiment does not have necessary ",
             "'Row' and 'Column' information in rowData to perform ",
             "spatial adjustment")
    }
    stopifnot("gpr" %in% assayNames(se))

    ## extract intensities for easier manipulation
    intensity <- assay(se, "gpr")
   
    ## if ID column present, only use de Bruijn probes (mask others)
    if ("ID" %in% names(rowData(se))) {
        intensity[grepl("^dBr_", rowData(se)$ID), ] <- NA_real_
    }

    ## add row/column indicies
    intensity <- cbind(rowData(se)[, c("Row", "Column")], intensity)

    ## compute spatial adjustment 
    intensity <- tibble::as_tibble(intensity)
    intensity <- tidyr::gather(intensity, sample, value, -Row, -Column)
    intensity <- dplyr::group_by(intensity, sample)
    intensity <- dplyr::do(intensity, spatialmedian = .wrapSA(., k))
    intensity <- unnest(intensity)

    ## spread back so samples are in separate columns
    intensity <- tidyr::spread(intensity, -Row, -Column)

    ## join to original rowData to get proper row orders 
    intensity <- dplyr::left_join(rowData(se), intensity, by = c("Row", "Column"))
    intensity <- dplyr::select_(intensity, paste0("-", names(rowData(se))))

    ## construct new SummarizedExperiment from input SummarizedExperiment
    new_se <- se
    
    ## add new assays
    if (returnBias) {
        new_assays <- list(gpr = , spatialbias = )
    } else {
        new_assays <- list(gpr = )
    }
    assays(new_se) <- new_assays
    
    ## add step to list
    if (! "steps" %in% names(metadata(new_se))) {
        metadata(new_se)$steps <- list()
    }
    metadata(new_se)$steps <- c(metadata(new_se)$steps, "spatial adjustment")
    
    return(new_se) 
}

## Helper function which takes a data.frame with 'value', 'Column', 'Row'
## columns and returns the spatial median within a 'k' by 'k' region
## surrounding each value, and returns the values as a similar data.frame
## with columns, 'value', 'Column', 'Row'.
##
## @param x data.frame with columns 'value', 'Column', 'Row'
## @param k size of local region for computing medians
##
## @return
## a data.frame with columns 'value', 'Column', 'Row'
## 
## @author Patrick Kimes
.wrapSA <- function(x, k) {
    y <- dplyr::select(x, value, Column, Row)
    y <- tidyr::spread(y, Column, value)
    y <- dplyr::select(y, -Row)
    y <- blockmedian(as.matrix(y), k, center = TRUE)
    y <- tibble::as_tibble(y) 
    y <- dplyr::mutate(y, Row = row_number())
    y <- tidyr::gather(y, Column, value, -Row)
    y <- dplyr::mutate(y, Column = as.integer(gsub("V", "", Column)))
    y
}



#' RMA Background Correction of Samples
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' performs RMA-style background correction on each sample, and returns the
#' same SummarizedExperiment object with the background corrected intensities.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data
#'
#' @return
#' SummarizedExperiment object with background corrected intensities.
#' 
#' @import preprocessCore
#' @export 
#' @author Patrick Kimes
rmaBackgroundCorrect <- function(se) {
    ## perform RMA background correction (normal, exponential mixture)
    new_assay <- rma.background.correct(as.matrix(assay(se, "gpr")))
    
    ## construct new SummarizedExperiment from input SummarizedExperiment
    new_se <- se

    ## replace assays
    assays(new_se) <- list(gpr = DataFrame(new_assay))
    
    ## add step to list
    if (! "steps" %in% names(metadata(new_se))) {
        metadata(new_se)$steps <- list()
    }
    metadata(new_se)$steps <- c(metadata(new_se)$steps, "background correction")

    return(new_se)
}


#' RMA Normalization of Samples
#'
#' Given PBM intensities stored as a SummarizedExperiment, this function
#' performs RMA-style normalization on each sample (currently, only quantile
#' normalization is supported) and returns thesame SummarizedExperiment object
#' with the normalized intensities.
#'
#' @param se SummarizedExperiment object containing PBM intensity data
#'
#' @return
#' SummarizedExperiment object with normalized intensities.
#'
#' @import preprocessCore
#' @export
#' @author Patrick Kimes
rmaNormalize <- function(se) {
    ## perform quantile normalization
    new_assay <- preprocessCore::normalize.quantiles(as.matrix(assay(se, "gpr")))
    
    ## construct new SummarizedExperiment from input SummarizedExperiment
    new_se <- se

    ## replace assays
    assays(new_se) <- list(gpr = DataFrame(new_assay))
    
    ## add step to list
    if (! "steps" %in% names(metadata(new_se))) {
        metadata(new_se)$steps <- list()
    }
    metadata(new_se)$steps <- c(metadata(new_se)$steps, "rma normalization")

    return(new_se)
}


#' Estimation of K-Mer Intensities
#'
#' This function uses stuff to do stuff.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data
#' @param kmers character vector of k-mers to predict
#' @param stdArray logical whether the array is a standard PBM array. If
#'        TRUE, then two processing steps will be applied prior to estimating
#'        k-mer intensities, namely: 1. subsetting on de Bruijn sequences as
#'        determined by the ID column of the SummarizedExperiment rowData, and
#'        2. trimming of sequences to only use the unique 36nt region.
#'        (default = TRUE)
#' @param ... parameters to be passed to glmnet call - will overwrite
#'        default paramaters
#' 
#' @return
#' SummarizedExperiment object with predicted affinities.
#'
#' @details
#' Default signature for the glmnet call is:
#' * family = 'gaussian'
#' * alpha = 0.5
#' * intercept = FALSE
#' * lambda = rev(0:50)
#' * lower.limits = 0
#' * thresh = 1e-10
#' 
#' @import Matrix
#' @import glmnet
#' @export
#' @author Patrick Kimes
predictkmers <- function(se, kmers = NULL, stdArray = TRUE, ...) {
    if (is.null(kmers)) {
        data(pbm_8mers)
        kmers <- pbm_8mers
    } else if (!is.vector(kmers, mode = "character")) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences to estimate.")
    }

    if (stdArray) {
        if (! "ID" %in% names(rowData(se))) {
            stop("Specified SummarizedExperiment object does not contain probe ID ",
                 "information in the rowData. This information must be provided as ",
                 "the 'ID' column in the rowData to filter on de Bruijn sequences.\n",
                 "If this is not a standard PBM array and no probe filtering or sequence",
                 "trimming is necessary, re-run the analysis with 'stdArray = FALSE'.")
        }
        ## subset on de Bruijn sequence probes
        se <- se[grepl("^dBr_", rowData(se)$ID), ]

        plens <- unique(nchar(rowData(se)$Sequence))
        if (length(plens) > 1 || plens != 60) {
            stop("Not all probe sequences in specified SummarizedExperiment object are 60nt.\n",
                 "If this is not a standard PBM array and no probe filtering or sequence",
                 "trimming is necessary, re-run the analysis with 'stdArray = FALSE'.")
        }
        tails <- Biostrings::subseq(rowData(se)$Sequence, 37, 60)
        tails <- unique(tails)
        if (length(tails) > 1) {
            stop("Not all probe sequences in specified SummarizedExperiment object have the ",
                 "same 24nt primer sequence.\n",
                 "If this is not a standard PBM array and no probe filtering or sequence",
                 "trimming is necessary, re-run the analysis with 'stdArray = FALSE'.")
        }
        rowData(se)$Sequence <- Biostrings::subseq(rowData(se)$Sequence, 1, 36)
    }
    
    ## find mapping between kmers and probes
    ovnames <- intersect(names(rowData(se)), c("Row", "Column", "ID", "Sequence"))
    kmermap <- mapkmers(rowData(se)[, ovnames], kmers)

    kmer_levels <- levels(kmermap$seq)
    
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
    vec_intensities <- as.matrix(assay(se, "gpr"))
    vec_intensities <- matrix(vec_intensities, ncol = 1)
    vec_intensities <- log2(vec_intensities)

    ## merge input parameters with default parameter for glmnet call
    glmnet_args <- list(x = designmat, y = vec_intensities,
                        family = 'gaussian', alpha = 0.5,
                        intercept = FALSE, lambda = rev(0:50),
                        lower.limits = 0, thresh = 1e-10)
    args_input <- list(...)
    glmnet_args <- replace(glmnet_args, names(args_input), args_input)

    ## fit elastic net model, constraint to non-negative
    efit <- do.call(glmnet, glmnet_args)
    
    ## extracting estimated affinities
    nlambda <- length(glmnet_args$lambda)
    coefs_all <- coef(efit)[-1, nlambda]
    coefs_est <- matrix(coefs_all[ 1:(ncol(se) * nkmers) ],
                        ncol = ncol(se), nrow = nkmers)
    coefs_probes <- coefs_all[(ncol(se) * nkmers) + 1:41944]

    ## turn coefficient estimates into table w/ rows = 8mers, cols = conditions
    coefs_est <- DataFrame(coefs_est)
    names(coefs_est) <- names(es)

    ## determine row data
    rowdat <- DataFrame(kmer = kmers)

    ## create new SummarizedExperiment
    SummarizedExperiment(assays = list(pred = coefs_est),
                         rowData = rowdat, colData = colData(se))
}


