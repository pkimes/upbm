#' @title Fit probe models
#'
#' @description
#' Given a PBMExperiment of normalized probe intensities, this function fits
#' robust probe-level models across samples using the \code{lmFit} and \code{eBayes}
#' functions from the \pkg{limma} package. The output is a PBMExperiment object with
#' probes in rows and conditions, rather than individual samples, in columns. The
#' condition labels of each sample should be specified as a column in the colData of
#' input PBMExperiment using the \code{stratify=} parameter. By default, the function
#' expects a \code{"condition"} column in the colData of the input object.
#' 
#' @param pe a PBMExperiment object containing PBM intensity data.
#' @param assay a string name of the assay to use for fitting probe models.
#'        (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param stratify a character string specifying a column in \code{colData(pe)} to
#'        use for grouping samples, e.g. into alleles. (default = \code{"condition"})
#' @param guardrail a logical value whether to stop function if any sample appears
#'        to be of clearly low quality. Currently, this only checks for whether
#'        more than 20\% of probes in any sample are NA. If any criteria is met,
#'        an error will be returned. (default = TRUE)
#' @param offset an integer offset to add to intensities before log2 scaling to
#'        prevent errors with zero intensities. If set to 0, probes with
#'        zero intensities are dropped/ignored in estimation. (default = 1)
#' @param verbose a logical value whether to print verbose output during analysis. (default = FALSE)
#' @param ... additional parameters to be passed to \code{limma::eBayes}. See
#'        Details for default parameters used in this function as they differ from
#'        the \code{limma::eBayes} default values.
#' 
#' @return
#' PBMExperiment object with probe-level log2 intensity information aggregated
#' across replicates for each unique value in \code{stratify}. Columns in the returned
#' PBMExperiment correspond to unique values of the \code{stratify} variable, and
#' rows correspond to the subset of probes passing the filtering criteria
#' of the original PBMExperiment object. The PBMExperiment includes three assays,
#' \code{"beta"}, \code{"sd"}, and \code{"df"}, storing the cross-replicate
#' mean, standard deviation, and degrees of freedom estimates for each condition and probe.
#'
#' @details
#' More specifically, \code{limma::lmFit} is called on the log2-transformed values of the
#' specified assay with no intercept and a binary design matrix with columns corresponding
#' to the unique values of \code{colData(pe)[[stratify]]}. Since PBM experiments are often
#' performed with few replicates per condition, the fit is then passed to \code{limma::eBayes}
#' with \code{robust = TRUE} and \code{trend = TRUE} by default. Parameters to
#' \code{limma::eBayes} can be modified by specifying them during the call to \code{probeFit}.
#' 
#' @references
#' The cross-probe robust and trended variance estimates are obtained using the
#' \code{eBayes} function from the limma package with \code{robust = TRUE} and
#' \code{trend = TRUE} by default. The procedure was developed and described
#' in the following paper:
#' \itemize{
#' \item Phipson, B., Lee, S., Majewski, I. J., Alexander, W. S., & Smyth, G. K. (2016). Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. Annals of Applied Statistics, 10(2), 946-963.
#' }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames assay
#' @importFrom dplyr left_join distinct
#' @importFrom limma lmFit eBayes
#' @importFrom statmod gauss.quad.prob
#' @importFrom purrr quietly
#' @export
#' @author Patrick Kimes
probeFit <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                     stratify = "condition", guardrail = TRUE,
                     offset = 1L, verbose = FALSE, ...) {

    stopifnot(is(pe, "PBMExperiment"))
    stopifnot(assay %in% SummarizedExperiment::assayNames(pe))
    stopifnot(is(offset, "numeric"))
    
    if (verbose) {
        cat("|| upbm::probeFit \n")
        cat("|| - Starting probe-level model fitting for", ncol(pe), "PBM scans.\n")
    }

    ## define eBayes parameters with defaults
    dots <- list(...)
    eb_args <- list(trend = TRUE, robust = TRUE)
    eb_args <- replace(eb_args, names(dots), dots)

    if (verbose) {
        cat("|| - Filtering probes according to", length(pe@probeFilter),
            "probeFilter rule(s).\n")
        ntotal <- nrow(pe)
    }

    ## filter using rules
    pe <- pbmFilterProbes(pe)
    
    if (verbose) {
        cat("|| - Data filtered from", ntotal, "probes to", nrow(pe), "probes.\n")
    }

    ## extract probe data matrix
    datp <- SummarizedExperiment::assay(pe, assay)
    datp <- log2(as.matrix(datp) + offset)
    datp[is.infinite(datp)] <- NA

    ## apply guardrails checks
    if (guardrail) {
        gr_nna <- colMeans(is.finite(datp))
        n_nna <- sum(gr_nna < .80)
        if (n_nna > 0L) {
            stop(paste0("- Data includes ", n_nna,
                        " ", ifelse(n_nna > 1L, "samples", "sample"),
                        " with more than 20% of probes having non-finite",
                        " values.\n",
                        " - ", ifelse(n_nna > 1L, "columns", "column"), ": ",
                        paste0(colnames(datp)[gr_nna < .80], collapse = ", "), "."))
        }
    }
    
    ## check if any probes have no non-NA values, drop from limma testing
    nProbes <- nrow(datp)
    dropProbes <- rowSums(is.finite(datp)) == 0L
    if (any(dropProbes)) {
        datp <- datp[!dropProbes, , drop = FALSE]
    }

    if (verbose && (sum(dropProbes) > 0)) {
        cat("|| - Skipping", sum(dropProbes), "probes with NAs across all samples.\n")
    }
    
    ## check stratification params
    coldat <- .pbmExpandStratify(pe, stratify)
    
    if (verbose) {
        cat("|| - Fitting probe-level model with conditions:\n")
        cat("||     -", paste0(colnames(coldat), collapse = ", "), "\n")
    }

    ## fit limma model
    fit <- purrr::quietly(limma::lmFit)(datp, coldat)
    if (verbose && (length(fit$warnings) > 0L)) {
        for (i in seq_len(length(fit$warnings))) {
            cat("|| -[limma::lmFit]", fit$warnings[i], "\n")
        }
    }
    fit <- fit$result

    if (verbose && (length(eb_args) > 0)) {
        cat("|| - Performing eBayes shrinkage with:\n")
        for (i in seq_len(length(eb_args))) {
            cat("||     -", names(eb_args)[i], ":", eb_args[[i]], "\n") 
        }
    }

    ## fit empirical Bayes adjustment with specified parameters
    eb_args$fit <- fit
    fit <- do.call(limma::eBayes, eb_args)
    
    ## return SummarizedExperiment with test results - one coef per assay
    assayl <- list(beta = fit$coefficients,
                   sd = sweep(fit$stdev.unscaled, 1, fit$s2.post^.5, `*`),
                   df = replicate(ncol(fit$coefficients), fit$df.total))

    ## return to full dimensions
    if (any(dropProbes)) {
        assayl <- lapply(assayl, function(x) {
            y <- matrix(NA, nrow = nProbes, ncol = ncol(x))
            colnames(y) <- colnames(x)
            y[!dropProbes, ] <- x
            y
        })
    }

    newpe <- PBMExperiment(assays = assayl,
                           rowData = rowData(pe),
                           probeCols = pe@probeCols,
                           probeTrim = pe@probeTrim,
                           probeFilter = pe@probeFilter)

    if (verbose) {
        cat("|| - Finished probe-level model fitting.\n")
        cat("|| - Returning PBMExperiment with", nrow(newpe), "rows and", ncol(newpe), "columns.\n")
    }
    return(newpe)
}


## helper to return design matrix for specified stratify column in PBMExperiment object
.pbmExpandStratify <- function(s, strat) {

    coldat <- data.frame(colData(s), check.names = FALSE,
                         check.rows = FALSE, stringsAsFactors = FALSE)
    coldat <- tibble::rownames_to_column(coldat, "sample")

    ## check validity of grouping colData column
    if (is.null(strat)) {
        coldat <- dplyr::mutate(coldat, Stratify = "allsamples")
    } else if (! strat %in% names(colData(s))) {
        stop("Specified \"stratify=\" argument should be a column in colData.\n",
             "")
    } else {
        coldat <- dplyr::rename(coldat, Stratify = I(strat))
    }
    
    coldat <- dplyr::select(coldat, sample, Stratify)
    coldat <- dplyr::mutate(coldat, z = 1)
    coldat <- tidyr::spread(coldat, Stratify, z, fill = 0)
    
    rord <- coldat$sample
    coldat <- dplyr::select(coldat, -sample)
    coldat <- as.matrix(coldat)
    coldat[match(colnames(s), rord), , drop = FALSE]
}
