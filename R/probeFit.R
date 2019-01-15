#' Fit Probe Models
#'
#' @description
#' This function uses limma to estimate probe-level mean and moderated
#' variance estimates across replicates within conditions.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data.
#' @param assay string name of the assay to use. (default = \code{assayNames(se)[1]})
#' @param groups unquoted name of column in colData of SummarizedExperiment (or
#'        '\code{sample}') for grouping samples; values in columns will be
#'        used as group labels. If no group labels are available, all samples
#'        can be pooled by specifying \code{NULL}. (default = condition)
#' @param offset integer offset to add to intensities before log2 scaling to
#'        prevent errors with zero intensities. If set to 0, probes with
#'        zero intensities are dropped/ignored in estimation. (default = 1)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to estimating affinities. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @param ... named arguments to be passed to \code{limma::eBayes}. See
#'        Details for default parameters used in this function (different from
#'        \code{limma::eBayes} default).
#' 
#' @return
#' SummarizedExperiment object with probe-level intensity information aggregated
#' across replicates. Each assay contains probe-level information for a single
#' group in the specified \code{groups} variable.
#'
#' @md
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames assay
#' @importFrom dplyr left_join distinct
#' @importFrom limma lmFit eBayes
#' @export
#' @author Patrick Kimes
probeFit <- function(se, assay = SummarizedExperiment::assayNames(se)[1],
                     groups = condition, offset = 1L, .filter = 1L, ...) {
    
    ## define eBayes parameters with defaults
    dots <- list(...)
    eb_args <- list(trend = TRUE, robust = TRUE)
    eb_args <- replace(eb_args, names(dots), dots)

    ## filter probes
    se <- pbmFilterProbes(se, .filter) 

    ## extract probe data matrix
    datp <- SummarizedExperiment::assay(se, assay)
    datp <- log2(as.matrix(datp) + offset)
    datp[is.infinite(datp)] <- NA
    
    ## check if any probes have no non-NA values, drop from limma testing
    n <- nrow(datp)
    p <- ncol(datp)
    dropProbes <- rowSums(is.finite(datp)) == 0L
    if (any(dropProbes)) {
        datp <- datp[!dropProbes, , drop = FALSE]
    }

    ## check stratification params
    groups <- rlang::enquo(groups)
    coldat <- .pbmCheckGroups(se, groups)
        
    ## fit limma model
    fit <- limma::lmFit(datp, coldat)

    ## fit empirical Bayes adjustment with specified parameters
    eb_args$fit <- fit
    fit <- do.call(limma::eBayes, eb_args)

    ## return SummarizedExperiment with test results - one coef per assay
    alist <- list(beta = fit$coefficients,
                  sd = sweep(fit$stdev.unscaled, 1, fit$s2.post^.5, `*`),
                  df = replicate(ncol(fit$coefficients), fit$df.total))

    ## return to full dimensions
    if (any(dropProbes)) {
        alist <- lapply(alist, function(x) {
            y <- matrix(NA, nrow = n, ncol = p)
            y[!dropProbes, ] <- x
            y
        })
    }
    
    SummarizedExperiment(assays = alist, rowData = rowData(se))
}


## helper to check specified groups column in SummarizedExperiment object
.pbmCheckGroups <- function(s, grp) {
    grp_str <- rlang::quo_name(grp)

    coldat <- data.frame(colData(s), check.names = FALSE,
                         check.rows = FALSE, stringsAsFactors = FALSE)
    coldat <- tibble::rownames_to_column(coldat, "sample")

    ## check validity of grouping colData column
    if (! grp_str %in% names(colData(s))) {
        if (grp_str != "NULL") {
            stop("Specified Group must be unquoted column in colData or NULL.")
        }
        coldat <- dplyr::mutate(coldat, Groups = "allsamples")
    } else {
        coldat <- dplyr::mutate(coldat, Groups = !!grp)
    }
    
    coldat <- dplyr::select(coldat, sample, Groups)
    coldat <- dplyr::mutate(coldat, z = 1)
    coldat <- tidyr::spread(coldat, Groups, z, fill = 0)
    
    rord <- coldat$sample
    coldat <- dplyr::select(coldat, -sample)
    coldat <- as.matrix(coldat)
    coldat[match(colnames(s), rord), , drop = FALSE]
}
