#' Estimate Probe Intensities and Variances
#'
#' @description
#' This function uses limma to estimate probe-level mean and moderated
#' variance estimates across replicates within conditions.
#' 
#' @param se SummarizedExperiment object containing PBM intensity data.
#' @param assay_name string name of the assay to use. (default = "fore")
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
#' 
#' @return
#' SummarizedExperiment object with probe-level intensity information aggregated
#' across replicates. Each assay contains probe-level information for a single
#' group in the specified \code{groups} variable.
#'
#' @md
#' @import SummarizedExperiment
#' @importFrom dplyr left_join distinct
#' @importFrom limma lmFit eBayes
#' @importFrom qvalue qvalue
#' @export
#' @author Patrick Kimes
probeEstimate <- function(se, assay_name = "fore", groups = condition,
                          offset = 1L, .filter = 1L) {
    
    ## filter probes
    se <- pbmFilterProbes(se, .filter) 

    ## extract probe data matrix
    datp <- assay(se, assay_name)
    datp <- log2(as.matrix(datp) + offset)
    datp[is.infinite(datp)] <- NA

    ## check stratification params
    groups <- rlang::enquo(groups)
    coldat <- .pbmCheckGroups(se, groups)

    fit <- lmFit(datp, coldat)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)

    ## return SummarizedExperiment with test results - one coef per assay
    alist <- list(t = fit$t,
                  coefs = fit$coefficients,
                  stdev = sweep(fit$stdev.unscaled, 1, fit$s2.post^.5, `*`),
                  df = replicate(ncol(fit$coefficients), fit$df.total))

    alist <- simplify2array(alist)
    alist <- lapply(seq_len(dim(alist)[2]), function(i) alist[, i, ])
    names(alist) <- colnames(fit$coefficients)

    SummarizedExperiment(assays = alist, rowData = rowData(se))
}


#' Estimate K-mers Probe Intensities
#'
#' @description
#' After performing probe-level aggregation across samples, this function applies probe
#' set aggregation to obtain K-mer level estimates of affinity. 
#'
#' @param se probe-level affinity summaries geneated by \code{probeEstimate} with each
#'        assay corresponding to a different group label. 
#' @param assay_name string name of the assay to aggregate. If NULL, all assays are
#'        aggregated. (default = NULL)
#' @param kmers character vector of k-mers to predict.
#' @param verbose logical whether to print extra messages during model fitting
#'        procedure. (default = FALSE)
#' @param .trim interger vector of length two specifying start and end
#'        of probe sequence to be used. Default is based on the universal
#'        PBM probe design where only leading 36nt should be used. 
#'        Ignored if NULL. (default = c(1, 36))
#'
#' @return
#' SummarizedExperiment of estimated K-mer intensities.
#'
#' @importFrom dplyr select_ group_by left_join ungroup do mutate
#' @importFrom tidyr unnest spread
#' @export
#' @author Patrick Kimes
kmerEstimate <- function(se, assay_name = NULL, kmers, .filter = 1L,
                         .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {

    stopifnot(is(se, "SummarizedExperiment"))
    stopifnot(is.null(assay_name) || assay_name %in% assayNames(se))
    
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
    
    adat <- dplyr::select_(adat, .dots = c(setdiff(ovnames, "Sequence"), "coefs", "stdev", "aname"))
    adat <- dplyr::left_join(kmermap, adat, by = setdiff(ovnames, "Sequence"))
    adat <- dplyr::group_by(adat, aname, seq)
    
    adat <- dplyr::do(adat, zt = sum(.$coefs/(.$stdev^2)) / sum(.$stdev^-2))
    adat <- dplyr::ungroup(adat)
    adat <- tidyr::unnest(adat)
    adat <- tidyr::spread(adat, aname, zt)

    zdat <- dplyr::select_(adat, .dots = paste0("`", assay_name, "`"))
    zdat <- as.matrix(zdat)
    
    rdat <- dplyr::select_(adat, .dots = paste0("-`", assay_name, "`"))
    rdat <- as.data.frame(rdat, optional = TRUE)
    
    SummarizedExperiment(assays = list(z = zdat), rowData = rdat)
}


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

