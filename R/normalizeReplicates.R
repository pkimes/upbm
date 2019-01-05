## Helper function to estimate slope of qq plot between two samples.
## Slope is computed as the median ratio of median centered quantiles
## between two samples. Lower and upper quantiles are trimmed as
## speicified. By default, the lower and upper 10% are excluded from
## the computation.
qqslope <- function(x, y, lq = .1, uq = .9, center = TRUE) {
    if (is.null(x) | is.null(y)) {
        return(NA)
    }
    zz <- qqplot(x$value, y$value, plot.it = FALSE)
    nlq <- floor(length(zz$x) * lq)
    nuq <- ceiling(length(zz$x) * uq)
    zz$x <- zz$x[nlq:nuq]
    zz$y <- zz$y[nlq:nuq]
    if (center) {
        zz$x <- zz$x - median(zz$x, na.rm = TRUE)
        zz$y <- zz$y - median(zz$y, na.rm = TRUE)
    }
    median(zz$x / zz$y, na.rm = TRUE)
}

#' Cross-Replicate Normalization
#'
#' @description
#' uPBM experimental replicates can often be on different scales, even
#' after log transformation. However, these differences are often consistent
#' across samples within replicates after \code{lowertailNormalization}.
#' Based on this observation, we perform cross-replicate normalization by
#' computing log-scale scaling factors at the level of replicates.
#'
#' @param se SummarizedExperiment object containing GPR intensity information.
#' @param assay_name string name of the assay to normalize. (default = "global")
#' @param pairwise logical whether scaling factors should be computed using all
#'        pairwise comparisons or just using a single comparison against a median
#'        quantile reference. See 'Details' for more information. (default = FALSE)
#' @param onlyref logical whether scaling factors should be computed using only
#'        reference samples. (default = FALSE)
#' @param ... additional parameters to be passed to \code{qqslope} to compute
#'        scaling factors.
#'
#' @return
#' SummarizedExperiment object with normalized intensities in
#' new assay, \code{repScaled}.
#' 
#' @details
#' Scaling factors for each replicate is computed by comparing the
#' probe intensity quantiles across matched conditions in replicate
#' experiments. The replicate scaling factor is then computed as the
#' average scaling factor across all conditions. This procedure requires
#' performing all pairwise comparisons of matched conditions across all
#' pairs of replicate experiments. While rare in uPBM data, if the
#' number of replicate experiments is large (e.g. > 5), this procedure
#' can be slow. This is the default procedure (\code{pairwise = TRUE}).
#'
#' Alternatively, the scaling factor can be computed using a single mean
#' quantile reference for each condition. This reference is taken as the
#' means of the quantile (and not the quantiles of the mean) probe intensities
#' for each condition across replicates. For each replicate, the scaling factor
#' is computed as the average scaling factor across all conditions compared
#' with the corresponding reference distribution. This can be faster when
#' the number of replicates is large, as it only requires a single comparison
#' for each sample in the data set. Specifying \code{pairwise = FALSE} will
#' use this approach.
#' 
#' @export
#' @author Patrick Kimes
normalizeReplicates <- function(se, assay_name = "scaled", pairwise = FALSE,
                                onlyref = FALSE, ...) {
    
    setidy <- assay2tidy(se, assay_name, long = TRUE)
    setidy <- dplyr::select(setidy, id, condition, value)
    setidy <- dplyr::mutate(setidy, value = log2(value))

    if (pairwise) {
        setidy <- tidyr::nest(setidy, value)
        tab <- tidyr::expand(dplyr::select(setidy, id, condition),
                             nesting(id, condition),
                             nesting(id, condition))
        tab <- dplyr::filter(tab, id != id1, condition == condition1, id < id1)
        tab <- dplyr::left_join(tab, setidy, by = c("id", "condition"))
        tab <- dplyr::left_join(tab, dplyr::rename_all(setidy, funs(paste0(., "1"))),
                                by = c("id1", "condition1"))
        tab <- dplyr::mutate(tab, sfactor = mapply(qqslope, x = data, y = data1, ...))
        tab <- dplyr::select(tab, -data, -data1, -condition1)
        tab2 <- tab
        tab2[, c("id", "id1")] <- tab2[, c("id1", "id")]
        tab <- dplyr::mutate(tab, sfactor = 1/sfactor)
        tab <- bind_rows(tab, tab2)
        if (onlyref) {
            tab <- dplyr::filter(tab, grepl("-REF", condition))
        }
    
        ## compute geometric means across conditions - arithmetic mean on log-scale
        tab <- dplyr::group_by(tab, id, id1)
        tab <- dplyr::summarize(tab, sfactor = exp(mean(log(sfactor), na.rm = TRUE)))
        tab <- dplyr::ungroup(tab)

        ## average scaling across pairwise comparisons with offset to get approximate same
        ## means on log-scale
        ## e.g. 1/2 and 2 --> 3/4 and 3/2
        tab <- dplyr::group_by(tab, id)
        tab <- dplyr::mutate(tab, sfactor = (1 + sum(sfactor, na.rm = TRUE)) / (1 + sum(!is.na(sfactor))))
        tab <- dplyr::ungroup(tab)

    } else {
        ## compute per-condition quantile reference (mean)
        seref <- dplyr::group_by(setidy, condition, id)
        seref <- dplyr::summarize(seref, value = list(value))

        seref <- dplyr::mutate(seref, vsort = lapply(value, sort), vsortn = sapply(vsort, length))
        seref <- dplyr::group_by(seref, condition)
        seref <- dplyr::mutate(seref, nmin = min(vsortn, na.rm = TRUE))

        seref <- dplyr::mutate(seref, vapprox = mapply(function(x, n1, n2) { approx(1L:n1, x, n = n2)$y },
                                                       x = vsort, n1 = vsortn, n2 = nmin, SIMPLIFY = FALSE))
        seref <- dplyr::select(seref, condition, id, vapprox)
        seref <- dplyr::summarize(seref, vapprox = list(rowMeans(do.call(cbind, vapprox))))
        seref <- dplyr::mutate(seref, vapprox = lapply(vapprox, function(x) { tibble(value = x) }))

        ## compute scaling factors based on qq plot against references
        tab <- tidyr::nest(setidy, value)
        tab <- dplyr::left_join(tab, seref, by = "condition")
        tab <- dplyr::mutate(tab, sfactor = mapply(qqslope, x = vapprox, y = data, ...))
        tab <- dplyr::select(tab, -data, -vapprox)

        if (onlyref) {
            tab <- dplyr::filter(tab, grepl("-REF", condition))
        }

        ## compute geometric means across conditions - arithmetic mean on log-scale
        tab <- dplyr::group_by(tab, id)
        tab <- dplyr::summarize(tab, sfactor = exp(mean(log(sfactor), na.rm = TRUE)))
        tab <- dplyr::ungroup(tab)
    }
    
    ## first, scale replicates to put on same variability
    colData(se)$repScale <- tab$sfactor[match(colData(se)$id, tab$id)]
    assay(se, "repScaled") <- sweep(as.matrix(assay(se, assay_name)), 2,
                                    colData(se)$repScale, `^`)

    ## second, shift replicates to put in same range

    ## compute REF log2 medians for scaled data
    seriesMed <- assay2tidy(se[, grepl("-REF$", colData(se)$condition)],
                            "repScaled", long = TRUE)
    seriesMed <- group_by(seriesMed, id)
    seriesMed <- summarize(seriesMed, med = median(log2(value), na.rm = TRUE))
    seriesMed <- ungroup(seriesMed)

    ## compute mean reference medians (pre-stretching)
    sourceMed <- assay2tidy(se[, grepl("-REF$", colData(se)$condition)],
                            assay_name, long = TRUE)
    sourceMed <- group_by(sourceMed, id)
    sourceMed <- summarize(sourceMed, med = median(log2(value), na.rm = TRUE))
    sourceMed <- mean(sourceMed$med, na.rm = TRUE)

    ## scale log2 medians by average of REFs, pre-stretching
    seriesMed <- dplyr::mutate(seriesMed, scaleMed = sourceMed - med)
    
    ## add scaled assay to data - non-log scale
    colData(se)$repShift <- seriesMed$scaleMed[match(colData(se)$id, seriesMed$id)]
    assay(se, "repScaled") <- sweep(as.matrix(assay(se, "repScaled")), 2,
                                     2^colData(se)$repShift, `*`)
    se
}
