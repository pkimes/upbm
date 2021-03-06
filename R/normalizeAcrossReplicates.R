#' @title Normalize across replicates
#'
#' @description
#' Universal PBM experiments are often performed with several conditions of interest,
#' e.g. allelic variants, assayed on separate arrays of the same plate with few replicates.
#' Within and across plates, probe intensities can vary for biologically
#' uninteresting reasons, such as concentration differences. To explicitly correct for these
#' differences, normalization is performed in two steps.
#'
#' First, normalization is performed within replicates (plates). More detail on this procedure
#' can be found in the \code{\link{normalizeWithinReplicates}} documentation.
#'
#' Second, normalization is performed across replicates (plates) with the assumption that
#' biologically uninteresting differences between replicates affect probe intensities
#' both multiplicatively and additively on the log-scale. A single log-scale multiplicative normalization
#' factor is first estimated for all samples within a replicate. Then, a log-scale additive normalization
#' is estimated such that the median intensities of the \code{baseline} samples in each replicate
#' are equal. More details on this calculation are provided below.
#'
#' @param pe SummarizedExperiment object containing GPR intensity information.
#' @param assay a string name of the assay to normalize.
#'        (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param group a character string specifying a column in \code{colData(pe)} to use for grouping replicates.
#'        (default = \code{"id"})
#' @param stratify a character string specifying a column in \code{colData(pe)} to use for determining
#'        the unique baseline scan within each \code{group} and to match samples across values of \code{group}.
#'        (default = \code{"condition"})
#' @param baseline a character string specifying the baseline condition in the \code{stratify} column to normalize
#'        other conditions against within each \code{group}. If not specified and set to NULL, the baseline
#'        value is guessed by looking for values in the \code{stratify} column ending in ``ref". If multiple
#'        unique matching values are found, a warning is thrown and the first matching sample is used.
#'        (default = NULL)
#' @param verbose a logical value whether to print verbose output during analysis. (default = FALSE)
#'
#' @return
#' Original PBMExperiment object with assay containing cross-replicate normalized intensities
#' (\code{"normalized"}) and new columns added to the colData, \code{"acrossRepMultScale"} and \code{"acrossRepAddScale"},
#' containing the inverse of the log-scale multiplicative and additive scaling factors used to normalize intensities.
#' If an assay with the same name is already included in the object, it will be overwritten.
#' 
#' @details
#' The following procedure is used to estimate the log-scale multiplicative factor
#' for each replicate. First, a cross-replicate reference is computed for each baseline condition
#' (specified by \code{stratify=} and \code{baseline=}) by taking the cross-replicate mean quantiles 
#' of the observed log2 intensities. Next, a \emph{per-replicate} log multiplicative scaling factor is
#' computed by taking the median ratio of the rank-ordered and median-centered log-probe intensities
#' between the baseline samples in each replicate and the reference distribution. Visually,
#' this can be interpreted as the approximate slope of the quantile-quantile (QQ) plot generated
#' using log-scale intensities. To reduce the impact of outlier probes, scaling factors are estimated
#' using only the middle 80% of probe intensities.
#'
#' After log-scale multiplicative factors have been estimated to correct for differences in
#' log-scale variance across replicates, a second log-scale additive factor is estimated
#' for each replicate to correct for differences in log-scale shift. A "global median" intensity
#' is first calculated across replicates by taking the geometric mean of the median
#' intensities in all \code{baseline} samples across replicates. This "global median" is computed
#' using the input probe intensities, i.e. without any cross-replicate normalization. 
#' The log-scale additive factor estimated as the difference between the median normalized probe
#' intensity of the \code{baseline} sample in each replicate and the "global median".
#' While the log-scale additive factor is estimated using only \code{baseline} samples, the normalization
#' is applied to all samples in the replicate.
#'
#' Cross-replicate normalization is first carried out for replicates containing a baseline sample as
#' described above. Replicates without a baseline sample are then normalized to already normalized
#' replicates using overlapping conditions in the \code{stratify=} column.
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom dplyr filter left_join mutate as_tibble select summarize rename_all group_by ungroup funs desc
#' @importFrom tidyr nest expand pivot_longer nesting
#' @importFrom stats approx qqplot
#' @importFrom tidyselect everything
#' @export
#' @author Patrick Kimes
normalizeAcrossReplicates <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                                      group = "id", stratify = "condition", baseline = NULL,
                                      verbose = FALSE) {

    stopifnot(is(pe, "PBMExperiment"))
    stopifnot(assay %in% SummarizedExperiment::assayNames(pe))
    
    ## check normalization groups
    stopifnot(group %in% names(colData(pe)))
    stopifnot(stratify %in% names(colData(pe)))
    
    if (verbose) {
        cat("|| upbm::normalizeAcrossReplicates \n")
        cat("|| - Starting cross-replicate normalization for", ncol(pe), "PBM scans.\n")
    }
    
    if (verbose) {
        cat("|| - Filtering probes according to", length(pe@probeFilter),
            "probeFilter rule(s).\n")
        ntotal <- nrow(pe)
    }

    ## filter probes - only for computing shift/scale factors (return original pe)
    fpe <- pbmFilterProbes(pe)
    
    if (verbose) {
        cat("|| - Data filtered from", ntotal, "probes to", nrow(fpe), "probes.\n")
    }

    ## check stratify params
    strats <- .pbmCheckStratify(s = fpe, strat = stratify, bl = baseline, gp = group,
                                needbl = FALSE, verb = verbose)
    coldat <- strats$coldat
    baseline <- strats$baseline

    ## check cases with only one replicate
    if (length(unique(coldat$Group)) == 1L) {
        if (! "normalized" %in% assayNames(pe)) {
            SummarizedExperiment::assay(pe, "normalized") <- SummarizedExperiment::assay(pe, assay)
        }
        if (verbose) {
            cat("|| - Sample only has one replicate group - no normalization needed.\n")
            cat("|| - Finished cross-replicate normalization.\n")
            cat("|| - Returning PBMExperiment with", nrow(pe), "rows and", ncol(pe), "columns.\n")
        }
        return(pe)
    }
    
    if (verbose) {
        cat("|| - Performing cross-replicate normalization with:\n")
        cat("||     - replicate group column:", group, "\n")
        cat("||     - condition column:", stratify, "\n")
        cat("||     - baseline condition:", baseline, "\n")
    }
    
    petidy <- as.data.frame(assay(pe, assay), optional = TRUE)
    petidy <- dplyr::as_tibble(petidy)
    petidy <- tidyr::pivot_longer(petidy, names_to = "sample", values_to = "value",
                                  tidyselect::everything())
    petidy <- dplyr::left_join(petidy, coldat, by = "sample")
    petidy <- dplyr::mutate(petidy, value = log2(value))
    
    ## compute median log2 intensities before any scaling
    pemed <- dplyr::group_by(petidy, Group, Stratify, sample, isBaseline)
    pemed <- dplyr::summarize(pemed, med = median(value, na.rm = TRUE))
    pemed <- dplyr::ungroup(pemed)
    pemed <- dplyr::left_join(pemed, dplyr::select(dplyr::filter(pemed, isBaseline), Group, med),
                              by = "Group", suffix = c("", "_bl"))

    ## compute per-condition quantile reference (mean)
    blref <- dplyr::filter(petidy, isBaseline)
    blref <- dplyr::group_by(blref, Group, sample)
    blref <- dplyr::summarize(blref, value = list(value))
    blref <- dplyr::ungroup(blref)
    
    blref <- dplyr::mutate(blref, vsort = lapply(value, sort), vsortn = sapply(vsort, length))
    blref <- dplyr::mutate(blref, nmin = min(vsortn, na.rm = TRUE))

    blref <- dplyr::mutate(blref, vapprox = mapply(function(x, n1, n2) { stats::approx(1L:n1, x, n = n2)$y },
                                                   x = vsort, n1 = vsortn, n2 = nmin, SIMPLIFY = FALSE))

    blref <- dplyr::select(blref, Group, sample, vapprox)
    blref <- dplyr::summarize(blref, vapprox = list(rowMeans(do.call(cbind, vapprox))))
    blref <- dplyr::mutate(blref, vapprox = lapply(vapprox, function(x) { tibble(value = x) }))
    blref <- dplyr::mutate(blref, Stratify = !!baseline)
    
    ## compute scaling factors based on qq plot against references
    tab <- dplyr::filter(petidy, isBaseline)
    tab <- tidyr::nest(tab, data = c(value))
    tab <- dplyr::left_join(tab, blref, by = "Stratify")
    
    ## (do the actual computing)
    tab <- dplyr::mutate(tab, sfactor = mapply(qqslope, x = vapprox, y = data))
    tab <- dplyr::select(tab, Group, sfactor)
    
    ## compute log-scale additive scaling factors
    tab <- dplyr::left_join(pemed, tab, by = "Group")
    med_mean <- mean(dplyr::filter(tab, isBaseline)$med)
    tab <- dplyr::mutate(tab, afactor = med_mean - sfactor * med_bl)
    tab <- dplyr::mutate(tab, med_new = (med - med_bl) * sfactor + med_mean)

    ## determine groups with or without baseline condition sample
    tab <- dplyr::group_by(tab, Group)
    tab <- dplyr::mutate(tab, hasBaseline = any(isBaseline))
    tab <- dplyr::ungroup(tab)
    
    ## loop through additional secondary baseline conditions for replicates
    ## without baseline conditions
    if (any(is.na(tab$sfactor))) {
        
        ## maximum number times to try to solve for scaling factors
        ## -- effectively, maximum number of steps from baseline-containing
        ##    replicate before giving up on normalizing replicate
        max_try <- 10
        for (i in seq_len(max_try)) { 
            altBaselines <- unique(tab$Stratify[is.na(tab$sfactor)])

            abltab <- dplyr::filter(tab, Stratify %in% altBaselines)

            ## choose order of alternative baseline conditions to try based
            altBaselines <- dplyr::group_by(abltab, Stratify)
            altBaselines <- dplyr::summarize(altBaselines, nbl = sum(hasBaseline), n = n())
            altBaselines <- dplyr::arrange(altBaselines, dplyr::desc(nbl), dplyr::desc(n))
            altBaselines <- altBaselines$Stratify

            ablref <- dplyr::select(petidy, sample, value)
            ablref <- dplyr::left_join(ablref, tab, by = "sample")

            iblref <- dplyr::filter(ablref, !is.na(sfactor), Stratify %in% altBaselines)

            ## skip if no normalized samples available for condition
            if (nrow(iblref) == 0L) {
                break
            }

            ## use normalized intensities
            iblref <- dplyr::mutate(iblref, value = afactor + value * sfactor)
            
            ## compute quantile reference for alternative baseline (mean)
            iblref <- dplyr::group_by(iblref, Group, Stratify, sample)
            iblref <- dplyr::summarize(iblref, value = list(value))
            iblref <- dplyr::ungroup(iblref)
            
            iblref <- dplyr::mutate(iblref, vsort = lapply(value, sort), vsortn = sapply(vsort, length))
            iblref <- dplyr::group_by(iblref, Stratify)
            iblref <- dplyr::mutate(iblref, nmin = min(vsortn, na.rm = TRUE))
            iblref <- dplyr::ungroup(iblref)
            
            iblref <- dplyr::mutate(iblref, vapprox = mapply(function(x, n1, n2) {
                stats::approx(1L:n1, x, n = n2)$y },
                x = vsort, n1 = vsortn, n2 = nmin, SIMPLIFY = FALSE))

            iblref <- dplyr::select(iblref, Group, Stratify, sample, vapprox)
            iblref <- dplyr::group_by(iblref, Stratify)
            iblref <- dplyr::summarize(iblref, vapprox = list(rowMeans(do.call(cbind, vapprox))))
            iblref <- dplyr::mutate(iblref, vapprox = lapply(vapprox, function(x) { tibble(value = x) }))
            ##iblref <- dplyr::mutate(iblref, Stratify = !!ibl)
            
            ## compute scaling factors based on qq plot against references
            itab <- dplyr::filter(ablref, is.na(sfactor))
            itab <- tidyr::nest(itab, data = c(value))
            itab <- dplyr::left_join(itab, iblref, by = "Stratify")

            ## (do the actual computing)
            itab <- dplyr::mutate(itab, sfactor = mapply(qqslope, x = vapprox, y = data))

            ## take geometric mean across conditions for each group
            itab <- dplyr::group_by(itab, Group)
            itab <- dplyr::summarize(itab, sfactor = 2^mean(log2(sfactor)))
            
            ## merge with old 
            tab <- dplyr::left_join(tab, itab, by = "Group", suffix = c("", ".i"))

            ## determine additive scaling factors
            imed_mean <- dplyr::filter(tab, Stratify %in% altBaselines, !is.na(sfactor))
            imed_mean <- dplyr::group_by(imed_mean, Stratify)
            imed_mean <- dplyr::summarize(imed_mean, med.i = mean(med_new))

            tab <- dplyr::left_join(tab, imed_mean, by = "Stratify")
            tab <- dplyr::mutate(tab, afactor.i = med.i - sfactor.i * med)

            ## replace NAs with new values
            tab <- dplyr::mutate(tab,
                                 sfactor = ifelse(is.na(sfactor), sfactor.i, sfactor),
                                 afactor = ifelse(is.na(afactor), afactor.i, afactor))
            tab <- dplyr::select(tab, -sfactor.i, -afactor.i, -med.i)

            ## stop if we don't need to try more times
            if (!any(is.na(tab$sfactor)))
                break
        }
    }
    
    if (any(is.na(tab$afactor) | is.na(tab$sfactor))) {
        warning("Cross-replicate normalization factors could not be estimated for some groups.\n",
                "This is likely due to the group not containing a baseline scan \n",
                "and not sharing any conditions with any other groups containing a baseline scan (within", max_try, "steps).\n",
                "In this case, the scans in the group cannot be normalized relative to other groups.\n",
                "Scans in following groups will have NA values: ",
                paste0(tab$Group[is.na(tab$afactor) | is.na(tab$sfactor)], collapse = ", "), ".")
    }
    
    ## perform normalization on input assay
    new_assay <- sweep(as.matrix(SummarizedExperiment::assay(pe, assay)), 2,
                       tab$sfactor[match(colData(pe)[[group]], tab$Group)], `^`)
    new_assay <- sweep(new_assay, 2,
                       2^tab$afactor[match(colData(pe)[[group]], tab$Group)], `*`)
    
    ## add to input PBMExperiment
    if ("normalized" %in% assayNames(pe)) {
        SummarizedExperiment::assay(pe, "normalized") <- NULL
        if (verbose) {
            cat("|| - Original PBMExperiment object contains \"normalized\" assay.\n")
            cat("|| - Existing \"normalized\" assay will be overwritten.\n")
        }
    }
    SummarizedExperiment::assays(pe) <- c(S4Vectors::SimpleList(normalized = new_assay),
                                          SummarizedExperiment::assays(pe))

    if (verbose) {
        cat("|| - Adding cross-replicate normalization factors to colData as \"acrossRepMultScale\" and \"acrossRepAddScale\".\n")
    }

    ## add scaling factors to colData
    colData(pe)$acrossRepMultScale <- tab$sfactor[match(colData(pe)[[group]], tab$Group)]
    colData(pe)$acrossRepAddScale <- tab$afactor[match(colData(pe)[[group]], tab$Group)]
    
    if (verbose) {
        cat("|| - Finished cross-replicate normalization.\n")
        cat("|| - Returning PBMExperiment with", nrow(pe), "rows and", ncol(pe), "columns.\n")
    }

    return(pe)
}

## Helper function to estimate slope of qq plot between two samples.
## Slope is computed as the median ratio of median centered quantiles
## between two samples. Lower and upper quantiles are trimmed as
## speicified. By default, the lower and upper 10% are excluded from
## the computation.
qqslope <- function(x, y, lq = .1, uq = .9, center = TRUE) {
    if (is.null(x) | is.null(y)) {
        return(NA)
    }
    zz <- stats::qqplot(x$value, y$value, plot.it = FALSE)
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
