#' @title Normalize across replicates
#'
#' @description
#' Universal PBM experiments are often performed with several conditions of interest,
#' e.g. various allelic variants of the same transcription factor, assayed on arrays of
#' the same plate, with few replicate plates (commonly 2 or 3). Within individual replicates (plates,
#' observed probe intensities can vary greatly across conditions for biologically
#' uninteresting reasons, such as concentration differences. To explicitly correct for these
#' differences, normalization is performed in two steps.
#'
#' First, normalization is performed within replicates (plates). More detail on this procedure
#' can be found in the \code{\link{normalizeWithinReplicates}} documentation.
#'
#' Second, normalization is performed across replicates (plates) with the assumption that
#' biologically uninteresting differences between replicates affect probe intensities
#' both multiplicatively and additively on the log-scale. A single log-scale multiplicative normalization
#' factor is first estimated for all samples within a replicate. Then, an log-scale additive normalization
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
#'        value is guessed by looking for ``ref" in any value of the \code{stratify} column. If multiple
#'        matching values are found, an error is thrown. If the baseline condition is missing from any
#'        \code{group}, an error is thrown. (default = NULL)
#' @param pairwise a logical value whether scaling factors should be computed using all
#'        pairwise comparisons or just using a single comparison against a median
#'        quantile reference. See 'Details' for more information. (default = FALSE)
#' @param onlybaseline a logical value whether scaling factors should be computed using only
#'        baseline samples. If any group is missing a baseline sample and set to TRUE, an error is thrown.
#'        (default = FALSE)
#' @param verbose a logical value whether to print verbose output during analysis. (default = FALSE)
#' @param ... additional parameters to be passed to \code{qqslope} to compute
#'        scaling factors.
#'
#' @return
#' Original PBMExperiment object with assay containing cross-replicate normalized intensities
#' (\code{"normalized"}) and new columns added to the colData, \code{"acrossRepMultScale"} and \code{"acrossRepAddScale"},
#' containing the inverse of the log-scale multiplicative and additive scaling factors used to normalize intensities.
#' If an assay with the same name is already included in the object, it will be overwritten.
#' 
#' @details
#' The following procedure is used to estimate the log-scale multiplicative factor
#' for each replicate. First, a cross-replicate reference is computed for each condition
#' (specified in the \code{stratify} column) by taking the cross-replicate mean quantiles 
#' of the observed log2 intensities. Next, a \emph{per-sample} log multiplicative scaling factor is
#' computed by taking the median ratio of the rank-ordered and median-centered log-probe intensities
#' between the sample and the corresponding reference distribution. Visually,
#' this can be interpreted as the approximate slope of the quantile-quantile (QQ) plot generated
#' using the log-scale intensities. The \emph{per-replicate} scaling factor is then computed
#' by taking the geometric mean of the per-sample factors across all samples
#' in the replicate (specified in the \code{group} column). To reduce the impact of outlier
#' probes, per-sample scaling factors are estimated using only the middle 80% of probe intensities.
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
#' While the default behavior is to compute log-scale multiplicative factors using all conditions
#' and by making comparisons to cross-replicate reference distributions, this can be changed.
#' To use only the \code{baseline} condition, set \code{onlybaseline = TRUE}. To use all pairwise comparisons
#' between replicates (rather than comparisons against the cross-replicate references), set \code{pairwise = TRUE}.
#' Note that specifying \code{pairwise = TRUE} may be substantially more computational expensive when the number of
#' replicates is large.
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom dplyr filter left_join mutate as_tibble select summarize rename_all group_by ungroup funs
#' @importFrom tidyr nest expand gather nesting
#' @importFrom stats approx qqplot
#' @export
#' @author Patrick Kimes
normalizeAcrossReplicates <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                                      group = "id", stratify = "condition", baseline = NULL,
                                      pairwise = FALSE, onlybaseline = FALSE, verbose = FALSE, ...) {

    stopifnot(is(pe, "PBMExperiment"))
    stopifnot(assay %in% SummarizedExperiment::assayNames(pe))
    stopifnot(is(pairwise, "logical"))
    stopifnot(is(onlybaseline, "logical"))
    
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
                                needbl = onlybaseline, verb = verbose)
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
    petidy <- tidyr::gather(petidy, sample, value)
    petidy <- dplyr::left_join(petidy, coldat, by = "sample")
    petidy <- dplyr::mutate(petidy, value = log2(value))
    
    ## compute median log2 intensities before any scaling
    pemed <- dplyr::group_by(petidy, Group, Stratify, sample, isBaseline)
    pemed <- dplyr::summarize(pemed, med = median(value, na.rm = TRUE))
    pemed <- dplyr::group_by(pemed, Group)
    pemed <- dplyr::left_join(pemed, dplyr::select(dplyr::filter(pemed, isBaseline), Group, med),
                              by = "Group", suffix = c("", "_bl"))
    
    ## compute log-scale multiplicative scaling factors
    if (pairwise) {
        if (verbose) {
            cat("|| - Performing cross-replicate normalization using all pairwise comparisons (pairwise = TRUE).\n")
        }
        petidy <- tidyr::nest(petidy, value)
        tab <- tidyr::expand(dplyr::select(petidy, Group, Stratify, sample),
                             tidyr::nesting(Group, Stratify, sample), tidyr::nesting(Group, Stratify, sample))
        
        tab <- dplyr::filter(tab, Group != Group1, Group < Group1)
        tab <- dplyr::filter(tab, Stratify == Stratify1)

        tab <- dplyr::left_join(tab, petidy, by = c("Group", "Stratify", "sample"))
        tab <- dplyr::left_join(tab, dplyr::rename_all(petidy, dplyr::funs(paste0(., "1"))),
                                by = c("Group1", "Stratify1", "sample1"))
        tab <- dplyr::mutate(tab, sfactor = mapply(qqslope, x = data, y = data1, ...))
        tab <- dplyr::select(tab, -data, -data1, -Stratify1, -sample1)
        
        tab2 <- tab
        tab2[, c("Group", "Group1")] <- tab2[, c("Group1", "Group")]
        tab <- dplyr::mutate(tab, sfactor = 1/sfactor)
        tab <- bind_rows(tab, tab2)
        
        if (onlybaseline) {
            if (verbose) {
                cat("|| - Performing cross-replicate normalization using only baseline condition (onlybaseline = TRUE).\n")
            }
            tab <- dplyr::filter(tab, isBaseline, isBaseline1)
        } else if (verbose) {
            cat("|| - Performing cross-replicate normalization using all conditions (onlybaseline = FALSE).\n")
        }

        ## compute geometric means across conditions - arithmetic mean on log-scale
        tab <- dplyr::group_by(tab, Group, Group1)
        tab <- dplyr::summarize(tab, sfactor = exp(mean(log(sfactor), na.rm = TRUE)))
        tab <- dplyr::ungroup(tab)

        ## average scaling across pairwise comparisons with offset to get approximate same
        ## means on log-scale
        ## e.g. 1/2 and 2 --> 3/4 and 3/2
        tab <- dplyr::group_by(tab, Group)
        tab <- dplyr::mutate(tab, sfactor = (1 + sum(sfactor, na.rm = TRUE)) / (1 + sum(!is.na(sfactor))))
        tab <- dplyr::ungroup(tab)

    } else {
        if (verbose) {
            cat("|| - Performing cross-replicate normalization using mean reference comparisons (pairwise = FALSE).\n")
        }
        ## compute per-condition quantile reference (mean)
        peref <- dplyr::group_by(petidy, Stratify, Group, sample)
        peref <- dplyr::summarize(peref, value = list(value))
        
        ## filter out conditions with only 1 replicate
        peref <- dplyr::group_by(peref, Stratify)
        peref <- dplyr::mutate(peref, nreps = dplyr::n_distinct(Group))
        peref <- dplyr::ungroup(peref)
        if (any(peref$nreps == 1L) && verbose) {
            cat("|| - Dropping following conditions with replicates in only one group from calculation of normalization factor:\n")
            cat("||     -", paste0(peref$Stratify[peref$nreps == 1L], collapse = ", "), "\n")
        }
        peref <- dplyr::filter(peref, nreps > 1L)
        peref <- dplyr::select(peref, -nreps)
        
        peref <- dplyr::mutate(peref, vsort = lapply(value, sort), vsortn = sapply(vsort, length))
        peref <- dplyr::group_by(peref, Stratify)
        peref <- dplyr::mutate(peref, nmin = min(vsortn, na.rm = TRUE))

        peref <- dplyr::mutate(peref, vapprox = mapply(function(x, n1, n2) { stats::approx(1L:n1, x, n = n2)$y },
                                                       x = vsort, n1 = vsortn, n2 = nmin, SIMPLIFY = FALSE))

        peref <- dplyr::select(peref, Stratify, Group, sample, vapprox)
        peref <- dplyr::summarize(peref, vapprox = list(rowMeans(do.call(cbind, vapprox))))
        peref <- dplyr::mutate(peref, vapprox = lapply(vapprox, function(x) { tibble(value = x) }))

        ## compute scaling factors based on qq plot against references
        tab <- tidyr::nest(petidy, value)
        tab <- dplyr::left_join(tab, peref, by = "Stratify")

        ## -- filter out conditions with only 1 replicate
        tab <- dplyr::filter(tab, vapply(vapprox, function(x) { !is.null(x) }, logical(1L))) 

        if (onlybaseline) {
            if (verbose) {
                cat("|| - Performing cross-replicate normalization using only baseline condition (onlybaseline = TRUE).\n")
            }
            tab <- dplyr::filter(tab, isBaseline)
        } else if (verbose) {
            cat("|| - Performing cross-replicate normalization using all conditions (onlybaseline = FALSE).\n")
        }

        ## (do the actual computing)
        tab <- dplyr::mutate(tab, sfactor = mapply(qqslope, x = vapprox, y = data, ...))
        tab <- dplyr::select(tab, -data, -vapprox)
    
        ## compute geometric means across conditions - arithmetic mean on log-scale
        tab <- dplyr::group_by(tab, Group)
        tab <- dplyr::summarize(tab, sfactor = exp(mean(log(sfactor), na.rm = TRUE)))
        tab <- dplyr::ungroup(tab)
    }

    ## compute log-scale additive scaling factors
    tab <- dplyr::left_join(tab, pemed, by = "Group")
    
    med_mean <- mean(dplyr::filter(tab, isBaseline)$med)
    tab <- dplyr::mutate(tab, afactor = med_mean - sfactor * med_bl)
    tab <- dplyr::mutate(tab, med_new = (med - med_bl) * sfactor + med_mean)

    ## fill in additive factors for any groups missing baseline using
    ## normalization factors for scans with same condition from other groups
    if (any(is.na(tab$med_bl))) {
        tab <- dplyr::group_by(tab, Stratify)
        tab <- dplyr::mutate(tab, afactor = ifelse(is.na(afactor),
                                                   mean(med_new, na.rm = TRUE) - sfactor * med,
                                                   afactor))
    }
    tab <- dplyr::group_by(tab, Group)
    tab <- dplyr::summarize(tab, sfactor = mean(sfactor), afactor = mean(afactor))
    tab <- dplyr::ungroup(tab)

    if (any(is.na(tab$afactor) | is.na(tab$sfactor))) {
        warning("Cross-replicate normalization factors could not be estimated for some groups.\n",
                "This is likely due to the group not containing a baseline scan \n",
                "and not sharing any conditions with any other groups containing a baseline scan.\n",
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

