#' Approximate Positional Bias
#'
#' Probe intensities in protein binding microarrays (PBMs) are
#' highly dependent not only on the sequence content of each
#' probe, but also the position of the contained sequences within
#' the probe. Typically, if a motif is closer to the glass slide,
#' the binding ability of the TF will be reduced, and the resulting
#' probe intensity will be lower than if the same motif were
#' positioned within the probe closer to the free end of the
#' sequence. This function provides a simple way to estimate
#' the positional bias of a TF based on the relative probe intensities
#' for probes containing the top N K-mers, ranked by median K-mer
#' affinity.
#' 
#' @param se SummarizedExperiment of probe intensities
#' @param assay_name string name of the assay to use. (default = "fore")
#' @param kmers character vector of k-mers to predict.
#' @param nk number of top k-mers to use to estimate positional
#'        trend for each sample. This should be larger for longer
#'        k-mers and smaller for shorter k-mers as the number of probes
#'        per k-mer depends strongly on k. (default = 100L)
#' @param smooth logical whether to return loess smoothed bias esimates along
#'        with median bias estimates. If \code{withse} is also TRUE, each
#'        position is inversely weighted by the estimated standard error.
#'        If TRUE, will return list of tables, "median" and "smooth".
#'        Loess smoothing can be adjusted using the \code{.smooth.span}
#'        parameter. (default = FALSE)
#' @param withse logical whether to return standard error estimates along
#'        with median trend estimates. If TRUE, will return list of tables,
#'        "median" and "se". The standard error is estimated by the scaled MAD
#'        divided by the square root of the number of observations at each
#'        position. (default = FALSE)
#' @param verbose logical whether to print extra messages during model fitting
#'        procedure. (default = FALSE)
#' @param log_scale logical whether to calculate bias on log2 scale
#'        intensities. (default = FALSE)
#' @param offset integer offset to add to intensities before performing
#'        log2 scale calculations. (default = 1L)
#' @param .smooth.span numeric span value to be passed to \code{stats::smooth}
#'        to adjust amount of smoothing. The default value is set to a smaller
#'        value than the default \code{stats::smooth} value. (default = 1/2) 
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to estimating weights. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @param .trim interger vector of length two specifying start and end
#'        of probe sequence to be used. Default is based on the universal
#'        PBM probe design where only leading 36nt should be used.
#'        Ignored if NULL. (default = c(1, 36))
#'
#' @return
#' A vector of weights of length equal to the M - K + 1, where M is the length
#' of the trimmed probe sequence, and K is the length of the K-mer sequences.
#' 
#' @export
#' @importFrom dplyr rename as_tibble top_n group_by ungroup mutate left_join select summarize arrange
#' @importFrom tidyr gather spread
#' @author Patrick Kimes
approxPositionWeights <- function(se, assay_name = "fore", kmers, nk = 100L, smooth = TRUE, withse = TRUE, 
                                  verbose = FALSE, log_scale = FALSE, offset = 1L, .smooth.span = 1/2,
                                  .filter = 1L, .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {
    ## check kmers specified
    kmers <- checkKmers(kmers, verbose)

    ## check Sequence info in rowData
    se <- checkProbeSequences(se, verbose)

    ## filter probes
    se <- pbmFilterProbes(se, .filter)

    ## trim probe sequences
    se <- trimProbeSequences(se, .trim)

    ## create mapping between kmers and probes
    kmmap <- mapKmers(rowData(se), kmers)
    kmmap <- dplyr::rename(kmmap, kmer = seq)

    ## compute median k-mer intensities
    kmsum <- summarizeKmers(se, assay_name = assay_name, kmers = kmers,
                            stat_set = "median", .trim = .trim, .filter = .filter)
    kmsum <- cbind(rowData(kmsum), assay(kmsum, "medianIntensity"))
    kmsum <- dplyr::as_tibble(as.data.frame(kmsum, optional = TRUE))

    ## determine top nk k-mers for each sample
    topkm <- tidyr::gather(kmsum, sample, value, one_of(colnames(se)))
    topkm <- dplyr::group_by(topkm, sample)
    topkm <- dplyr::top_n(topkm, nk, value)
    topkm <- dplyr::ungroup(topkm)

    ## extract probe-level intensities
    sevals <- dplyr::as_tibble(as.data.frame(assay(se, assay_name), optional = TRUE))
    sevals <- dplyr::mutate(sevals, probe_idx = 1:n())
    sevals <- tidyr::gather(sevals, sample, value, one_of(colnames(se)))

    ## get log2 ratios relative to median k-mer intensity
    kmall <- dplyr::left_join(topkm, kmmap, by = "kmer")
    kmall <- dplyr::left_join(kmall, sevals, by = c("probe_idx", "sample"), suffix = c(".km", ".pr"))

    if (log_scale) {
        kmall <- dplyr::mutate(kmall, l2r = log2( log2(value.pr + offset) /
                                                  log2(value.km + offset) ))
    } else {
        kmall <- dplyr::mutate(kmall, l2r = log2(value.pr / value.km))
    }
    
    kmall <- dplyr::select(kmall, sample, pos, l2r)

    ## compute position bias and SE estimates 
    ptrend <- dplyr::group_by(kmall, sample, pos)
    if (withse) {
        ptrend <- dplyr::summarize(ptrend, l2r.median = median(l2r, na.rm = TRUE),
                                   l2r.mad = mad(l2r, na.rm = TRUE) / sqrt(n()))
    } else {
        ptrend <- dplyr::summarize(ptrend, l2r.median = median(l2r, na.rm = TRUE),
                                   l2r.mad = 1L)
    }
    ptrend <- dplyr::ungroup(ptrend)

    ## compute smoothed position bias estimates
    if (smooth) {
        ptrend <- dplyr::group_by(ptrend, sample)
        ptrend <- dplyr::do(ptrend,
                            l2r.smooth = predict(loess(l2r.median ~ pos, data = .,
                                                       weights = 1 / .$l2r.mad,
                                                       span = .smooth.span),
                                                 newdata = .$pos),
                            l2r.median = .$l2r.median,
                            l2r.mad = .$l2r.mad,
                            pos = .$pos)
        ptrend <- tidyr::unnest(ptrend)
    }

    ## clean up tables and return
    pmed <- dplyr::select(ptrend, sample, pos, l2r.median)
    pmed <- tidyr::spread(pmed, sample, l2r.median)
    pmed <- dplyr::arrange(pmed, pos)
    res <- list(median = pmed)
    if (withse) {
        pse <- dplyr::select(ptrend, sample, pos, l2r.mad)
        pse <- tidyr::spread(pse, sample, l2r.mad)
        pse <- dplyr::arrange(pse, pos)
        res <- c(res, list(stderr = pse))
    }
    if (smooth) {
        psm <- dplyr::select(ptrend, sample, pos, l2r.smooth)
        psm <- tidyr::spread(psm, sample, l2r.smooth)
        psm <- dplyr::arrange(psm, pos)
        res <- c(res, list(smooth = psm))
    }
    
    if (length(res) == 1L) {
        return(res[[1]])
    } else {
        return(res)
    }
}


#' Simple Function to Plot Position Weights
#'
#' @description
#' The \code{approxPositionWeights} function returns either a tibble or
#' list of tibbles contianing estimated approximate probe sequence
#' position biases. These are biases are expressed as the log2 ratio of
#' the observed to median probe intensities of fixed K-mer sequences
#' along the probe sequence. Typically, the bias is characterizerd by a
#' sharp drop in affinity near the end of the probe and a gradual decline
#' as the sequence is positioned closer to the glass slide, with the
#' highest affinity between 5 and 10 nucleotides from the end of the probe.
#'
#' This function is a simple wrapper to plot the log2 ratio position bias
#' for samples as a function of the "position" from the end of the probe,
#' with 1 corresponding to the end of the probe, and moving towards the
#' glass slide as values increase.
#' 
#' @param weights tibble of position weights estimated, e.g. using
#'        \code{approxPositionWeights}.
#' @param se optional specification of SummarizedExperiment used to
#'        generate weights. This is only useful if a column in the
#'        colData of the SummarizedExperiment should be used to label
#'        samples rather than the default column names (samples). If
#'        a different colData column should be used instead, this should
#'        be specified to \code{color} unquoted. (default = NULL)
#' @param color column in \code{colData(se)} that should be used to
#'        group and color samples in the plot. This should be specified
#'        unquoted. By default, the column names are labeled "sample"
#'        and used. (default = sample)
#'
#' @return
#' ggplot object of position bias
#'
#' @examples
#' \dontrun{
#' ## estimate approximate weights 
#' we <- approxPositionWeights(se)
#'
#' ## plot smoothed weights
#' pbmPlotPositionWeights(we$smooth)
#' }
#' 
#' @importFrom dplyr left_join as_tibble
#' @importFrom rlang enquo
#' @importFrom tidyr gather
#' @import SummarizedExperiment
#' @export
#' @author Patrick Kimes
pbmPlotPositionWeights <- function(weights, se = NULL, color = sample) {
    color <- rlang::enquo(color)

    stopifnot(is(weights, "data.frame"))
    stopifnot("pos" %in% names(weights))
    
    wedat <- tidyr::gather(weights, sample, value, -pos)

    if (!is.null(se)) {
        stopifnot(is(se, "SummarizedExperiment"))
        stopifnot(ncol(se) >= ncol(weights) - 1L) 
        coldat <- as.data.frame(colData(gpr), optional = TRUE)
        coldat <- as_tibble(coldat, rownames = "sample")
        wedat <- dplyr::left_join(wedat, coldat, by = "sample")
    }
    
    ggplot(wedat, aes(x = pos, y = value, color = !!color)) +
        geom_line(alpha = 1/2) +
        geom_point() +
        theme_bw() +
        xlab("position (nt)") +
        ylab("bias (log2 ratio)")
}
