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
approxPositionWeights <- function(se, kmers, nk = 100L, smooth = TRUE, withse = TRUE,
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
    kmsum <- summarizeKmers(se, kmers = kmers, stat_set = "median",
                            .trim = .trim, .filter = .filter)
    kmsum <- cbind(rowData(kmsum), assay(kmsum, "medianIntensity"))
    kmsum <- dplyr::as_tibble(as.data.frame(kmsum, optional = TRUE))

    ## determine top nk k-mers for each sample
    topkm <- tidyr::gather(kmsum, sample, value, one_of(colnames(se)))
    topkm <- dplyr::group_by(topkm, sample)
    topkm <- dplyr::top_n(topkm, nk, value)
    topkm <- dplyr::ungroup(topkm)

    ## extract probe-level intensities
    sevals <- dplyr::as_tibble(as.data.frame(assay(se), optional = TRUE))
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
