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
#'        per k-mer depends strongly on k. (default = 10L)
#' @param verbose logical whether to print extra messages during model fitting
#'        procedure. (default = FALSE)
#' @param log_scale logical whether to calculate bias on log2 scale
#'        intensities. (default = FALSE)
#' @param offset integer offset to add to intensities before performing
#'        log2 scale calculations. (default = 1L)
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
approxPositionWeights <- function(se, kmers, nk = 10L, verbose = FALSE, .filter = 1L,
                                  log_scale = FALSE, offset = 1L,
                                  .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {
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
        kmall <- dplyr::mutate(kmall, log2ratio = log2( log2(value.pr + offset) /
                                                        log2(value.km + offset) ))
    } else {
        kmall <- dplyr::mutate(kmall, log2ratio = log2(value.pr / value.km))
    }
    
    ptrend <- dplyr::select(kmall, sample, pos, log2ratio)
    ptrend <- dplyr::group_by(ptrend, sample, pos)
    ptrend <- dplyr::summarize(ptrend, log2ratio = median(log2ratio, na.rm = TRUE))
    ptrend <- dplyr::ungroup(ptrend)
    ptrend <- tidyr::spread(ptrend, sample, log2ratio)
    dplyr::arrange(ptrend, pos)
}
