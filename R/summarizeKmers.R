#' Summarize K-mer Intensities
#'
#' This function calculates k-mer intensity information from
#' PBM data. 
#' 
#' @param se SummarizedExperiment object containing PBM intensity data.
#' @param assay_name string name of the assay to use. (default = "fore")
#' @param kmers character vector of k-mers to predict.
#' @param offset integer offset to add to intensities before log2 scaling to
#'        prevent errors with zero intensities. (default = 1)
#' @param verbose logical whether to print extra messages during model fitting
#'        procedure. (default = FALSE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to estimating affinities. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @param .trim interger vector of length two specifying start and end
#'        of probe sequence to be used. Default is based on the universal
#'        PBM probe design where only leading 36nt should be used. 
#'        Ignored if NULL. (default = c(1, 36))
#' 
#' @return
#' SummarizedExperiment object with summarized intensity information in
#' separate assays.
#' 
#' @md
#' @importFrom dplyr as_tibble mutate left_join select
#' @importFrom tidyr nest
#' @importFrom matrixStats colMedians colMeans2 colSds colMads
#' @export
#' @author Patrick Kimes
summarizeKmers <- function(se, assay_name = "fore", kmers = NULL, offset = 1,
                           verbose = FALSE, .filter = 1L,
                           .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {

    ## check kmers specified
    kmers <- checkKmers(kmers, verbose)

    ## check Sequence info in rowData
    se <- checkProbeSequences(se, verbose)
    
    ## filter probes
    se <- pbmFilterProbes(se, .filter) 

    ## trim probe sequences
    se <- trimProbeSequences(se, .trim)

    ## find mapping between kmers and probes
    ovnames <- intersect(names(rowData(se)), c("Row", "Column", "ID", "Sequence"))
    kmermap <- mapKmers(rowData(se)[, ovnames], kmers)
    
    ## use ordering from input 'kmers'
    kmermap$seq <- factor(kmermap$seq, levels = kmers)
    
    ## extract intensities
    pdat <- assay(se, assay_name)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- dplyr::as_tibble(pdat)
    pdat <- dplyr::mutate(pdat,
                          Row = rowData(se)[, "Row"],
                          Column = rowData(se)[, "Column"])

    pdat <- dplyr::left_join(pdat, dplyr::select(kmermap, "Row", "Column", "seq"),
                             by = c("Row", "Column"))

    ## group by k-mer sequence
    pdat_sets <- dplyr::select(pdat, -Row, -Column)
    pdat_sets <- tidyr::nest(pdat_sets, -seq)
    pdat_samples <- names(pdat_sets$data[[1]])
    pdat_sets <- dplyr::mutate(pdat_sets, data = lapply(data, as.matrix))

    ## calculate median intensities
    pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMedians, na.rm = TRUE))
    median_vals <- DataFrame(do.call(rbind, pdatm$m))
    names(median_vals) <- pdat_samples
    median_seqs <- pdatm$seq

    ## calculates mean intensities
    pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMeans2, na.rm = TRUE))
    mean_vals <- DataFrame(do.call(rbind, pdatm$m))
    names(mean_vals) <- pdat_samples
    mean_seqs <- pdatm$seq
    
    ## calculates mad intensities
    pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMads, na.rm = TRUE))
    mad_vals <- DataFrame(do.call(rbind, pdatm$m))
    names(mad_vals) <- pdat_samples
    mad_seqs <- pdatm$seq
    
    ## calculates SD intensities
    pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colSds, na.rm = TRUE))
    sd_vals <- DataFrame(do.call(rbind, pdatm$m))
    names(sd_vals) <- pdat_samples
    sd_seqs <- pdatm$seq

    ## calculates log2 mean intensities
    pdatm <- dplyr::mutate(pdat_sets,
                           m = lapply(data, function(x) {
                               matrixStats::colMeans2(log2(x + offset), na.rm = TRUE)
                           }))
    log2mean_vals <- DataFrame(do.call(rbind, pdatm$m))
    names(log2mean_vals) <- pdat_samples
    log2mean_seqs <- pdatm$seq

    ## calculates log2 mad intensities
    pdatm <- dplyr::mutate(pdat_sets,
                           m = lapply(data, function(x) {
                               matrixStats::colMads(log2(x + offset), na.rm = TRUE)
                               }))
    log2mad_vals <- DataFrame(do.call(rbind, pdatm$m))
    names(log2mad_vals) <- pdat_samples
    log2mad_seqs <- pdatm$seq
    
    ## calculates log2 SD intensities
    pdatm <- dplyr::mutate(pdat_sets,
                           m = lapply(data, function(x) {
                               matrixStats::colSds(log2(x + offset), na.rm = TRUE)
                           }))
    log2sd_vals <- DataFrame(do.call(rbind, pdatm$m))
    names(log2sd_vals) <- pdat_samples
    sd_seqs <- pdatm$seq
    
    ## calculates number of probes
    pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, nrow))
    n_vals <- unlist(pdatm$m)
    n_seqs <- pdatm$seq

    ## calculates number of NAs
    pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, function(x) { colSums(is.na(x)) }))
    na_vals <- DataFrame(do.call(rbind, pdatm$m))
    names(na_vals) <- pdat_samples
    na_seqs <- pdatm$seq

    ## check all in same order
    stopifnot(all(median_seqs == mean_seqs))
    stopifnot(all(median_seqs == mad_seqs))
    stopifnot(all(median_seqs == sd_seqs))
    stopifnot(all(median_seqs == n_seqs))
    stopifnot(all(median_seqs == na_seqs))
    stopifnot(all(median_seqs == log2mean_seqs))
    stopifnot(all(median_seqs == log2mad_seqs))
    stopifnot(all(median_seqs == log2sd_seqs))
    
    ## determine row data
    rowdat <- DataFrame(kmer = median_seqs,
                        nprobes = n_vals)
    
    ## create new SummarizedExperiment
    SummarizedExperiment(assays = list(medianIntensity = median_vals,
                                       meanIntensity = mean_vals,
                                       madIntensity = mad_vals,
                                       sdIntensity = sd_vals,
                                       naProbes = na_vals),
                                       log2MeanIntensity = log2mean_vals,
                                       log2MadIntensity = log2mad_vals,
                                       log2SdIntensity = log2sd_vals,
                         rowData = rowdat, colData = colData(se))
}
