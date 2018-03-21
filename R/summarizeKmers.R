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
#' @param stat_set caracter vector of statistics to calculate for each sample.
#'        The set of supported statistics are listed in the details. By default,
#'        only a subset of the total possible statistics are computed.
#'        If all supported statistics should be calculated, specify NULL.
#'        (default = \code{c("median", "log2mad", "lo2sd", "na")})
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
#' @details
#' The following statistics are currently supported.
#' * \code{"median"}: median probe intensity
#' * \code{"mean"}: mean probe intensity
#' * \code{"mad"}: MAD probe intensity (scaled by 1.4826)
#' * \code{"sd"}: SD probe intensity
#' * \code{"log2mean"}: mean probe log2(intensity + shift)
#' * \code{"log2mad"}: MAD probe log2(intensity + shift)
#' * \code{"log2sd"}: SD probe log(intensity + shift)
#' * \code{"na"}: number of NA probes
#' 
#' @return
#' SummarizedExperiment object with intensity information summarized for
#' k-mers given in separate assays, with each row corresponding to a separate
#' k-mer. k-mer sequences and number of corresponding probes are given
#' in the rowdata of the SummarizedBenchmark.
#' 
#' @md
#' @importFrom dplyr as_tibble mutate left_join select
#' @importFrom tidyr nest
#' @importFrom matrixStats colMedians colMeans2 colSds colMads
#' @export
#' @author Patrick Kimes
summarizeKmers <- function(se, assay_name = "fore", kmers = NULL, offset = 1,
                           stat_set = c("median", "log2mad", "lo2sd", "na"),
                           verbose = FALSE, .filter = 1L,
                           .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {

    ## check stat statistics are valid
    valid_stats <- c("median", "mean", "mad", "sd",
                     "log2mean", "log2mad", "log2sd", "na")
    if (is.null(stat_set)) {
        stat_set <- valid_stats
    } else {
        stat_set <- valid_stats[pmatch(stat_set, valid_stats)]
        stat_set <- stat_set[!is.na(stat_set)]
    }
    stopifnot(length(stat_set) > 0)
    
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
    pdat_seqs <- pdat_sets$seq
    pdat_sets <- dplyr::mutate(pdat_sets, data = lapply(data, as.matrix))

    ## calculates number of probes
    pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, nrow))
    n_vals <- unlist(pdatm$m)
    n_seqs <- pdatm$seq

    ## store statistics in separate assays
    assay_list <- list()
    
    ## calculate median intensities
    if ("median" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMedians, na.rm = TRUE))
        median_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(median_vals) <- pdat_samples
        stopifnot(all(pdat_seqs == pdatm$seq))
        assay_list$medianIntensity <- median_vals
    }
    
    ## calculates mean intensities
    if ("mean" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMeans2, na.rm = TRUE))
        mean_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(mean_vals) <- pdat_samples
        stopifnot(all(pdat_seqs == pdatm$seq))
        assay_list$meanIntensity <- mean_vals
    }
    
    ## calculates mad intensities
    if ("mad" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMads, na.rm = TRUE))
        mad_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(mad_vals) <- pdat_samples
        stopifnot(all(pdat_seqs == pdatm$seq))
        assay_list$madIntensity <- mad_vals
    }

    ## calculates SD intensities
    if ("sd" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colSds, na.rm = TRUE))
        sd_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(sd_vals) <- pdat_samples
        stopifnot(all(pdat_seqs == pdatm$seq))
        assay_list$sdIntensity <- sd_vals
    }
    
    ## calculates log2 mean intensities
    if ("log2mean" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets,
                               m = lapply(data, function(x) {
                                   matrixStats::colMeans2(log2(x + offset), na.rm = TRUE)
                               }))
        log2mean_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(log2mean_vals) <- pdat_samples
        stopifnot(all(pdat_seqs == pdatm$seq))
        assay_list$log2meanIntensity <- log2mean_vals
    }

    ## calculates log2 mad intensities
    if ("log2mad" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets,
                               m = lapply(data, function(x) {
                                   matrixStats::colMads(log2(x + offset), na.rm = TRUE)
                               }))
        log2mad_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(log2mad_vals) <- pdat_samples
        stopifnot(all(pdat_seqs == pdatm$seq))
        assay_list$log2madIntensity <- log2mad_vals
    }
    
    ## calculates log2 SD intensities
    if ("log2sd" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets,
                               m = lapply(data, function(x) {
                                   matrixStats::colSds(log2(x + offset), na.rm = TRUE)
                               }))
        log2sd_vals <- DataFrame(do.call(rbind, pdatm$m))
        stopifnot(all(pdat_seqs == pdatm$seq))
        sd_seqs <- pdatm$seq
        assay_list$log2sdIntensity <- log2sd_vals
    }
        
    ## calculates number of NAs
    if ("na" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, function(x) { colSums(is.na(x)) }))
        na_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(na_vals) <- pdat_samples
        stopifnot(all(pdat_seqs == pdatm$seq))
        assay_list$naProbes <- na_vals
    }
    
    ## determine row data
    rowdat <- DataFrame(kmer = median_seqs,
                        nprobes = n_vals)
    
    ## create new SummarizedExperiment
    SummarizedExperiment(assays = assay_list,
                         rowData = rowdat, colData = colData(se))
}