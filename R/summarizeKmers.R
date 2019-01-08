#' Summarize K-mer Intensities
#'
#' This function calculates k-mer intensity information from
#' PBM data. 
#' 
#' @param se SummarizedExperiment object containing PBM intensity data.
#' @param assay string name of the assay to use. (default = \code{SummarizedExperiment::assayNames(se)[1]})
#' @param kmers character vector of k-mers to predict.
#' @param stat_set character vector of statistics to calculate for each sample.
#'        The set of supported statistics are listed in the details. By default,
#'        all possible statistics are computed. Details on available statistics are
#'        given in details. (default = \code{c("median", "mean", "mad", "sd", "log2mean", "log2mad", "log2sd", "na", "quantile")})
#' @param weights tibble of sequence position weights (log2) to use when computing
#'        metrics. Can be calculated using \code{approxPositionWeights}. Ignored
#'        if NULL. (default = NULL)
#' @param offset integer offset to add to intensities before log2 scaling to
#'        prevent errors with zero intensities. If set to 0, probes with
#'        zero intensities are dropped/ignored for log-scaled metrics. (default = 1)
#' @param q probability values or vector of values which should be used for computing
#'        quantiles if "quantile" is specified as part of \code{stat_set}.
#'        (default = 0.25)
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
#' * \code{"log2mean"}: mean probe log2(intensity + offset)
#' * \code{"log2mad"}: MAD probe log2(intensity + offset)
#' * \code{"log2sd"}: SD probe log(intensity + offset)
#' * \code{"na"}: number of NA probes
#' * \code{"quantile"}: q-quantile probe intensity
#' 
#' @return
#' SummarizedExperiment object with intensity information summarized for
#' k-mers given in separate assays, with each row corresponding to a separate
#' k-mer. k-mer sequences and number of corresponding probes are given
#' in the rowdata of the SummarizedBenchmark.
#' 
#' @md
#' @importFrom stats mad median
#' @importFrom dplyr as_tibble mutate left_join select
#' @importFrom tidyr nest
#' @importFrom matrixStats colMedians colMeans2 colSds colMads colQuantiles
#' @export
#' @author Patrick Kimes
summarizeKmers <- function(se, assay = SummarizedExperiment::assayNames(se)[1], kmers = NULL, 
                           stat_set = c("median", "mean", "mad", "sd", "log2mean",
                                        "log2mad", "log2sd", "na", "quantile"),
                           offset = 1, weights = NULL, q = 0.25, verbose = FALSE, .filter = 1L,
                           .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {

    ## check stat statistics are valid
    stat_set <- match.arg(stat_set, several.ok = TRUE)
    if (length(stat_set) == 0) {
        stop("Please specify at least one valid statistic.")
    }

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
    kmermap <- mapKmers(rowData(se)[, ovnames, drop = FALSE], kmers)

    ## use ordering from input 'kmers'
    kmermap$seq <- factor(kmermap$seq, levels = kmers)
    
    ## extract intensities
    pdat <- SummarizedExperiment::assay(se, assay)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- dplyr::as_tibble(pdat)

    ## check whether row/column indices were available

    if (all(c("Row", "Column") %in% names(kmermap))) {
        pdat <- dplyr::mutate(pdat,
                              Row = rowData(se)[, "Row"],
                              Column = rowData(se)[, "Column"])
        
        if (is.null(weights)) {
            pdat <- dplyr::left_join(pdat, dplyr::select(kmermap, "Row", "Column", "seq"),
                                     by = c("Row", "Column"))
        } else {
            pdat <- dplyr::left_join(pdat, dplyr::select(kmermap, "Row", "Column", "seq", "pos"),
                                     by = c("Row", "Column"))
            pdat <- dplyr::left_join(tidyr::gather(pdat, sample, value, -seq, -pos, -Row, -Column),
                                     tidyr::gather(weights, sample, l2r, -pos),
                                     by = c("sample", "pos"))
            pdat <- dplyr::mutate(pdat, value = value / 2^(l2r))
            pdat <- dplyr::select(pdat, -l2r)
            pdat <- tidyr::spread(pdat, sample, value)
            pdat <- dplyr::select(pdat, -pos)
        }
        pdat <- dplyr::select(pdat, -Row, -Column)
    } else {
        pdat <- dplyr::mutate(pdat, probe_idx = 1:n())
        if (is.null(weights)) {
            pdat <- dplyr::left_join(pdat, dplyr::select(kmermap, "probe_idx", "seq"),
                                     by = "probe_idx")
        } else {
            pdat <- dplyr::left_join(pdat, dplyr::select(kmermap, "probe_idx", "seq", "pos"),
                                     by = "probe_idx")
            pdat <- dplyr::left_join(tidyr::gather(pdat, sample, value, -seq, -pos, -probe_idx),
                                     tidyr::gather(weights, sample, l2r, -pos),
                                     by = c("sample", "pos"))
            pdat <- dplyr::mutate(pdat, value = value / 2^(l2r))
            pdat <- dplyr::select(pdat, -l2r)
            pdat <- tidyr::spread(pdat, sample, value)
            pdat <- dplyr::select(pdat, -pos)
        }
        pdat <- dplyr::select(pdat, -probe_idx)
    }
    
    ## group by k-mer sequence
    pdat_sets <- tidyr::nest(pdat, -seq)
    pdat_samples <- names(pdat_sets$data[[1]])
    pdat_seqs <- pdat_sets$seq
    pdat_sets <- dplyr::mutate(pdat_sets, data = lapply(data, as.matrix))

    ## calculate number of probes
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
        median_vals <- median_vals[colnames(se)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$medianIntensity <- median_vals
    }
    
    ## calculate mean intensities
    if ("mean" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMeans2, na.rm = TRUE))
        mean_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(mean_vals) <- pdat_samples
        mean_vals <- mean_vals[colnames(se)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$meanIntensity <- mean_vals
    }
    
    ## calculate mad intensities
    if ("mad" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMads, na.rm = TRUE))
        mad_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(mad_vals) <- pdat_samples
        mad_vals <- mad_vals[colnames(se)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$madIntensity <- mad_vals
    }

    ## calculate SD intensities
    if ("sd" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colSds, na.rm = TRUE))
        sd_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(sd_vals) <- pdat_samples
        sd_vals <- sd_vals[colnames(se)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$sdIntensity <- sd_vals
    }
    
    ## calculate log2 mean intensities
    if ("log2mean" %in% stat_set) {
        if (offset <= 0) {  
            pdatm <- dplyr::mutate(pdat_sets,
                                   m = lapply(data, function(x) {
                                       x[x <= 0] <- NA
                                       matrixStats::colMeans2(log2(x), na.rm = TRUE)
                                   }))
        } else {
            pdatm <- dplyr::mutate(pdat_sets,
                                   m = lapply(data, function(x) {
                                       matrixStats::colMeans2(log2(x + offset), na.rm = TRUE)
                                   }))
        }
        log2mean_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(log2mean_vals) <- pdat_samples
        log2mean_vals <- log2mean_vals[colnames(se)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$log2meanIntensity <- log2mean_vals
    }

    ## calculate log2 mad intensities
    if ("log2mad" %in% stat_set) {
        if (offset <= 0) {
            pdatm <- dplyr::mutate(pdat_sets,
                                   m = lapply(data, function(x) {
                                       x[x <= 0] <- NA
                                       matrixStats::colMads(log2(x), na.rm = TRUE)
                                   }))
        } else {
            pdatm <- dplyr::mutate(pdat_sets,
                                   m = lapply(data, function(x) {
                                       matrixStats::colMads(log2(x + offset), na.rm = TRUE)
                                   }))
        }
        log2mad_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(log2mad_vals) <- pdat_samples
        log2mad_vals <- log2mad_vals[colnames(se)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$log2madIntensity <- log2mad_vals
    }
    
    ## calculate log2 SD intensities
    if ("log2sd" %in% stat_set) {
        if (offset <= 0) {
            pdatm <- dplyr::mutate(pdat_sets,
                                   m = lapply(data, function(x) {
                                       matrixStats::colSds(log2(x[x > 0]), na.rm = TRUE)
                                   }))
        } else {
            pdatm <- dplyr::mutate(pdat_sets,
                                   m = lapply(data, function(x) {
                                       matrixStats::colSds(log2(x + offset), na.rm = TRUE)
                                   }))
        }
        log2sd_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(log2sd_vals) <- pdat_samples
        log2sd_vals <- log2sd_vals[colnames(se)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$log2sdIntensity <- log2sd_vals
    }
        
    ## calculate number of NAs
    if ("na" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, function(x) { colSums(is.na(x)) }))
        na_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(na_vals) <- pdat_samples
        na_vals <- na_vals[colnames(se)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$naProbes <- na_vals
    }

    ## calculate quantile intensities
    if ("quantile" %in% stat_set) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colQuantiles,
                                                     probs = q, drop = FALSE, na.rm = TRUE))
        qtiles <- colnames(pdatm$m[[1]])
        qtiles <- paste0("q", gsub("%", "Intensity", qtiles))
        for (i in seq_len(length(qtiles))) {
            qt_vals <- lapply(pdatm$m, function(x) { x[, i] })
            qt_vals <- DataFrame(do.call(rbind, qt_vals))
            names(qt_vals) <- pdat_samples
            qt_vals <- qt_vals[colnames(se)]
            stopifnot(pdat_seqs == pdatm$seq)
            assay_list[[qtiles[[i]]]] <- qt_vals
        }
    }
        
    ## determine row data
    rowdat <- DataFrame(kmer = as.character(pdat_seqs),
                        nprobes = n_vals)
    
    ## create new SummarizedExperiment
    SummarizedExperiment(assays = assay_list,
                         rowData = rowdat, colData = colData(se))
}

