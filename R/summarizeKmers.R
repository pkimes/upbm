#' @title Compute K-mer summary metrics from probe intensities
#'
#' @description
#' This function calculates k-mer intensity information from
#' PBM data. 
#' 
#' @param pe a PBMExperiment object containing PBM intensity data.
#' @param assay a string name of the assay to adjust. (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param kmers a character vector of k-mers to predict. (default = \code{uniqueKmers(8L)})
#' @param metrics a character vector of statistics to calculate for each scan.
#'        The set of supported statistics are listed in the details. By default,
#'        all possible statistics are computed, but this is not recommended as
#'        it can be computational expensive. (default = \code{c("median", "mean", "mad", "sd", "log2mean", "log2mad", "log2sd", "na", "quantile")})
#' @param offset an integer offset to add to intensities before log2 scaling to
#'        prevent errors with zero intensities. If set to 0, probes with
#'        zero intensities are dropped/ignored for log-scaled metrics. (default = 1)
#' @param q a vector of quantiles which should be used if \code{"quantile"} is specified as part
#'        of \code{metrics}. (default = 0.25)
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
#' k-mers. Each metric is returned as a separate assay, with each rows corresponding
#' to k-mers. k-mer sequences and the number of corresponding probes are given
#' in the rowdata of the object.
#' 
#' @importFrom stats mad median
#' @importFrom dplyr as_tibble mutate left_join select
#' @importFrom tidyr nest
#' @importFrom matrixStats colMedians colMeans2 colSds colMads colQuantiles
#' @export
#' @author Patrick Kimes
summarizeKmers <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                           kmers = uniqueKmers(8L), 
                           metrics = c("median", "mean", "mad", "sd", "log2mean",
                                       "log2mad", "log2sd", "na", "quantile"),
                           offset = 1, q = 0.25) {
    stopifnot(is(pe, "PBMExperiment")) 

    ## check stat statistics are valid
    metrics <- match.arg(metrics, several.ok = TRUE)
    if (length(metrics) == 0) {
        stop("Please specify at least one valid statistic.")
    }

    ## check kmers specified
    if (!is.vector(kmers, mode = "character")) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences as character strings.")
    }
    if (length(unique(nchar(kmers))) != 1L) {
        stop("If specified, 'kmers' must be a vector of nucleotide sequences of equal length.")
    }
    
    ## filter probes
    pe <- pbmFilterProbes(pe)

    ## trim probe sequences
    pe <- pbmTrimProbes(pe)

    ## filtered/trimmed PBMDesign
    pdes <- PBMDesign(pe)

    pdest <- as.data.frame(pdes@design, optional = TRUE)
    pdest <- dplyr::as_tibble(pdest)

    ## find mapping between kmers and probes
    kmermap <- mapKmers(pdes, kmers)

    ## extract intensities
    pdat <- SummarizedExperiment::assay(pe, assay)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- dplyr::as_tibble(pdat)

    ## check whether row/column indices were available
    pdat <- cbind(pdat, pdest)
    pdat <- dplyr::inner_join(pdat, dplyr::select_at(kmermap, c(pe@probeCols, "seq")),
                              by = pe@probeCols)
    pdat <- dplyr::select(pdat, -(!! pe@probeCols))
    
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
    if ("median" %in% metrics) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMedians, na.rm = TRUE))
        median_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(median_vals) <- pdat_samples
        median_vals <- median_vals[colnames(pe)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$medianIntensity <- median_vals
    }
    
    ## calculate mean intensities
    if ("mean" %in% metrics) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMeans2, na.rm = TRUE))
        mean_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(mean_vals) <- pdat_samples
        mean_vals <- mean_vals[colnames(pe)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$meanIntensity <- mean_vals
    }
    
    ## calculate mad intensities
    if ("mad" %in% metrics) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colMads, na.rm = TRUE))
        mad_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(mad_vals) <- pdat_samples
        mad_vals <- mad_vals[colnames(pe)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$madIntensity <- mad_vals
    }

    ## calculate SD intensities
    if ("sd" %in% metrics) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colSds, na.rm = TRUE))
        sd_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(sd_vals) <- pdat_samples
        sd_vals <- sd_vals[colnames(pe)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$sdIntensity <- sd_vals
    }
    
    ## calculate log2 mean intensities
    if ("log2mean" %in% metrics) {
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
        log2mean_vals <- log2mean_vals[colnames(pe)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$log2meanIntensity <- log2mean_vals
    }

    ## calculate log2 mad intensities
    if ("log2mad" %in% metrics) {
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
        log2mad_vals <- log2mad_vals[colnames(pe)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$log2madIntensity <- log2mad_vals
    }
    
    ## calculate log2 SD intensities
    if ("log2sd" %in% metrics) {
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
        log2sd_vals <- log2sd_vals[colnames(pe)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$log2sdIntensity <- log2sd_vals
    }
        
    ## calculate number of NAs
    if ("na" %in% metrics) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, function(x) { colSums(is.na(x)) }))
        na_vals <- DataFrame(do.call(rbind, pdatm$m))
        names(na_vals) <- pdat_samples
        na_vals <- na_vals[colnames(pe)]
        stopifnot(pdat_seqs == pdatm$seq)
        assay_list$naProbes <- na_vals
    }

    ## calculate quantile intensities
    if ("quantile" %in% metrics) {
        pdatm <- dplyr::mutate(pdat_sets, m = lapply(data, matrixStats::colQuantiles,
                                                     probs = q, drop = FALSE, na.rm = TRUE))
        qtiles <- colnames(pdatm$m[[1]])
        qtiles <- paste0("q", gsub("%", "Intensity", qtiles))
        for (i in seq_len(length(qtiles))) {
            qt_vals <- lapply(pdatm$m, function(x) { x[, i] })
            qt_vals <- DataFrame(do.call(rbind, qt_vals))
            names(qt_vals) <- pdat_samples
            qt_vals <- qt_vals[colnames(pe)]
            stopifnot(pdat_seqs == pdatm$seq)
            assay_list[[qtiles[[i]]]] <- qt_vals
        }
    }
        
    ## determine row data
    rowdat <- DataFrame(kmer = as.character(pdat_seqs),
                        nprobes = n_vals)
    
    ## create new SummarizedExperiment
    SummarizedExperiment(assays = assay_list,
                         rowData = rowdat,
                         colData = colData(pe))
}

