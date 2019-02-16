#' @title Compute Cy3 scaling factors using regression
#'
#' @description
#' PBM arrays are scanned twice, once for Cy3-tagged dUTPs to quantify
#' dsDNA abundance at each probe, and again for the Alexa488-tagged 
#' protein. The Cy3 scans can be used to first filter out probes which
#' appear to have poor dsDNA enrichment, and second, to scale Alexa488
#' intensities to account for differences in dsDNA abundance between
#' probes.
#'
#' \emph{Both this function and the \code{cy3FitEmpirical} function can
#' be used to estimate scaling factors for filtering and scaling
#' Alexa488 probe intensities using Cy3 information. The \code{cy3FitEmpirical}
#' approach is recommended when an appropriate reference set of Cy3 scans
#' are available.}
#'
#' Since probe intensities in the Cy3 scans are roughly
#' proportional to the number of adenines in the probe sequence, the
#' original Universal PBM Analysis Suite proposed modeling the observed Cy3
#' intensities with linear regression using the counts of all tri-nucleotide sequences
#' starting with an adenine as covariates. The observed-to-expected intensity ratios
#' are then used as the scaling factors.
#'
#' Given a PBMExperiment of Cy3 scans, this function fits the tri-nucleotide
#' regression models and returns the expected Cy3 intensities,
#' the residuals from the fits, and the corresponding observed-to-expected
#' ratios for each probe as new assays added to the original Cy3
#' PBMExperiment object.
#'
#' The returned PBMExperiment object can be passed to \code{cy3Normalize}
#' to filter low quality probes and/or normalize Alexa488 intensities by the
#' computed ratios.
#' 
#' @param pe a PBMExperiment object containing Cy3 intensity data.
#' @param assay a numeric index or string name specifying the assay to use.
#'        (default = \code{SummarizedExperiment::assayNames(pe)[1]})
#' @param refit a logical value whether to filter outliers and refit trinucleotide
#'        linear regression model. (default = TRUE)
#' @param threshold a numeric threshold on absolute value of log2 ratio between
#'        observed and expected Cy3 intensities. (default = 1/2)
#' @param verbose a logical value whether to print verbose output during
#'        analysis. (default = FALSE)
#'
#' @return
#' Original Cy3 PBMExperiment object with additional assays corresponding
#' to ratio of observed to expected probe intensities, and whether probes were
#' flagged as low-quality based on \code{abs(log2(ratio)) > threshold}.
#' Cy3 models are stored in the metadata of the returned object.
#'
#' @references
#' \itemize{
#' \item Berger, M. F., & Bulyk, M. L. (2009). Universal protein-binding microarrays for the comprehensive characterization of the DNA-binding specificities of transcription factors. Nature Protocols, 4(3), 393-411.
#' }
#'
#' @seealso cy3Normalize cy3FitEmpirical
#' @importFrom stats lm na.exclude predict
#' @importFrom dplyr as_tibble select select_at bind_cols group_by do ungroup left_join mutate starts_with
#' @importFrom tidyr gather spread
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' @export
#' @author Patrick Kimes
cy3FitModel <- function(pe, assay = SummarizedExperiment::assayNames(pe)[1],
                        refit = TRUE, threshold = 1/2, verbose = FALSE) {
    stopifnot(is(pe, "PBMExperiment")) 

    if (verbose) {
        cat("|| upbm::cy3FitModel \n")
        cat("|| - Starting calculation of Cy3 deviations from OLS fit",
            "for", ncol(pe), "Cy3 PBM scans.\n")
    }

    if (verbose) {
        cat("|| - Filtering probes according to", length(pe@probeFilter),
            "probeFilter rule(s).\n")
        ntotal <- nrow(pe)
    }

    ## filter probes
    npe <- pbmFilterProbes(pe)

    if (verbose) {
        cat("|| - Data filtered from", ntotal, "probes to", nrow(npe), "probes.\n")
    }

    if (verbose && length(pe@probeTrim > 0L)) {
        cat("|| - Trimming probes according to probeTrim settings:",
            paste0("[", paste0(npe@probeTrim, collapse = ", "), "]"), "\n")
    }

    ## trim probe sequences
    npe <- pbmTrimProbes(npe)

    ## determine row, sequence metadata
    rowdat <- as.data.frame(rowData(npe), optional = TRUE)
    rowdat <- dplyr::as_tibble(rowdat)
    rowdat <- dplyr::select_at(rowdat, npe@probeCols)

    nt_freq <- Biostrings::DNAStringSet(rowdat$Sequence)
    nt_freq <- Biostrings::oligonucleotideFrequency(nt_freq, 3)
    nt_freq <- dplyr::as_tibble(nt_freq)
    nt_freq <- dplyr::select(nt_freq, dplyr::starts_with("A"))

    rowdat <- dplyr::bind_cols(rowdat, nt_freq)

    ## extract intensities
    pdat <- SummarizedExperiment::assay(npe, assay)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- dplyr::as_tibble(pdat)
    pdat <- dplyr::bind_cols(pdat, rowdat)
    pdat <- tidyr::gather(pdat, condition, intensity, colnames(pe))

    if (verbose) {
        cat("|| - Fitting trinucleotide models to probes.\n")
    }

    ## fit models
    pdat <- dplyr::group_by(pdat, condition)
    pfits <- dplyr::do(pdat,
                       fit = lm(intensity ~ AAA + AAC + AAG + AAT + ACA +
                                    ACC + ACG + ACT + AGA + AGC + AGG +
                                    AGT + ATA + ATC + ATG + ATT, data = .,
                                na.action = na.exclude))
    pfits <- dplyr::ungroup(pfits)
    pdat <- dplyr::ungroup(pdat)
    
    ## compute expected values
    pexps <- dplyr::left_join(pfits, tidyr::nest(pdat, -condition), by = "condition")
    pexps <- dplyr::mutate(pexps, preds = mapply(predict, object = fit, newdata = data,
                                                 SIMPLIFY = FALSE),
                           probecols = lapply(data, `[`, npe@probeCols))
    pexps <- dplyr::select(pexps, condition, preds, probecols)
    pexps <- tidyr::unnest(pexps)
    pexps <- tidyr::spread(pexps, condition, preds)

    ## compute ratios
    pratios <- tidyr::gather(pexps, condition, preds, -(!! npe@probeCols))
    pratios <- dplyr::left_join(pratios, dplyr::select_at(pdat, c("condition", "intensity", npe@probeCols)),
                                by = c("condition", npe@probeCols))
    pratios <- dplyr::mutate(pratios, ratio = intensity / preds)
    pratios <- dplyr::select_at(pratios, c("condition", "ratio", npe@probeCols))
    pratios <- tidyr::spread(pratios, condition, ratio)

    ## compute low quality probes (outside 2-fold difference)
    pdrop <- tidyr::gather(pratios, condition, ratio, -(!! npe@probeCols))
    pdrop <- dplyr::mutate(pdrop, lowq = ratio > 2^threshold | ratio < 2^(-threshold))
    pdrop <- dplyr::select_at(pdrop, c("condition", "lowq", npe@probeCols))
    
    ## refit data excluding outliers
    if (refit) {
        ## add initial fits to metadata for reference
        initial_fits <- pfits$fit
        names(initial_fits) <- pfits$condition
        metadata(pe)$cy3models_init <- initial_fits

        if (verbose) {
            cat("|| - Refitting trinucleotide models to probes passing filtering rules.\n")
        }
        
        ## refit models after setting outliers to NA
        pdat <- dplyr::left_join(pdat, pdrop, by = c("condition", npe@probeCols))
        pdat <- dplyr::mutate(pdat, intensity = ifelse(lowq, NA, intensity))
        pdat <- dplyr::group_by(pdat, condition)
        pfits <- dplyr::do(pdat,
                           fit = lm(intensity ~ AAA + AAC + AAG + AAT + ACA +
                                        ACC + ACG + ACT + AGA + AGC + AGG +
                                        AGT + ATA + ATC + ATG + ATT, data = .,
                                    na.action = na.exclude))
        pfits <- dplyr::ungroup(pfits)
        pdat <- dplyr::ungroup(pdat)
        
        ## re-compute expected values
        pexps <- dplyr::left_join(pfits, tidyr::nest(pdat, -condition), by = "condition")
        pexps <- dplyr::mutate(pexps, preds = mapply(predict, object = fit, newdata = data,
                                                     SIMPLIFY = FALSE),
                               probecols = lapply(data, `[`, npe@probeCols))
        pexps <- dplyr::select(pexps, condition, preds, probecols)
        pexps <- tidyr::unnest(pexps)
        pexps <- tidyr::spread(pexps, condition, preds)

        ## re-compute ratios
        pratios <- tidyr::gather(pexps, condition, preds, -(!! npe@probeCols))
        pratios <- dplyr::left_join(pratios, dplyr::select_at(pdat, c("condition", "intensity", npe@probeCols)),
                                    by = c("condition", npe@probeCols))
        pratios <- dplyr::mutate(pratios, ratio = intensity / preds)
        pratios <- dplyr::select_at(pratios, c("condition", "ratio", npe@probeCols))
        pratios <- tidyr::spread(pratios, condition, ratio)
    }
    pdrop <- tidyr::spread(pdrop, condition, lowq)
    
    ## left join to original rowData to get full set
    full_rowdat <- as.data.frame(rowData(pe), optional = TRUE)
    full_rowdat <- dplyr::as_tibble(full_rowdat)
    full_rowdat <- dplyr::select_at(full_rowdat, npe@probeCols)
    pexps <- dplyr::left_join(dplyr::select_at(full_rowdat, npe@probeCols), pexps, by = npe@probeCols)
    pratios <- dplyr::left_join(dplyr::select_at(full_rowdat, npe@probeCols), pratios, by = npe@probeCols)
    pdrop <- dplyr::left_join(dplyr::select_at(full_rowdat, npe@probeCols), pdrop, by = npe@probeCols)
    
    pexps <- DataFrame(pexps, check.names = FALSE)
    pexps <- pexps[, rownames(colData(pe)), drop = FALSE]
    pratios <- DataFrame(pratios, check.names = FALSE)
    pratios <- pratios[, rownames(colData(pe)), drop = FALSE]
    pdrop <- DataFrame(pdrop, check.names = FALSE)
    pdrop <- pdrop[, rownames(colData(pe)), drop = FALSE]

    ## add to SE object
    SummarizedExperiment::assay(pe, "expected") <- pexps
    SummarizedExperiment::assay(pe, "ratio") <- pratios
    SummarizedExperiment::assay(pe, "lowq") <- pdrop

    ## add fits to metadata
    fits <- pfits$fit
    names(fits) <- pfits$condition
    metadata(pe)$cy3models <- fits

    if (verbose) {
        cat("|| - Finished calculation of Cy3 deviations.\n")
        cat("|| - Returning PBMExperiment with", nrow(pe), "rows and", ncol(pe), "columns.\n")
    }
    return(pe)
}
