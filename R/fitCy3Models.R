#' Compute Cy3 Scaling Factors
#'
#' PBM arrays are scanned twice, once for Cy3-tagged dUTPs to quantify
#' dsDNA abundance at each probe, and again for the Alexa488-tagged 
#' protein. The Cy3 scans can be used to first filter out probes which
#' appear to have poor dsDNA enrichment, and second, to scale Alexa488
#' intensities to account for differences in dsDNA abundance between
#' probes. Since probe intensities in the Cy3 scans are roughly
#' proportional to the number of adenines in the probe sequence, the
#' original Universal PBM Analysis Suite proposed estimating scaling factors
#' by looking at the residuals of a linear regression fit using all
#' tri-nucleotide sequences starting with an adenine. Given a
#' SummarizedExperiment of Cy3 scans, this function fit the tri-nucleotide
#' linear regression models and returns the expected Cy3 intensities,
#' the residuals from the fits, and the corresponding observed-expected
#' ratios for each probe as new assays added to the original Cy3
#' SummarizedExperiment object.
#'
#' The returned SummarizedExperiment object can be passed to \code{cy3Normalize}
#' to filter low quality probes and/or normalize Alexa488 intensities by the
#' computed ratios.
#' 
#' @param se SummarizedExperiment object containing PBM Cy3 intensity data.
#' @param assay_name string name of the assay to use. (default = "fore")
#' @param refit logical whether to filter outliers and refit trinucleotide
#'        linear regression model. (default = TRUE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to estimating affinities. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @param .trim interger vector of length two specifying start and end
#'        of probe sequence to be used. Default is based on the universal
#'        PBM probe design where only leading 36nt should be used. 
#'        Ignored if NULL. (default = c(1, 36))
#'
#' @return
#' original SummarizedExperiment object with additional assays corresponding
#' to the ratio of observed vs. expected probe intensities, and whether probes were
#' flagged as low-quality based on log2(ratio) > 1 or < -1. Cy3 models are stored
#' in the metadata of the returned object.
#'
#' @seealso cy3Normalize
#' @importFrom dplyr as_tibble select bind_cols group_by do ungroup left_join mutate
#' @importFrom tidyr gather spread
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' @export
#' @author Patrick Kimes
fitCy3Models <- function(se, assay_name = "fore", refit = TRUE, verbose = FALSE, .filter = 1L, 
                         .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {

    ## check Sequence info in rowData
    nse <- checkProbeSequences(se, FALSE)

    ## filter probes
    nse <- pbmFilterProbes(nse, .filter)

    ## trim probe sequences
    nse <- trimProbeSequences(nse, c(1, 36))

    ## determine row, sequence metadata
    rowdat <- as.data.frame(rowData(nse), optional = TRUE)
    rowdat <- dplyr::as_tibble(rowdat)
    rowdat <- dplyr::select(rowdat, Row, Column, Sequence)

    nt_freq <- Biostrings::DNAStringSet(rowdat$Sequence)
    nt_freq <- Biostrings::oligonucleotideFrequency(nt_freq, 3)
    nt_freq <- dplyr::as_tibble(nt_freq)
    nt_freq <- dplyr::select(nt_freq, starts_with("A"))

    rowdat <- dplyr::bind_cols(rowdat, nt_freq)

    ## extract intensities
    pdat <- assay(nse, assay_name)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- dplyr::as_tibble(pdat)
    pdat <- dplyr::bind_cols(pdat, rowdat)
    pdat <- tidyr::gather(pdat, condition, intensity, colnames(se))

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
                           coords = lapply(data, `[`, c("Row", "Column")))
    pexps <- dplyr::select(pexps, condition, preds, coords)
    pexps <- tidyr::unnest(pexps)
    pexps <- tidyr::spread(pexps, condition, preds)

    ## compute ratios
    pratios <- tidyr::gather(pexps, condition, preds, -Row, -Column)
    pratios <- dplyr::left_join(pratios, dplyr::select(pdat, Row, Column, condition, intensity),
                                by = c("Row", "Column", "condition"))
    pratios <- dplyr::mutate(pratios, ratio = intensity / preds)
    pratios <- dplyr::select(pratios, Row, Column, condition, ratio)
    pratios <- tidyr::spread(pratios, condition, ratio)

    ## compute low quality probes (outside 2-fold difference)
    pdrop <- tidyr::gather(pratios, condition, ratio, -Row, -Column)
    pdrop <- dplyr::mutate(pdrop, lowq = ratio > 2 | ratio < 1/2)    
    pdrop <- dplyr::select(pdrop, Row, Column, condition, lowq)
    
    ## refit data excluding outliers
    if (refit) {
        ## add initial fits to metadata for reference
        initial_fits <- pfits$fit
        names(initial_fits) <- pfits$condition
        metadata(se)$cy3models_init <- initial_fits

        ## refit models after setting outliers to NA
        pdat <- dplyr::left_join(pdat, pdrop, by = c("Row", "Column", "condition"))
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
                               coords = lapply(data, `[`, c("Row", "Column")))
        pexps <- dplyr::select(pexps, condition, preds, coords)
        pexps <- tidyr::unnest(pexps)
        pexps <- tidyr::spread(pexps, condition, preds)

        ## re-compute ratios
        pratios <- tidyr::gather(pexps, condition, preds, -Row, -Column)
        pratios <- dplyr::left_join(pratios, dplyr::select(pdat, Row, Column, condition, intensity),
                                    by = c("Row", "Column", "condition"))
        pratios <- dplyr::mutate(pratios, ratio = intensity / preds)
        pratios <- dplyr::select(pratios, Row, Column, condition, ratio)
        pratios <- tidyr::spread(pratios, condition, ratio)
    }
    pdrop <- tidyr::spread(pdrop, condition, lowq)
    
    ## left join to original rowData to get full set
    full_rowdat <- as.data.frame(rowData(se), optional = TRUE)
    full_rowdat <- dplyr::as_tibble(full_rowdat)
    full_rowdat <- dplyr::select(full_rowdat, Row, Column, Sequence)
    pexps <- dplyr::left_join(dplyr::select(full_rowdat, Row, Column), pexps, by = c("Row", "Column"))
    pratios <- dplyr::left_join(dplyr::select(full_rowdat, Row, Column), pratios, by = c("Row", "Column"))
    pdrop <- dplyr::left_join(dplyr::select(full_rowdat, Row, Column), pdrop, by = c("Row", "Column"))

    pexps <- DataFrame(pexps, check.names = FALSE)
    pexps <- pexps[, rownames(colData(se))]
    pratios <- DataFrame(pratios, check.names = FALSE)
    pratios <- pratios[, rownames(colData(se))]
    pdrop <- DataFrame(pdrop, check.names = FALSE)
    pdrop <- pdrop[, rownames(colData(se))]

    ## add to SE object
    assay(se, "expected") <- pexps
    assay(se, "ratio") <- pratios
    assay(se, "lowq") <- pdrop

    ## add fits to metadata
    fits <- pfits$fit
    names(fits) <- pfits$condition
    metadata(se)$cy3models <- fits

    return(se)
}
