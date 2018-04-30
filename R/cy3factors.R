#' Compute Cy3 Scaling Factors
#'
#' PBM arrays are scanned twice, once for Cy3-tagged dUTPs to quantify
#' dsDNA abundance at each probe, and again for the Alexa488-tagged 
#' protein.
#' 
#' @param se SummarizedExpierment object containing PBM Cy3 intensity data.
#' @param assay_name string name of the assay to use. (default = "fore")
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
#' to residuals, expected intensities, and ratio of observed vs. expected
#' probe intensities. Cy3 models are stored in the metadata of the returned
#' object.
#'
#' @importFrom dplyr as_tibble select bind_cols group_by do ungroup left_join mutate
#' @importFrom tidyr gather spread
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' @export
#' @author Patrick Kimes
cy3factors <- function(se, assay_name = "fore", .filter = 1L, verbose = FALSE,
                       .trim = if (.filter > 0L) { c(1, 36) } else { NULL }) {

    ## check Sequence info in rowData
    nse <- checkProbeSequences(se, FALSE)

    ## filter probes
    nse <- pbmFilterProbes(nse, 1)

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
    pdat <- assay(nse, "fore")
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

    ## compute residuals
    presids <- dplyr::mutate(pfits, fit = lapply(fit, resid))
    presids <- tidyr::unnest(presids)
    presids <- dplyr::bind_cols(presids, dplyr::select(ungroup(pdat), Row, Column))
    presids <- tidyr::spread(presids, condition, fit)

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

    ## left join to original rowData to get full set
    full_rowdat <- as.data.frame(rowData(se), optional = TRUE)
    full_rowdat <- dplyr::as_tibble(full_rowdat)
    full_rowdat <- dplyr::select(full_rowdat, Row, Column, Sequence)
    presids <- dplyr::left_join(dplyr::select(full_rowdat, Row, Column), presids, by = c("Row", "Column"))
    pexps <- dplyr::left_join(dplyr::select(full_rowdat, Row, Column), pexps, by = c("Row", "Column"))
    pratios <- dplyr::left_join(dplyr::select(full_rowdat, Row, Column), pratios, by = c("Row", "Column"))

    presids <- DataFrame(presids, check.names = FALSE)
    presids <- presids[, rownames(colData(se))]
    pexps <- DataFrame(pexps, check.names = FALSE)
    pexps <- pexps[, rownames(colData(se))]
    pratios <- DataFrame(pratios, check.names = FALSE)
    pratios <- pratios[, rownames(colData(se))]
    
    assay(se, "resid") <- presids
    assay(se, "expected") <- pexps
    assay(se, "ratio") <- pratios

    ## add fits to metadata
    fits <- pfits$fit
    names(fits) <- pfits$condition
    metadata(se)$cy3models <- fits

    return(se)
}
