#' Lower-Tail Normalization (TESTING)
#' 
#' This is a collection of experimental cross-sample normalization procedures which
#' all attempt to align samples based on lower-tail of probe intensities (or log2-intensities).
#' Inherent to all of these cross-sample normalization approaches is the assumption that the
#' lower-tail probes of all samples are primarily background noise, and should therefore be
#' similiarly distributed, regardless of specific binding affinity or specificity for each
#' protein. These methods assume that the primary reason cross-sample normalization is needed is
#' due to differences in protein concentration, which can be corrected by scale and shift
#' transformation when concentration differences are small.
#' 
#' Currently, 3 appraoches for scale and shift transformation are supported: "regression",
#' "pca" and "normal."  The "regression" approach takes the lower \code{q}-quantile probes of the
#' reference sample and fits an ordinary least squares (OLS) model for each non-reference sample (Y)
#' on the reference sample (X). The fitted intercept and slope parameters are used as the shift and
#' scaling parameters. The "pca" approach applies principal components analysis (PCA) to the lower-tail
#' probes for each pair of the reference and a non-reference sample, and the first PC is used as the
#' scaling parameter. The shift parameter is calculated such that the first PC passes through the
#' means of both samples. For the "normal" method, a background (normal) plus signal is assumed
#' for probe-level intensities, and the normal background mean is estimated using
#' the \code{bg.parameters} function from the \code{affy} package. The corresponding
#' standard deviation is estimated separately from the lower tail for each sample by
#' fitting a truncated normal distribution to the tail.
#' 
#' \emph{This code is still under testing and not a finalized normalization procedure}.
#'
#' @param se SummarizedExperiment object containing GPR intensity information.
#' @param assay_name string name of the assay to normalize. (default = "fore")
#' @param q percentile between 0 and 1 specifying lower tail to use for
#'        fitting truncated normal distribution. (default = 0.4)
#' @param stratify unquoted name of column in colData of SummarizedExperiment (or
#'        '\code{sample}') to use for comparing samples; values in column must be
#'        unique for each sample. (default = condition)
#' @param baseline string name of baseline condition to compare other conditions
#'        against; if not specified, guessed by looking for
#'        'ref' in any value of the stratifying variable. (default = NULL)
#' @param log_scale logical whether cross-sample normalization should be applied
#'        using log2-scaled intensities rather than raw intensities. (default = FALSE)
#' @param method character string specifying method to use for normalization.
#'        Must be one of "regression", "pca", or "normal". (default = "regression")
#' @param .fits logical whether to just return a table of fits rather
#'        than the normalized SummarizedExperiment object. (default = FALSE)
#'
#' @return
#' SummarizedExperiment object with normalized intensities in
#' new assay, \code{scaled}.
#' 
#' @import truncnorm fitdistrplus
#' @importFrom stats quantile
#' @importFrom dplyr as_tibble tibble mutate group_by filter do select ungroup left_join
#' @importFrom tidyr gather spread nest
#' @importFrom rlang enquo quo_name
#' @importFrom affy bg.parameters
#' @export
#' @author Dongyuan Song, Patrick Kimes    
lowertailNormalization <- function(se, assay_name = "fore", q = 0.4, stratify = condition,
                                   baseline = NULL, log_scale = FALSE, 
                                   method = c("regression", "pca", "normal"),
                                   regtype = 1L, .fits = FALSE) {
    stopifnot(q > 0, q < 1)
    stopifnot(assay_name %in% assayNames(se))
    match.arg(method)
    
    ## check normalization stratification settings
    stratify <- rlang::enquo(stratify)
    strats <- .pbmCheckStratify(se, stratify, baseline)
    coldat <- strats$coldat
    baseline <- strats$baseline
    
    new_assay <- as.data.frame(assay(se, assay_name), optional = TRUE)
    new_assay <- dplyr::as_tibble(new_assay)
    new_assay <- dplyr::mutate(new_assay,
                               Row = rowData(se)[, "Row"],
                               Column = rowData(se)[, "Column"])
    new_assay <- tidyr::gather(new_assay, sample, value, -Row, -Column)
    new_assay <- dplyr::left_join(new_assay, coldat, by = "sample")
    new_assay <- dplyr::group_by(new_assay, sample, Stratify)
    new_assay <- dplyr::mutate(new_assay, ul = quantile(value, probs = q, na.rm = TRUE))
    new_assay <- dplyr::ungroup(new_assay)

    assay_fits <- dplyr::filter(new_assay, !is.na(value), value > 0)
    if (log_scale) {
        assay_fits <- dplyr::mutate(assay_fits, value = log2(value), ul = log2(ul))
    }
    
    if (method == "regression" || method == "pca") {
        bl_assay <- dplyr::filter(assay_fits, Stratify == baseline, value < ul)
        bl_assay <- dplyr::select(bl_assay, Row, Column, value)
        assay_fits <- dplyr::left_join(assay_fits, bl_assay, by = c("Row", "Column"),
                                       suffix = c("", ".bl"))
        assay_fits <- dplyr::filter(assay_fits, !is.na(value.bl))
        ## if (.filter_both) {
        ##     assay_fits <- dplyr::filter(assay_fits, value < ul)
        ## }
        assay_fits <- tidyr::nest(assay_fits, -sample, -Stratify, -ul)
        assay_ref <- dplyr::filter(assay_fits, Stratify == baseline)
        assay_fits <- dplyr::filter(assay_fits, Stratify != baseline)
        if (method == "regression") {
            assay_fits <- dplyr::mutate(assay_fits,
                                        fits = lapply(data, function(x) lm(value ~ value.bl, data = x)),
                                        est_shift = sapply(fits, function(x) { coef(x)[1] }),
                                        est_scale = sapply(fits, function(x) { coef(x)[2] }),
                                        r2adj = sapply(fits, function(x) { summary(x)$adj.r.squared }))
        } else if (method == "pca") {
            assay_fits <- dplyr::mutate(assay_fits,
                                        fits = lapply(data, function(x) { prcomp(as.matrix(x[, c("value", "value.bl")])) }),
                                        est_scale = sapply(fits, function(x) { x$rotation[1, 1] / x$rotation[2, 1] }),
                                        est_shift = mapply(function(x, y) { x$center[1] - y * x$center[2] },
                                                           x = fits, y = est_scale),
                                        varexpl = sapply(fits, function(x) { x$sdev[1]^2 / sum(x$sdev^2) }))
        }
        assay_fits <- dplyr::bind_rows(assay_fits, dplyr::mutate(assay_ref, est_shift = 0L, est_scale = 1L))
        assay_fits <- dplyr::select(assay_fits, -data)
    } else if (method == "normal") {
        assay_fits <- tidyr::nest(assay_fits, -sample, -Stratify, -ul)
        assay_fits <- dplyr::mutate(assay_fits,
                                    est_shift = sapply(data, function(x) affy::bg.parameters(x$value)$mu),
                                    data = mapply(function(x, y) { x[x$value <= y, ] }, data, ul, SIMPLIFY = FALSE),
                                    re_fit = mapply(function(x, y, z) { .fit_tnorm(x$value, y, z) },
                                                    data, ul, est_shift, SIMPLIFY = FALSE),
                                    est_scale = sapply(re_fit, function(x) x$estimate["sd"]))
        assay_fits <- dplyr::select(assay_fits, -re_fit)
    }

    if (.fits) {
        return(assay_fits)
    }
    
    ## reference params for scaling
    ref_shift <- assay_fits$est_shift[assay_fits$Stratify == baseline]
    ref_scale <- assay_fits$est_scale[assay_fits$Stratify == baseline]

    ## adjust to reference
    new_assay <- dplyr::left_join(new_assay, assay_fits, by = c("sample", "Stratify"))
    if (log_scale) {
        new_assay <- dplyr::mutate(new_assay, value = ifelse(value > 0, value, NA))
        new_assay <- dplyr::mutate(new_assay, value = log2(value))
    }
    new_assay <- dplyr::mutate(new_assay, value = (value - est_shift) / est_scale)
    new_assay <- dplyr::mutate(new_assay, value = ref_shift + value * ref_scale)
    if (log_scale) {
        new_assay <- dplyr::mutate(new_assay, value = 2^value)
    }
    
    ## return to square assay shape
    new_assay <- dplyr::select(new_assay, sample, value, Row, Column)
    new_assay <- tidyr::spread(new_assay, sample, value)

    ## match row order to rowData
    c_order <- paste(rowData(se)$Row, rowData(se)$Column, sep = "-")
    new_order <- match(paste(new_assay$Row, new_assay$Column, sep = "-"), c_order)
    stopifnot(!duplicated(new_order), length(new_order) == nrow(se))
    new_assay <- new_assay[new_order, ]
    new_assay <- dplyr::select(new_assay, -Row, -Column)

    ## match column order to colData
    stopifnot(colnames(new_assay) %in% colnames(se))
    new_assay <- new_assay[, colnames(se)]

    ## add new assay to se object
    assay(se, "scaled") <- DataFrame(new_assay, check.names = FALSE)

    ## store reference mean, sd information
    metadata(se)$ref_shift <- ref_shift
    metadata(se)$ref_scale <- ref_scale
    
    return(se)
}

.fit_tnorm <- function(vals, upper, vbar) {
    fitdist(vals, distr = dtruncnorm,
            fix.arg = list("a" = -20, "b" = upper, "mean" = vbar), 
            start = list(sd = sd(vals)))
}
