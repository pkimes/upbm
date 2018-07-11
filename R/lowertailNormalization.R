#' Lower-Tail Normalization (TESTING)
#' 
#' This is an experimental normalization procedure based on fitting
#' a truncated normal distribution to the lower tail of log2-scale probe
#' intensities to register intensities across experiments. The mean of
#' the distribution is estimated using the \code{bg.parameters} function
#' from the \code{affy} package, and only the standard deviation is estimated
#' separately for the lower tail. The intensities are shift/scale transformed
#' on the log2-scale to match a specified reference sample.
#' \emph{This code is still under testing and not a finalized normalization procedure}.
#'
#' @param se SummarizedExperiment object containing GPR intensity information.
#' @param assay_name string name of the assay to normalize. (default = "fore")
#' @param q percentile between 0 and 1 specifying lower tail to use for
#'        fitting truncated normal distribution. (default = 0.4)
#' @param stratify unquoted name of column in colData of SummarizedExperiment (or
#'        '\code{sample}') to use for comparing samples; values in column must be
#'        unique for each sample. (default = condition)
#' @param baseline string name of baseline condition to
#'        compare other conditions against; if not specified, guessed by looking for
#'        'ref' in any value of the stratifying variable. (default = NULL)
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
#' @author Patrick Kimes
lowertailNormalization <- function(se, assay_name = "fore", q = 0.4, stratify = condition,
                                   baseline = NULL, .fits = FALSE) {
    stopifnot(q > 0, q < 1)
    stopifnot(assay_name %in% assayNames(se))

    ## check normalization stratification settings
    if (!.fits) {
        stratify <- rlang::enquo(stratify)
        strats <- .pbmCheckStratify(se, stratify, baseline)
        coldat <- strats$coldat
        baseline <- strats$baseline
    } else {
        coldat <- dplyr::tibble(sample = colnames(se),
                                Stratify = colnames(se))
    }
    
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
    assay_fits <- dplyr::mutate(assay_fits, value = log2(value), ul = log2(ul))
    assay_fits <- tidyr::nest(assay_fits, -sample, -Stratify, -ul)

    assay_fits <- dplyr::mutate(assay_fits,
                                rma_fit = lapply(data, function(x) affy::bg.parameters(x$value)),
                                rma_mean = sapply(rma_fit, function(x) x$mu),
                                rma_sd = sapply(rma_fit, function(x) x$sigma))
    assay_fits <- dplyr::mutate(assay_fits, data = mapply(function(x, y) { x[x$value <= y, ] }, data, ul, SIMPLIFY = FALSE))
    assay_fits <- dplyr::mutate(assay_fits,
                                re_fit = mapply(function(x, y, z) { .fit_tnorm(x$value, y, z) },
                                                data, ul, rma_mean, SIMPLIFY = FALSE),
                                est_mean = rma_mean,
                                est_sd = sapply(re_fit, function(x) x$estimate["sd"]))
    
    if (.fits) {
        return(assay_fits)
    }

    ## reference params for scaling
    ref_mean <- assay_fits$est_mean[assay_fits$Stratify == baseline]
    ref_sd <- assay_fits$est_sd[assay_fits$Stratify == baseline]

    ## adjust to reference
    assay_fits <- dplyr::select(assay_fits, -rma_fit, -re_fit)
    new_assay <- dplyr::left_join(new_assay, assay_fits, by = c("sample", "Stratify"))
    new_assay <- dplyr::mutate(new_assay, value = ifelse(value > 0, value, NA))
    new_assay <- dplyr::mutate(new_assay, value = log2(value))
    new_assay <- dplyr::mutate(new_assay, value = (value - est_mean) / est_sd)
    new_assay <- dplyr::mutate(new_assay, value = ref_mean + value * ref_sd)
    new_assay <- dplyr::mutate(new_assay, value = 2^value)

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
    metadata(se)$ref_mean <- ref_mean
    metadata(se)$ref_sd <- ref_sd
    
    return(se)
}

.fit_tnorm <- function(vals, upper, vbar) {
    fitdist(vals, distr = dtruncnorm,
            fix.arg = list("a" = -20, "b" = upper, "mean" = vbar), 
            start = list(sd = sd(vals)))
}
