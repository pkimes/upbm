#' TMM Normalization
#' 
#' This is the function using trimmed mean of M-values (TMM) method
#' to perform cross sample normalization. TMM normalization is very popular
#' for RNA-seq data normalization and was included
#' in package \code{edgeR}. However, here a modified TMM methods is used.
#' The difference is that weights are not introduced when calculating trimmed
#' mean since unlike RNA-seq, PBM data shows almost constant variance of M-value across A-value.
#' Also, for PBM data there is no library size.
#' 
#' For TMM normalization in \code{upbm}, the first step is trmming A-values
#' (one side by \code{q} or two sides by \code{q} and \code{qb}). The second step
#' is trimming M-values.
#' @param se SummarizedExperiment object containing GPR intensity information.
#' @param assay string name of the assay to normalize.
#'        (default = \code{SummarizedExperiment::assayNames(se)[1]})
#' @param q percentile between 0 and 1 specifying lower tail A-values to use for
#'        normalization. (default = 0.6)
#' @param qb percentile either between 0 and 1-\code{q} or equals \code{Auto}. If \code{qb}
#'        is numeric then it specifies the proportion for lowest A-values that will
#'        be removed; if \code{qb = "Auto"} then this cutoff will be the A-value of
#'        the medians of two background intensities. (defalut = 0)
#' @param q0 percentile between 0 and 0.5 specifying the fraction of M-values to be trimmed
#'        from each side. If \code{q0 = 0.5} it yields median.  (default = 0.2)
#' @param stratify string name of column in colData of SummarizedExperiment to
#'        use for comparing samples; values in column must be
#'        unique for each sample. Alternatively, can specify '\code{"sample"}' to
#'        use column names. (default = "condition")
#' @param baseline string name of baseline condition to compare other conditions
#'        against; if not specified, guessed by looking for
#'        'ref' in any value of the stratifying variable. (default = NULL)
#' @param .fits logical whether to just return a table of fits rather
#'        than the normalized SummarizedExperiment object. (default = FALSE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to normalization. See \code{pbmFilterProbes}
#'        for more details on probe filter levels. (default = 1)
#' @return
#' SummarizedExperiment object with normalized intensities in
#' new assay, \code{scaled}.
#' 
#' @importFrom stats quantile
#' @importFrom dplyr as_tibble tibble mutate group_by filter do select ungroup left_join
#' @importFrom tidyr gather spread nest
#' @importFrom rlang enquo quo_name
#' @export
#' @author Dongyuan Song, Patrick Kimes
tmmNormalization <- function(se, assay = SummarizedExperiment::assayNames(se)[1],
                             q = 0.6, q0 = 0.2, qb = 0, stratify = "condition",
                             baseline = NULL,
                             .filter = TRUE,
                             .fits = FALSE) {
  stopifnot(q >= 0, q < 1)
  stopifnot(qb < 1 - q || qb == "Auto")
  stopifnot(q0 >= 0, q0 <= 0.5)
  stopifnot(assay %in% SummarizedExperiment::assayNames(se))

  ## filter probes - only for computing shift/scale factors (return original se)
  fse <- pbmFilterProbes(se, .filter)
  
  ## check normalization stratification settings
  strats <- .pbmCheckStratify(fse, stratify, baseline)
  coldat <- strats$coldat
  baseline <- strats$baseline
  
  ## tidy up data for computing factors
  scale_assay <- as.data.frame(assay(fse, assay), optional = TRUE)
  scale_assay <- dplyr::as_tibble(scale_assay)
  scale_assay <- dplyr::mutate(scale_assay,
                               Row = rowData(fse)[, "Row"],
                               Column = rowData(fse)[, "Column"])
  scale_assay <- tidyr::gather(scale_assay, sample, value, -Row, -Column)
  scale_assay <- dplyr::left_join(scale_assay, coldat, by = "sample")
  scale_assay <- dplyr::group_by(scale_assay, sample, Stratify)
  
  scale_assay <- dplyr::ungroup(scale_assay)
  
  assay_fits <- dplyr::filter(scale_assay, !is.na(value), value > 0)
  
  bl_assay <- dplyr::filter(assay_fits, Stratify == baseline)
  bl_assay <- dplyr::select(bl_assay, Row, Column, value)
  assay_fits <- dplyr::left_join(assay_fits, bl_assay, by = c("Row", "Column"),
                                 suffix = c("", ".bl"))
  assay_fits <- dplyr::filter(assay_fits, !is.na(value.bl))
  assay_fits <- dplyr::mutate(assay_fits,
                              M.value = (log2(value) - log2(value.bl)),
                              A.value = (log2(value) + log2(value.bl))/2)  
  
  assay_fits <- tidyr::nest(assay_fits, -sample, -Stratify)
  assay_ref <- dplyr::filter(assay_fits, Stratify == baseline)
  assay_fits <- dplyr::filter(assay_fits, Stratify != baseline)
  
  ## extract background intensities
  bdat <- assay(se, "back")
  bdat <- as.data.frame(bdat, optional = TRUE)
  bdat <- tibble::as_tibble(bdat)
  bdat <- dplyr::mutate(bdat,
                        Row = rowData(se)[, "Row"],
                        Column = rowData(se)[, "Column"])
  bdat <- tidyr::gather(bdat, sample, value, -Column, -Row)
  bdat <- dplyr::left_join(bdat, coldat, by = "sample")
  bdat <- dplyr::select(bdat, -sample)
  
  ## manipulate to get reference condition in separate column
  bdat <- tidyr::spread(bdat, Stratify, value)
  
  b_ref <- dplyr::select(bdat, !!baseline)
  b_ref <- median(as.matrix(b_ref), na.rm = TRUE)
  
  bdat <- dplyr::rename(bdat, Baseline = !!baseline)
  bdat <- tidyr::gather(bdat, Stratify, value, -Baseline, -Row, -Column)
  
  bdat <- dplyr::group_by(bdat, Stratify)
  bdat <- dplyr::summarise(bdat, median = median(value, na.rm = TRUE))
  
  bdat <- dplyr::mutate(bdat, A_b = (log2(median) + log2(b_ref)) / 2)
  bdat <- dplyr::select(bdat, -median)
  
  ## Join bdat with assay_fits
  assay_fits <- dplyr::left_join(assay_fits, bdat, by = "Stratify")
  
  ## TMM method
  assay_fits <- dplyr::mutate(assay_fits,
                              est_shift = 0L,
                              est_scale = mapply(function(x, y) {
                                  z <- dplyr::filter(x, A.value < quantile(A.value, q))
                                  
                                  if (qb == "Auto") z <- dplyr::filter(z, A.value > y)
                                  else z <- dplyr::filter(z, A.value > quantile(A.value, qb))
                                  
                                  z <- dplyr::summarise(z, est_scale = 2^mean(M.value, trim = q0, na.rm = TRUE))
                                  z <- as.numeric(z)
                                  return(z)
                              }, data, A_b)
                              )
  
  assay_fits <- dplyr::bind_rows(assay_fits, dplyr::mutate(assay_ref, est_shift = 0L, est_scale = 1L))
  assay_fits <- dplyr::select(assay_fits, -data)
  
  if (.fits) {
    return(assay_fits)
  }
  
  ## reference params for scaling
  ref_shift <- assay_fits$est_shift[assay_fits$Stratify == baseline]
  ref_scale <- assay_fits$est_scale[assay_fits$Stratify == baseline]
  
  ## tidy up original data
  new_assay <- as.data.frame(assay(se, assay), optional = TRUE)
  new_assay <- dplyr::as_tibble(new_assay)
  new_assay <- dplyr::mutate(new_assay,
                             Row = rowData(se)[, "Row"],
                             Column = rowData(se)[, "Column"])
  new_assay <- tidyr::gather(new_assay, sample, value, -Row, -Column)
  new_assay <- dplyr::left_join(new_assay, coldat, by = "sample")
  
  ## adjust to reference
  new_assay <- dplyr::left_join(new_assay, assay_fits, by = c("sample", "Stratify"))
  
  new_assay <- dplyr::mutate(new_assay, value = (value - est_shift) / est_scale)
  new_assay <- dplyr::mutate(new_assay, value = ref_shift + value * ref_scale)
  
  ## return to square assay shape
  new_assay <- dplyr::select(new_assay, sample, value, Row, Column)
  new_assay <- tidyr::spread(new_assay, sample, value)
  
  ## match row order to rowData
  c_order <- paste(rowData(se)$Row, rowData(se)$Column, sep = "-")
  new_order <- match(c_order, paste(new_assay$Row, new_assay$Column, sep = "-"))
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
  
  ## store scaling parameters
  assay_fits <- dplyr::select(assay_fits, -Stratify)
  if ("fits" %in% names(assay_fits)) {
    assay_fits <- dplyr::select(assay_fits, -fits)
  }
  coldat <- merge(colData(se), data.frame(assay_fits, row.names = "sample"),
                  by = 0, all = TRUE)
  rownames(coldat) <- coldat$Row.names
  coldat$Row.names <- NULL
  
  ## match colData row order with SE col order
  coldat <- coldat[match(colnames(se), rownames(coldat)), , drop = FALSE]
  
  stopifnot(all(rownames(coldat) == colnames(se)))
  colData(se) <- coldat
  
  return(se)
}
