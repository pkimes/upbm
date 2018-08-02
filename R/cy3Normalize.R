#' Perform Cy3 Normalization
#'
#' Uses Cy3 intensities to quantify relative dsDNA material at each probe based
#' on observed vs. expected intensity ratios.
#' 
#' @param se SummarizedExperiment object containing GPR intensity information.
#' @param cy3se SummarizedExperiment object containing Cy3 GPR model output from
#'        \code{fitCy3Models}.
#' @param assay_name string name of the assay to normalize in \code{se}. (default = "fore")
#' @param match_by unquoted name of column in colData of SummarizedExperiments
#'        to use for matching samples across the two experiments; values of
#'        column must be unique for each sample in each experiment.
#'        (default = condition)
#' @param filter logical whether to filter "low quality" flagged probes by setting to NA.
#'        (default = TRUE)
#' @param scale logical whether to perform scaling by ratio of observed vs. expected
#'        Cy3 intensity at each probe. (default = FALSE)
#'
#' @return
#' SummarizedExperiment object.
#'
#' @seealso fitCy3Models
#' @importFrom dplyr as_tibble mutate select left_join
#' @importFrom tidyr gather spread
#' @export
#' @author Patrick Kimes
cy3Normalize <- function(se, cy3se, assay_name = "fore", match_by = condition, filter = TRUE, scale = FALSE) {

    stopifnot("ratio" %in% assayNames(cy3se))
    
    if (!filter & !scale) {
        return(se)
    }

    match_by <- rlang::enquo(match_by)
    match_by_str <- rlang::quo_name(match_by)

    ## check validity of matching colData column
    stopifnot(match_by_str %in% names(colData(se)))
    stopifnot(match_by_str %in% names(colData(cy3se)))
    match_vals1 <- colData(se)[[match_by_str]]
    match_vals2 <- colData(cy3se)[[match_by_str]]
    if (any(duplicated(match_vals2))) {
        stop("Matching variable '", match_by_str, "' is not unique across samples in 'cy3se'.\n",
             "Specify a different column in colData.")
    }

    match_overlap <- intersect(match_vals1, match_vals2)
    match_diff1 <- setdiff(match_vals1, match_vals2)
    match_diff2 <- setdiff(match_vals2, match_vals1)
    
    if (length(match_overlap) == 0) {
        stop("No samples matched across the specified two experiments.")
    }
    if (length(match_diff1) > 0) {
        warning("Dropping samples from se:\n  ", paste(match_diff1, collapse = ", "))
    }
    if (length(match_diff2) > 0) {
        warning("Dropping samples from cy3se:\n  ", paste(match_diff2, collapse = ", "))
    }

    ## filter samples
    se <- se[, colData(se)[[match_by_str]] %in% match_overlap, drop = FALSE]
    cy3se <- cy3se[, colData(cy3se)[[match_by_str]] %in% match_overlap, drop = FALSE]

    
    ## condition must be a unique column for faceting plot
    coldat1 <- data.frame(colData(se), check.names = FALSE,
                          check.rows = FALSE, stringsAsFactors = FALSE)
    coldat1 <- tibble::rownames_to_column(coldat1, "sample")
    coldat1 <- dplyr::mutate(coldat1, Match = rlang::UQ(match_by))
    coldat1 <- dplyr::select(coldat1, sample, Match)

    coldat2 <- data.frame(colData(cy3se), check.names = FALSE,
                          check.rows = FALSE, stringsAsFactors = FALSE)
    coldat2 <- tibble::rownames_to_column(coldat2, "sample")
    coldat2 <- dplyr::mutate(coldat2, Match = rlang::UQ(match_by))
    coldat2 <- dplyr::select(coldat2, sample, Match)

    
    new_assay <- as.data.frame(assay(se, assay_name), optional = TRUE)
    new_assay <- dplyr::as_tibble(new_assay)
    new_assay <- dplyr::mutate(new_assay,
                               Row = rowData(se)[, "Row"],
                               Column = rowData(se)[, "Column"])
    new_assay <- tidyr::gather(new_assay, sample, value, -Row, -Column)
    new_assay <- dplyr::left_join(new_assay, coldat1, by = "sample")
    ##new_assay <- dplyr::select(new_assay, -sample)

    
    if (scale) {
        scale_assay <- as.data.frame(assay(cy3se, "ratio"), optional = TRUE)
        scale_assay <- dplyr::as_tibble(scale_assay)
        scale_assay <- dplyr::mutate(scale_assay,
                                   Row = rowData(cy3se)[, "Row"],
                                   Column = rowData(cy3se)[, "Column"])
        scale_assay <- tidyr::gather(scale_assay, sample, scalar, -Row, -Column)
        scale_assay <- dplyr::left_join(scale_assay, coldat2, by = "sample")
        scale_assay <- dplyr::select(scale_assay, -sample)
        
        new_assay <- dplyr::left_join(new_assay, scale_assay, by = c("Row", "Column", "Match"))
        new_assay <- dplyr::mutate(new_assay, value = value / scalar)
        new_assay <- dplyr::select(new_assay, -scalar)
    }
        
    if (filter) {
        filter_assay <- as.data.frame(assay(cy3se, "lowq"), optional = TRUE)
        filter_assay <- dplyr::as_tibble(filter_assay)
        filter_assay <- dplyr::mutate(filter_assay,
                                      Row = rowData(cy3se)[, "Row"],
                                      Column = rowData(cy3se)[, "Column"])
        filter_assay <- tidyr::gather(filter_assay, sample, flag, -Row, -Column)
        filter_assay <- dplyr::left_join(filter_assay, coldat2, by = "sample")
        filter_assay <- dplyr::select(filter_assay, -sample)

        ## also filter anything NA in Cy3
        filter_assay$flag[is.na(filter_assay$flag)] <- FALSE
        
        new_assay <- dplyr::left_join(new_assay, filter_assay, by = c("Row", "Column", "Match"))
        new_assay <- dplyr::mutate(new_assay, value = ifelse(flag, NA, value))
        new_assay <- dplyr::select(new_assay, -flag)
    }

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

    ## replace assay in se object
    assay(se, assay_name) <- DataFrame(new_assay, check.names = FALSE)

    return(se)
}
