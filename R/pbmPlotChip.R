#' Plot of Intensities on Grid
#'
#' Probe-level intensities plotted by row and column indices
#' defined in the \code{rowData} of the SummarizedExperiment
#' object. This can be useful for visually inspecting spatial
#' trends and biases across the scanned chip.
#' 
#' @param se SummarizedExperiment object containing GPR
#'        intensity information.
#' @param assay string name of the assay to plot.
#'        (default = \code{SummarizedExperiment::assayNames(se)[1]})
#' @param log_scale logical whether to plot the intensities
#'        on the log-scale. (default = TRUE)
#' @param relative_scale string name of column to in \code{colData}
#'        or numeric vector to use for scaling intensities in each
#'        sample. This can be useful for plotting spatial bias of
#'        samples relative to the total intensity in the sample. 
#'        Ignored if NULL. (default = NULL)
#' @param bound non-negative numeric value to use for bounding outliers.
#'        When specified, values are bound to \code{(-bound, bound)}.
#'        If set to 0L, results are binarized to positive, negative.
#'        Ignored if \code{Inf}. (default = Inf)
#' @param .facet logical whether plot should be faceted using
#'        the default 'condition' column from the colData of the
#'        SummarizedExperiment. (default = TRUE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to plotting. (default = 1)
#' 
#' @return
#' ggplot object
#'
#' @details
#' If the values in the specified assay include negative
#' values, \code{log_scale = FALSE} will be used. Note that
#' this behavior is different from the other plotting functions
#' in this package. To get around this, manually set all negative
#' values to NA prior to plotting.
#'
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom tidyr gather
#' @importFrom dplyr left_join
#' @importFrom rlang sym
#' @import ggplot2
#' @importFrom SummarizedExperiment assay assayNames rowData colData
#' @export
#' @author Patrick Kimes
pbmPlotChip <- function(se, assay = SummarizedExperiment::assayNames(se)[1],
                        log_scale = TRUE, relative_scale = NULL,
                        bound = Inf, .facet = TRUE) {
    stopifnot(assay %in% SummarizedExperiment::assayNames(se))
    stopifnot("Row" %in% names(rowData(se)))
    stopifnot("Column" %in% names(rowData(se)))

    ## condition must be a unique column for faceting plot
    if (.facet && any(duplicated(colData(se)$condition))) {
        stop("'condition' column of colData must be unique.")
    }

    ## filter probes
    se <- pbmFilterProbes(se) 

    ## condition must be a unique column for faceting plot
    coldat <- data.frame(colData(se), check.names = FALSE,
                         check.rows = FALSE, stringsAsFactors = FALSE)
    coldat <- tibble::rownames_to_column(coldat, "sample")

    ## extract intensities
    pdat <- SummarizedExperiment::assay(se, assay)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- tibble::as_tibble(pdat)
    pdat <- dplyr::mutate(pdat,
                          Row = rowData(se)[, "Row"],
                          Column = rowData(se)[, "Column"])
    pdat <- tidyr::gather(pdat, sample, value, -Column, -Row)
    pdat <- dplyr::left_join(pdat, coldat, by = "sample")

    ## apply scaling and bounding
    if (!is.null(relative_scale)) {
        if (length(relative_scale) != 1 || !(relative_scale %in% names(coldat))) {
            warning("'relative_scale' should be length = 1 column name in colData.\n",
                    "valid names: ", paste(names(coldat), collapse = ", "), ".")
            relative_scale <- NULL
        } else {
            pdat <- dplyr::mutate(pdat, value = value / !! rlang::sym(relative_scale))
        }
    }

    ## check for negative values
    pdat_negative <- (min(pdat$value, na.rm = TRUE) < 0)
    
    ## check if log-scaling is valid
    if (log_scale && !is.null(relative_scale)) {
        warning("Ignoring log-scaling since relative plotting also specified.")
        log_scale <- FALSE
    }
    if (log_scale && pdat_negative) {
        warning("Ignoring log-scaling since min value is negative.")
        log_scale <- FALSE
    }

    ## handle log-scale plotting if requested 
    if (log_scale) {
        ## direct-scale values since 'fill' doesn't accept 'trans='
        pdat <- dplyr::mutate(pdat, value = log2(value))
        ## re-check for negative values
        pdat_negative <- (min(pdat$value, na.rm = TRUE) < 0)
        ptitle <- "PBM Intensity (log2)"
    } else {
        ptitle <- "PBM Intensity"
    }

    ## bound values if 'bound' is finite
    if (bound == 0L) {
        pdat <- dplyr::mutate(pdat, value = ifelse(value > 0, 1, -1))
    } else if (bound < Inf) {
        pdat <- dplyr::mutate(pdat,
                              value = pmin(value, abs(bound)),
                              value = pmax(value, -abs(bound)))
    }
    
    ## handle fill-scales if negative values present 
    if (pdat_negative) {
        pfill <- scale_fill_gradient2(assay)
        if (bound == 0L) {
            pfill$limits <- c(-1, 1)
        } else if (bound < Inf) {
            pfill$limits <- c(-abs(bound), abs(bound))
        }
    } else {
        if (log_scale) {
            pfill <- scale_fill_distiller(assay, direction = 1,
                                          labels = function(x) {2^as.numeric(x)})
        } else {
            pfill <- scale_fill_distiller(assay, direction = 1)
        }
        if (bound < Inf) {
            pfill$limits <- c(0, abs(bound))
        }
    }
    
    gp <- ggplot(pdat, aes(x = Column, y = Row, fill = value)) +
        geom_tile() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw() +
        pfill + 
        ggtitle(ptitle)

    ## facet plot using default
    if (.facet) {
        gp <- gp + facet_grid(. ~ condition)
    }
    gp
}
