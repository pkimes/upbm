#' Paired Scatterplot of Intensities
#'
#' Probe-level intensities plotted pairwise between two replicates which
#' are matched using the \code{match_by} column of the colData. Only pairs of
#' samples which are present in both experiments are plotted. 
#'
#' As with \code{pbmPlotScatter}, if \code{maplot} is TRUE, a MA plot is created
#' with the same data, plotting A (log-scale mean intensity) against M (log-scale
#' difference in intensities).
#' 
#' @param se1 first SummarizedExperiment object containing GPR
#'        intensity information.
#' @param se2 second SummarizedExperiment object containing GPR
#'        intensity information to be compared with first experiment.
#' @param assay_name string name of the assay to plot.
#'        (default = "fore")
#' @param match_by unquoted name of column in colData of SummarizedExperiments
#'        to use for matching samples across the two experiments; values of
#'        column must be unique for each sample in each experiment.
#'        (default = condition)
#' @param log_scale logical whether to plot the intensities
#'        on the log-scale. (default = TRUE)
#' @param maplot logical whether to plot MA plot rather than standard
#'        scatter plot. (default = FALSE)
#' @param .method value passed to \code{method = } parameter of the
#'        \code{ggplot2::geom_smooth} function for adding a smoothed
#'        fit to each scatter plot; to prevent any line, set NULL.
#'        (default = "auto")
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to plotting. (default = 1)
#'
#' @return
#' ggplot object
#'
#' @details
#' Currently, function only includes capability of plotting
#' pairwise comparison of two experiment sets. This was initially
#' implemented with the understanding that rarely more than
#' two replicates of a set of experiments are performed with
#' the PBM technology. If the need becomes clear, this function
#' will be extended to accept an arbitrary number of
#' SummarizedExperiment objects.
#'
#' @importFrom tidyr spread gather
#' @importFrom dplyr mutate select left_join rename bind_rows
#' @importFrom rlang enquo quo_name UQ
#' @import ggplot2 SummarizedExperiment
#' @export
#' @author Patrick Kimes
pbmPlotComparison <- function(se1, se2, assay_name  = "fore", match_by = condition,
                              log_scale = TRUE,  maplot = FALSE,
                              .method = "auto", .filter = 1) {
    stopifnot(assay_name %in% assayNames(se1))
    if (! "Row" %in% names(rowData(se1)) || ! "Column" %in% names(rowData(se1))) {
        if ("kmer" %in% names(rowData(se1))) {
            rowData(se1)$Row <- rowData(se1)$kmer
            rowData(se1)$Column <- rowData(se1)$kmer
        } else {
            stop("rowData must either contain 'Row','Column' columns or 'kmer' column.")
        }
    }

    stopifnot(assay_name %in% assayNames(se2))
    if (! "Row" %in% names(rowData(se2)) || ! "Column" %in% names(rowData(se2))) {
        if ("kmer" %in% names(rowData(se2))) {
            rowData(se2)$Row <- rowData(se2)$kmer
            rowData(se2)$Column <- rowData(se2)$kmer
        } else {
            stop("rowData must either contain 'Row','Column' columns or 'kmer' column.")
        }
    }
    
    match_by <- rlang::enquo(match_by)
    match_by_str <- rlang::quo_name(match_by)
    
    ## check validity of matching colData column
    stopifnot(match_by_str %in% names(colData(se1)))
    stopifnot(match_by_str %in% names(colData(se2)))
    match_vals1 <- colData(se1)[[match_by_str]]
    match_vals2 <- colData(se2)[[match_by_str]]
    if (any(duplicated(match_vals1))| any(duplicated(match_vals2))) {
        stop("Matching variable '", match_by_str, "' is not unique across samples.\n",
             "Specify a different column in colData.")
    }

    match_overlap <- intersect(match_vals1, match_vals2)
    match_diff1 <- setdiff(match_vals1, match_vals2)
    match_diff2 <- setdiff(match_vals2, match_vals1)
    
    if (length(match_overlap) == 0) {
        stop("No samples matched across the specified two experiments.")
    }
    if (length(match_diff1) > 0) {
        warning("Dropping samples from se1:\n  ", paste(match_diff1, collapse = ", "))
    }
    if (length(match_diff2) > 0) {
        warning("Dropping samples from se2:\n  ", paste(match_diff2, collapse = ", "))
    }

    ## filter samples
    se1 <- se1[, colData(se1)[[match_by_str]] %in% match_overlap, drop = FALSE]
    se2 <- se2[, colData(se2)[[match_by_str]] %in% match_overlap, drop = FALSE]
    
    ## filter probes
    se1 <- pbmFilterProbes(se1, .filter) 
    se2 <- pbmFilterProbes(se2, .filter) 

    ## condition must be a unique column for faceting plot
    coldat1 <- data.frame(colData(se1), check.names = FALSE,
                          check.rows = FALSE, stringsAsFactors = FALSE)
    coldat1 <- tibble::rownames_to_column(coldat1, "sample")
    coldat1 <- dplyr::mutate(coldat1, Match = rlang::UQ(match_by))
    coldat1 <- dplyr::select(coldat1, sample, Match)

    coldat2 <- data.frame(colData(se2), check.names = FALSE,
                          check.rows = FALSE, stringsAsFactors = FALSE)
    coldat2 <- tibble::rownames_to_column(coldat2, "sample")
    coldat2 <- dplyr::mutate(coldat2, Match = rlang::UQ(match_by))
    coldat2 <- dplyr::select(coldat2, sample, Match)
    
    ## extract intensities
    pdat1 <- assay(se1, assay_name)
    pdat1 <- as.data.frame(pdat1, optional = TRUE)
    pdat1 <- tibble::as_tibble(pdat1)
    pdat1 <- dplyr::mutate(pdat1,
                           Row = rowData(se1)[, "Row"],
                           Column = rowData(se1)[, "Column"])
    pdat1 <- tidyr::gather(pdat1, sample, value, -Column, -Row)
    pdat1 <- dplyr::left_join(pdat1, coldat1, by = "sample")
    pdat1 <- dplyr::select(pdat1, -sample)

    pdat2 <- assay(se2, assay_name)
    pdat2 <- as.data.frame(pdat2, optional = TRUE)
    pdat2 <- tibble::as_tibble(pdat2)
    pdat2 <- dplyr::mutate(pdat2,
                           Row = rowData(se2)[, "Row"],
                           Column = rowData(se2)[, "Column"])
    pdat2 <- tidyr::gather(pdat2, sample, value, -Column, -Row)
    pdat2 <- dplyr::left_join(pdat2, coldat2, by = "sample")
    pdat2 <- dplyr::select(pdat2, -sample)
    
    ## check for negative values
    pdat1_pnegative <- mean(pdat1$value < 0, na.rm = TRUE)
    pdat2_pnegative <- mean(pdat2$value < 0, na.rm = TRUE)

    ## warn user about dropped values
    if (log_scale && pdat1_pnegative + pdat2_pnegative > 0) {
        warning("Log-scaling will drop ",
                round(pdat1_pnegative * 100, 2), "% and ",
                round(pdat2_pnegative * 100, 2), "% of values.")
    }
    
    ## handle log-scale plotting if requested
    if (maplot) {
        if (log_scale) {
            ptitle <- paste0("PBM Intensity Comparison MA Plot (log2)")
        } else {
            warning("Note that MA plots are typically drawn on the log2 scale.\n",
                    "Consider re-plotting with log_scale = TRUE if the original ",
                    "values were not on the log2 scale.")
            ptitle <- paste0("PBM Intensity Comparison MA Plot")
        }
        pxaxis <- scale_x_continuous("A; mean (rep1 + rep2) / 2")
        pyaxis <- scale_y_continuous("M; difference (rep1 - rep2)")
        pline <- geom_hline(yintercept = 0, color = 'dodgerblue3')
    } else {
        if (log_scale) {
            ptitle <- "PBM Intensity Comparison (log2)"
            pxaxis <- scale_x_continuous("rep1 intensity (log)", breaks = 2^(0:100), trans = "log10")
            pyaxis <- scale_y_continuous("rep2 intensity (log)", breaks = 2^(0:100), trans = "log10")
        } else {
            ptitle <- "PBM Intensity Comparison"
            pxaxis <- scale_x_continuous("rep1 intensity")
            pyaxis <- scale_y_continuous("rep2 intensity")
        }
        pline <- geom_abline(intercept = 0, slope = 1, color = 'dodgerblue3')
    }

    ## Spread gather to get reference condition in separate column
    pdat <- dplyr::bind_rows(list(rep1 = pdat1, rep2 = pdat2), .id = "replicate")
    pdat <- tidyr::spread(pdat, replicate, value)

    ## calculate MA values if necessary and create base of plot
    if (maplot) {
        if (log_scale) {
            pdat <- dplyr::mutate(pdat,
                                  rep1 = log2(rep1),
                                  rep2 = log2(rep2))
        }
        pdat <- dplyr::mutate(pdat,
                              M = rep1 - rep2,
                              A = (rep1 + rep2) / 2)
        gp <- ggplot(pdat, aes(x = A, y = M)) +
            theme_bw()
    } else {
        gp <- ggplot(pdat, aes(x = rep1, y = rep2)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 1/2, hjust = 1))
    }
    
    ## add common parts of plot
    gp <- gp +
        geom_point(alpha = 1/10) +
        expand_limits(y = 0) +
        facet_grid(. ~ Match, scales = "free_x") +
        pxaxis +
        pyaxis +
        pline + 
        ggtitle(ptitle)

    if (!is.null(.method)) {
        gp <- gp + geom_smooth(method = .method, se = FALSE)
    }

    gp
}
