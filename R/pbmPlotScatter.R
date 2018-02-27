#' Scatterplot of Intensities
#'
#' Probe-level intensities plotted pairwise against a reference
#' condition. Given a baseline condition specified with \code{baseline=}
#' (which must map uniquely in the \code{stratify=} column of the colData),
#' pairwise probe intensities are plotted with the identity line (x = y).
#'
#' Alternatively, if \code{maplot} is TRUE, a MA plot is created with the
#' same data, plotting A (log-scale mean intensity) against M (log-scale
#' difference in intensities).
#' 
#' @param se SummarizedExperiment object containing GPR
#'        intensity information.
#' @param assay_name string name of the assay to plot.
#'        (default = "fore")
#' @param stratify unquoted name of column in colData of SummarizedExperiment (or
#'        '\code{sample}') to use for comparing samples; values in column must be
#'        unique for each sample. (default = condition)
#' @param baseline string name of baseline condition to
#'        compare other conditions against; if not specified, guessed by looking for
#'        'ref' in any value of the stratifying variable. (default = NULL)
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
#' MA plots are drawn on the log scale by definition.
#'
#' @importFrom tidyr spread gather
#' @importFrom dplyr mutate select left_join rename
#' @importFrom rlang enquo quo_name UQ 
#' @import ggplot2 SummarizedExperiment
#' @export
#' @author Patrick Kimes
pbmPlotScatter <- function(se, assay_name = "fore", stratify = condition, baseline = NULL,
                           log_scale = TRUE, maplot = FALSE, .method = "auto", .filter = 1) {
    stopifnot(assay_name %in% assayNames(se))

    if (! "Row" %in% names(rowData(se)) || ! "Column" %in% names(rowData(se))) {
        if ("kmer" %in% names(rowData(se))) {
            rowData(se)$Row <- rowData(se)$kmer
            rowData(se)$Column <- rowData(se)$kmer
        } else {
            stop("rowData must either contain 'Row','Column' columns or 'kmer' column.")
        }
    }
        
    stratify <- rlang::enquo(stratify)
    stratify_str <- rlang::quo_name(stratify)
    
    ## check validity of stratifying colData column
    stopifnot(stratify_str %in% names(colData(se)))
    strat_vals <- colData(se)[[stratify_str]]
    if (any(duplicated(strat_vals))) {
        stop("Stratifying variable '", stratify_str, "' is not unique across samples.\n",
             "Specify a different column in colData.")
    }
        
    ## determine baseline if necessary; check validity
    if (is.null(baseline)) {
        baseline <- grep("ref", strat_vals, value = TRUE, ignore.case = TRUE)
        if (length(baseline) > 1) {
            stop("Too many candidate baseline states in '", stratify, "' column: ",
                 paste0(baseline, collapse = ", "), ".\n",
                 "Specify correct baseline condition w/ 'baseline'.")
        }
    } else {
        if (! baseline %in% strat_vals) {
            stop(baseline, " is not a value in '", stratify, "' column.\n",
                 "Specify correct baseline condition w/ 'baseline'.")
        }
    } 
    
    ## filter probes
    se <- pbmFilterProbes(se, assay_name, .filter) 

    ## condition must be a unique column for faceting plot
    coldat <- as.data.frame(colData(se))
    coldat <- tibble::rownames_to_column(coldat, "sample")
    coldat <- dplyr::mutate(coldat, Stratify = rlang::UQ(stratify))
    coldat <- dplyr::select(coldat, sample, Stratify)
        
    ## extract intensities
    pdat <- cbind(assay(se, assay_name), rowData(se)[, c("Column", "Row")])
    pdat <- data.frame(pdat, stringsAsFactors = FALSE)
    pdat <- tibble::as_tibble(pdat)
    pdat <- tidyr::gather(pdat, sample, value, -Column, -Row)
    pdat <- dplyr::left_join(pdat, coldat, by = "sample")
    pdat <- dplyr::select(pdat, -sample)
    
    ## check for negative values
    pdat_pnegative <- mean(pdat$value < 0, na.rm = TRUE)

    ## warn user about dropped values
    if (log_scale && pdat_pnegative > 0) {
        warning("Log-scaling will drop ", round(pdat_pnegative * 100, 2),
                "% of values.")
    }
    
    ## handle log-scale plotting if requested
    if (maplot) {
        ptitle <- paste0("PBM Intensity MA Plot (log2; reference = '", baseline, "')")
        pxaxis <- scale_x_continuous("A; mean (condition + reference) / 2")
        pyaxis <- scale_y_continuous("M; difference (condition - reference)")
        pline <- geom_hline(yintercept = 0, color = 'dodgerblue3')
    } else {
        if (log_scale) {
            ptitle <- "PBM Intensity (log2)"
            pxaxis <- scale_x_continuous(paste0("reference (", baseline, ") intensity (log)"),
                                         breaks = 2^(0:100), trans = "log10")
            pyaxis <- scale_y_continuous("intensity (log)", breaks = 2^(0:100), trans = "log10")
        } else {
            ptitle <- "PBM Intensity"
            pxaxis <- scale_x_continuous(paste0("reference (", baseline, ") intensity"))
            pyaxis <- scale_y_continuous("intensity")
        }
        pline <- geom_abline(intercept = 0, slope = 1, color = 'dodgerblue3')
    }
    
    ## manipulate to get reference condition in separate column
    pdat <- tidyr::spread(pdat, Stratify, value)
    pdat <- dplyr::select(pdat, -Row, -Column)
    pdat <- dplyr::rename(pdat, Baseline = rlang::UQ(baseline))
    pdat <- tidyr::gather(pdat, Stratify, value, -Baseline)

    ## calculate MA values if necessary and create base of plot
    if (maplot) {
        pdat <- dplyr::mutate(pdat,
                              M = log2(value / Baseline),
                              A = .5*log2(value * Baseline))
        gp <- ggplot(pdat, aes(x = A, y = M)) +
            theme_bw()
    } else {
        gp <- ggplot(pdat, aes(x = Baseline, y = value)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 1/2, hjust = 1))
    }

    ## add common parts of plot
    gp <- gp + 
        geom_point(alpha = 1/10) +
        expand_limits(y = 0) +
        facet_grid(. ~ Stratify, scales = "free_x") +
        pxaxis +
        pyaxis +
        pline + 
        ggtitle(ptitle)

    if (!is.null(.method)) {
        gp <- gp + geom_smooth(method = .method, se = FALSE)
    }

    gp
}
