#' Scatterplot of Intensities
#'
#' @description
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
#' @param assay string name of the assay to plot.
#'        (default = \code{SummarizedExperiment::assayNames(se)[1]})
#' @param stratify string name of column in colData of SummarizedExperiment to
#'        use for comparing samples; values in column must be
#'        unique for each sample. Alternatively, can specify '\code{"sample"}' to
#'        use column names. (default = "condition")
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
pbmPlotScatter <- function(se, assay = SummarizedExperiment::assayNames(se)[1],
                           stratify = "condition", baseline = NULL,
                           log_scale = TRUE, maplot = FALSE, .method = "auto") {
    stopifnot(assay %in% SummarizedExperiment::assayNames(se))

    if (! "Row" %in% names(rowData(se)) || ! "Column" %in% names(rowData(se))) {
        if ("kmer" %in% names(rowData(se))) {
            rowData(se)$Row <- rowData(se)$kmer
            rowData(se)$Column <- rowData(se)$kmer
        } else {
            stop("rowData must either contain 'Row','Column' columns or 'kmer' column.")
        }
    }

    ## check plot stratification params
    strats <- .pbmCheckStratify(se, stratify, baseline)
    coldat <- strats$coldat
    baseline <- strats$baseline

    ## filter probes
    se <- pbmFilterProbes(se) 
    
    ## extract intensities
    pdat <- SummarizedExperiment::assay(se, assay)
    pdat <- as.data.frame(pdat, optional = TRUE)
    pdat <- tibble::as_tibble(pdat)
    pdat <- dplyr::mutate(pdat,
                          Row = rowData(se)[, "Row"],
                          Column = rowData(se)[, "Column"])
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
        if (log_scale) {
            ptitle <- paste0("PBM Intensity MA Plot (log2; reference = '", baseline, "')")
        } else {
            warning("Note that MA plots are typically drawn on the log2 scale.\n",
                    "Consider re-plotting with log_scale = TRUE if the original ",
                    "values were not on the log2 scale.")
            ptitle <- paste0("PBM Intensity MA Plot (reference = '", baseline, "')")
        }
        pxaxis <- scale_x_continuous("A; mean (condition + reference) / 2")
        pyaxis <- scale_y_continuous("M; difference (condition - reference)")
        pline <- geom_hline(yintercept = 0, color = 'dodgerblue3')
    } else {
        if (log_scale) {
            ptitle <- "PBM Intensity (log2)"
            pxaxis <- scale_x_continuous(paste0("reference (", baseline, ") intensity (log)"),
                                         breaks = 2^(0:100), trans = "log2")
            pyaxis <- scale_y_continuous("intensity (log)", breaks = 2^(0:100), trans = "log2")
        } else {
            ptitle <- "PBM Intensity"
            pxaxis <- scale_x_continuous(paste0("reference (", baseline, ") intensity"))
            pyaxis <- scale_y_continuous("intensity")
        }
        pline <- geom_abline(intercept = 0, slope = 1, color = 'dodgerblue3')
    }
    
    ## manipulate to get reference condition in separate column
    pdat <- tidyr::spread(pdat, Stratify, value)
    pdat <- dplyr::rename(pdat, Baseline = !!baseline)
    pdat <- tidyr::gather(pdat, Stratify, value, -Baseline, -Row, -Column)

    ## calculate MA values if necessary and create base of plot
    if (maplot) {
        if (log_scale) {
            pdat <- dplyr::mutate(pdat,
                                  value = log2(value),
                                  Baseline = log2(Baseline))
        }
        pdat <- dplyr::mutate(pdat,
                              M = value - Baseline,
                              A = (value + Baseline) / 2)
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


.pbmCheckStratify <- function(s, strat, bl, gp = NULL) {
    
    ## check validity of stratifying colData column
    stopifnot(strat %in% names(colData(s)))
    
    strat_vals <- colData(s)[[strat]]
    if (!is.null(gp)) {
        group_vals <- colData(s)[[gp]]
    } else {
        group_vals <- rep("noGroups", length(strat_vals))
    }
    sgtab <- table(strat_vals, group_vals)
    uniq_strat <- rownames(sgtab)
    
    if (any(colMaxs(sgtab, na.rm = TRUE) > 1L)) {
        stop("Stratifying variable '", strat, "' is not unique across samples within groups.\n",
             "Specify a different column in colData.")
    }

    ## determine baseline if necessary; check validity
    if (is.null(bl)) {
        bl <- grep("ref", uniq_strat, value = TRUE, ignore.case = TRUE)
        if (length(bl) > 1) {
            warning("Too many candidate baseline states in '", strat, "' column: ",
                    paste0(bl, collapse = ", "), ".\n",
                    "Using first match: ", bl[1], ".\n",
                    "If other value should be used, specify correct baseline condition w/ 'baseline'.")
            bl <- bl[1]
        }
    } else {
        if (! bl %in% uniq_strat) {
            stop(bl, " is not a value in '", strat, "' column.\n",
                 "Specify correct baseline condition w/ 'baseline'.")
        }
    } 

    if (any(sgtab[bl, ] < 1L)) {
        stop("Baseline condition [", bl, "] is not present in all groups.\n",
             "Missing from groups: ",
             paste(colnames(sgtab)[sgtab[bl, ] < 1L], collapse = ", "))
    }
    
    ## condition must be a unique column for faceting plot
    coldat <- data.frame(colData(s), check.names = FALSE,
                         check.rows = FALSE, stringsAsFactors = FALSE)
    coldat <- tibble::rownames_to_column(coldat, "sample")
    coldat <- dplyr::rename(coldat, Stratify = I(strat))
    if (is.null(gp)) {
        coldat <- dplyr::select(coldat, sample, Stratify)
    } else {
        coldat <- dplyr::rename(coldat, Group = I(gp))
        coldat <- dplyr::select(coldat, sample, Stratify, Group)
    }
    return(list(coldat = coldat, baseline = bl))
}
