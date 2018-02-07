#' Plot of Intensity Densities
#'
#' Densities of all probe-level intensities are plotted.
#'
#' @param se SummarizedExperiment object containing GPR
#'        intensity information.
#' @param assay_name string name of the assay to plot.
#'        (default = "gpr")
#' @param log_scale logical whether to plot the intensities
#'        on the log-scale. (default = TRUE)
#' @param .facet logical whether plot should be faceted using
#'        the default 'condition' column from the colData of the
#'        SummarizedExperiment. (default = TRUE)
#' @param .filter integer specifying level of probe filtering to
#'        perform prior to plotting. (default = 0)
#'
#' @return
#' ggplot object
#'
#' @importFrom dplyr select
#' @import ggplot2 SummarizedExperiment
#' @export
#' @author Patrick Kimes
pbmPlotDensity <- function(se, assay_name = "gpr", log_scale = TRUE,
                           .facet = TRUE, .filter = 0) {
    stopifnot(assay_name %in% assayNames(se))

    ## condition must be a unique column for faceting plot
    if (.facet && any(duplicated(colData(se)$condition))) {
        stop("'condition' column of colData must be unique.")
    }

    ## filter probes
    se <- pbmFilterProbes(se, assay_name, .filter) 

    ## extract sample metadata
    coldat <- as.data.frame(colData(se))
    coldat <- tibble::rownames_to_column(coldat, "sample")

    ## extract intensities
    pdat <- assay(se, assay_name)
    pdat <- data.frame(pdat, stringsAsFactors = FALSE)
    pdat <- tibble::as_tibble(pdat)
    pdat <- tidyr::gather(pdat, sample, value)
    pdat <- dplyr::left_join(pdat, coldat, by = "sample")

    ## check for negative values
    pdat_pnegative <- mean(pdat$value < 0, na.rm = TRUE)

    ## warn user about dropped values
    if (log_scale && pdat_pnegative > 0) {
        warning("Log-scaling will drop ", round(pdat_pnegative * 100, 2),
                "% of values.")
    }
    
    ## handle log-scale plotting if requested
    if (log_scale) {
        ptitle <- "PBM Intensity (log2)"
        pxaxis <- scale_x_continuous("intensity (log)", breaks = 2^(0:100), trans = "log10")
    } else {
        ptitle <- "PBM Intensity"
        pxaxis <- scale_x_continuous("intensity")
    }
    
    gp <- ggplot(pdat, aes(x = value)) +
        geom_density(color = 'black', fill = 'black', alpha = 1/4) +
        theme_bw() +
        pxaxis +
        ggtitle(ptitle)

    ## facet plot using default
    if (.facet) {
        gp <- gp + facet_grid(. ~ condition)
    } 
    gp
}
