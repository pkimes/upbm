#' Filter PBM Probes
#'
#' Not all probes in a PBM experiment should be used for analysis.
#' Given a SummarizedExperiment containing PBM intensities,
#' this functions provides multiple levels of filtering to remove
#' uninformative and low quality probes based on quality control
#' metrics.
#' 
#' @param se SummarizedExperiment object containing GPR
#'        intensity information.
#' @param assay_name string name of the assay to plot.
#'        (default = "gpr")
#' @param level integer specifying level of probe filtering to
#'        perform prior to plotting. (default = 0)
#'
#' @return
#' SummarizedExperiment object
#'
#' @details
#' The function supports the following levels of filtering.
#' All filtering is progressive, i.e. setting level = 2 will
#' also include filtering at level = 1.
#' * **0**: no filtering.
#' * **1**: filter on de Bruijn probes.
#'
#' @md
#' @export
#' @author Patrick Kimes
pbmFilterProbes <- function(se, assay_name, level = 0L) {

    if (level > 0L && ! "ID" %in% names(rowData(se))) {
        warning("Must have 'ID' column in rowData to use filter level > 0")
        level <- 0L
    }
    
    if (level > 0L) {
        se <- se[grepl("^dBr_", rowData(se)$ID), ]
    }
    
    return(se)
}
