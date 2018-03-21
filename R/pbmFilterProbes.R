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
#' @param level integer specifying level of probe filtering to
#'        perform prior to plotting. (default = 0)
#'
#' @return
#' SummarizedExperiment object.
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
pbmFilterProbes <- function(se, level = 0L) {

    if (level > 0L && ! "ID" %in% names(rowData(se))) {
        stop("Must have 'ID' column in rowData to use filter level > 0.\n",
             "If **all** probes should be used, set .filter = 0.\n",
             "If only a subset of probes (i.e. de Bruijn probes) should be used,",
             "add 'ID' to rowData.")
        level <- 0L
    }
    
    if (level > 0L) {
        se <- se[grepl("^dBr_", rowData(se)$ID), ]
    }
    
    return(se)
}
