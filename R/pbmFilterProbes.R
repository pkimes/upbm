#' Filter PBM Probes
#'
#' Not all probes in a PBM experiment should be used for analysis.
#' Given a SummarizedExperiment containing PBM intensities,
#' this functions provides multiple levels of filtering to remove
#' uninformative and low quality probes based on quality control
#' metrics.
#' 
#' @param se SummarizedExperiment object containing GPR
#'        intensity information or DataFrame/data.frame object
#'        containing the probe design.
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

    if (is(se, "SummarizedExperiment")) {
        pd <- rowData(se)
    } else if (is(se, "DataFrame") | is(se, "data.frame")) {
        pd <- se
    } else {
        stop("se must be a SummarizedExperiment or DataFrame/data.frame")
    }

    if (level > 0L && ! "ID" %in% names(pd)) {
        stop(paste0("Must have 'ID' column ", ifelse(is(se, "SummarizedExperiment"), "in rowData ", ""),
                    "to use filter level > 0.\n"),
             "If **all** probes should be used, set filter level to 0.\n",
             "If only a subset of probes (i.e. de Bruijn probes) should be used, ",
             "add 'ID' column.")
        level <- 0L
    }
    
    if (level > 0L) {
        se <- se[grepl("^dBr_", pd$ID), , drop = FALSE]
    }
    
    return(se)
}
