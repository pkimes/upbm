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
#' This function will soon support multiple levels of filter.
#' Currently, the function is simply a placeholder.
#' 
#' @export
#' @author Patrick Kimes
pbmFilterProbes <- function(se, assay_name, level = 0) {

    ## placeholder

    return(se)
}
