#' PBMExperiment class
#'
#' This class is an extension of the \code{\link{SummarizedExperiment}} class to
#' store protein binding microarray data.
#' 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import S4Vectors SummarizedExperiment
#' @exportClass PBMExperiment
#' @export
#' @author Patrick Kimes
setClass("PBMExperiment",
         contains = "SummarizedExperiment")


#' Constructor function for PBMExperiment objects.
#'
#' Function to construct \code{PBMExperiment} objects.
#'
#' @param assays a list containing.
#' @param colData a \code{\link{DataFrame}}.
#' @param rowData a \code{\link{DataFrame}}.
#' @param ... additional parameters passed to \code{\link{SummarizedExperiment}}.
#'
#' @return
#' a \code{\link{PBMExperiment}} object.
#'
#' @aliases PBMExperiment
#' @export
#' @author Patrick Kimes
PBMExperiment <- function(assays, colData = NULL, rowData = NULL, ...) {
    se <- SummarizedExperiment(assays, rowData = rowData, colData = colData, ... )
    new("PBMExperiment", se)
}


#' Read Standard GPR File
#'
#' Simple wrapper for reading standard GPR files.
#' GPR files are tab-delimited with 35 header lines.
#' 
#' @param f path to gpr file
#'
#' @return
#' a tibble.
#'
#' @keywords internal
#' @importFrom readr read_tsv
#' @author Patrick Kimes
read_gpr <- function(f) {
    read_tsv(f, skip = 35)
}
