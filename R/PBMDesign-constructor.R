#' Create a new PBMDesign object
#'
#' @description
#' Create a new PBMDesign object of protein binding microarray probe design information.
#' Alternatively, the function can be called on a PBMExperiment to extract the
#' probe design information associated with experimental data.
#'
#' @param object a data.frame with each row corresponding to a probe on the array.
#'        Must include `Sequence' and (unique) `probeID' columns, along with any
#'        other metadata for probes, e.g. array `Row' or `Column' spatial coordinates.
#'        Alternatively, a \code{\link[=PBMExperiment-class]{PBMExperiment}} object to
#'        return the associated \code{PBMDesign} object.
#' @param ... optional probe design parameters to be defined as part of the \code{PBMDesign}
#'        object. See the \code{\link[=PBMDesign-class]{PBMDesign}} class definition for
#'        a list of probe design parameters. Important parameters are as described in
#'        the Details section below.
#'
#' @details
#' Probe design parameters can be specified by name. The following are a couple important
#' parameters which are defined by default for universal PBM designs in the \pkg{upbmAux}
#' package.
#' \enumerate{
#' \item \code{probeFilter}: optional named list of probe filters to be used to subset
#' probes during data analysis steps. List names must correspond to columns in `design'.
#' List entries must be single-parameter functions to be called on the corresponding column
#' to return a logical vector of probes to keep (TRUE) and drop (FALSE) during analysis.
#' \item \code{probeTrim}: optional integer vector of length 2 specifying start and end
#' positions in probe `Sequence' to use in analysis steps.
#' }
#'
#' @return
#' \code{PBMDesign} object.
#'
#' @references
#' \itemize{
#' \item Berger, M. F., & Bulyk, M. L. (2017). Universal Protein Binding Microarray (PBM) Analysis Suite Homepage. Retrieved October 16, 2020, from \url{http://thebrain.bwh.harvard.edu/PBMAnalysisSuite/indexSep2017.html}
#' }
#
#' @examples
#' ## Universal array designs included with the Universal PBM Analysis Suite software
#' ## available at the referenced link can be read in as data frames or tibbles (here
#' ## as an object 'mydesign') and converted to a PBMDesign object.
#' ## The 'probeFilter=' and 'probeTrim=' settings here filter to de Bruijn sequences
#' ## and use only the first 36 bases of each probe sequence for downstream analysis.
#' \dontrun{
#' PBMDesign(
#'     object = mydesign, 
#'     probeFilter = list(probeID = function(x) { grepl("^dBr", x) }), 
#'     probeTrim = c(1, 36)
#' )
#' }
#' 
#' @seealso \code{\link{PBMDesign-class}}
#' @name PBMDesign
#' @export
#' @author Patrick Kimes
NULL

.PBMDesign.table <- function(object, ...) {
    new("PBMDesign", design = object, ...)
}

.PBMDesign.PBMExperiment <- function(object) {
    new("PBMDesign", design = as.data.frame(rowData(object)[, object@probeCols, drop = FALSE], optional = TRUE),
        probeFilter = object@probeFilter, probeTrim = object@probeTrim)
}

#' @rdname PBMDesign
#' @exportMethod "PBMDesign"
setMethod("PBMDesign", signature(object = "data.frame"), .PBMDesign.table)

#' @rdname PBMDesign
#' @exportMethod "PBMDesign"
setMethod("PBMDesign", signature(object = "DataFrame"), .PBMDesign.table)

#' @rdname PBMDesign
#' @exportMethod "PBMDesign"
setMethod("PBMDesign", signature(object = "PBMExperiment"), .PBMDesign.PBMExperiment)
