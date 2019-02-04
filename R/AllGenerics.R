#' @rdname PBMDesign
#' @export
setGeneric("PBMDesign", valueClass = "PBMDesign",
           function(x, ...) standardGeneric("PBMDesign"),
           useAsDefault = function(x, ...) {
               if (missing(x)) {
                   warning("PBMDesign should be constructed by specifying table of probes with 'x ='.")
               } else {
                   warning("PBMDesign should be constructed by specifying table of probes with 'x ='.\n",
                           "Ignoring specified 'x='.")
               }
               .PBMDesign(design = data.frame(Sequence = character(),
                                              probeID = character()), ...)
           })

#' Set design in PBMExperiment object
#'
#' @description
#' Add, remove or place \code{\link[=PBMDesign-class]{PBMDesign}} in
#' \code{\link[=PBMExperiment-class]{PBMExperiment}} object.
#' Design can be removed by setting the value to \code{NULL}.
#' 
#' @param x \code{\link[=PBMExperiment-class]{PBMExperiment}} object.
#' @param value \code{\link[=PBMDesign-class]{PBMDesign}} or \code{NULL}.
#' 
#' @return
#' modified PBMExperiment object
#'
#' @seealso \code{\link{PBMDesign}}
#' @rdname PBMDesign-setter
#' @author Patrick Kimes
#' @export
setGeneric("PBMDesign<-", 
           function(x, value) standardGeneric("PBMDesign<-"))


