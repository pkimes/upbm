#' @rdname PBMDesign
#' @export
setGeneric("PBMDesign", valueClass = "PBMDesign",
           function(object, ...) standardGeneric("PBMDesign"),
           useAsDefault = function(object, ...) {
               if (missing(object)) {
                   warning("PBMDesign should be constructed by specifying table of probes with 'object ='.")
               } else {
                   warning("PBMDesign should be constructed by specifying table of probes with 'object ='.\n",
                           "Ignoring specified 'object ='.")
               }
               .PBMDesign(design = data.frame(Sequence = character(),
                                              probeID = character()), ...)
           })

#' @title Set design in PBMExperiment object
#'
#' @description
#' Add or replace \code{\link[=PBMDesign-class]{PBMDesign}} in
#' \code{\link[=PBMExperiment-class]{PBMExperiment}} object.
#' 
#' @param object a \code{\link[=PBMExperiment-class]{PBMExperiment}} object.
#' @param value a \code{\link[=PBMDesign-class]{PBMDesign}} object.
#' 
#' @return
#' modified PBMExperiment object
#'
#' @seealso \code{\link{PBMDesign}}
#' @name PBMDesign-replace
#' @aliases PBMDesign<-
#' @export
#' @author Patrick Kimes
setGeneric("PBMDesign<-", 
           function(object, value) standardGeneric("PBMDesign<-"))


