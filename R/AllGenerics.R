#' @rdname PBMDesign
#' @export
setGeneric("PBMDesign", valueClass = "PBMDesign",
           function(x, ...) standardGeneric("PBMDesign"))

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


