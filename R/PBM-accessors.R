#' @title PBM slot accessors and setters
#'
#' @description
#' Simple functions for accessing and setting specialized slots
#' for PBM associated classes (PBMDesign, PBMExperiment).
#' 
#' @param object a PBMDDesign or PBMExperiment object.
#' @param ... further arguments, perhaps used by methods.
#' @param value a new value to be assigned to the corresponding components
#'        of object.
#'
#' @return
#' slot value for accessors, modified object for setters.
#'
#' @name PBMclass-accessors
#' @export
#' @author Patrick Kimes
probeFilter <- function(object) {
    stopifnot(is(object, "PBMDesign") || is(object, "PBMExperiment"))
    object@probeFilter
}

#' @rdname PBMclass-accessors
#' @export
probeTrim <- function(object) {
    stopifnot(is(object, "PBMDesign") || is(object, "PBMExperiment"))
    object@probeTrim
}

#' @rdname PBMclass-accessors
#' @export
probeCols <- function(object) {
    stopifnot(is(object, "PBMDesign") || is(object, "PBMExperiment"))
    if (is(object, "PBMExperiment")) {
        return(object@probeCols)
    } else {
        return(colnames(object@design))
    }
}

#' @rdname PBMclass-accessors
#' @export
`probeFilter<-` <- function(object, value) {
    stopifnot(is(object, "PBMDesign") || is(object, "PBMExperiment"))
    object@probeFilter <- value
    validObject(object)
    object
}

#' @rdname PBMclass-accessors
#' @export
`probeTrim<-` <- function(object, value) {
    stopifnot(is(object, "PBMDesign") || is(object, "PBMExperiment"))
    object@probeTrim <- value
    validObject(object)
    object
}

#' @rdname PBMclass-accessors
#' @importMethodsFrom BiocGenerics "design<-"
#' @exportMethod "design<-"
setReplaceMethod("design",
                 signature(object = "PBMDesign", value = "data.frame"),
                 function (object, ..., value) {
                     object@design <- value
                     object
                 })

#' @rdname PBMclass-accessors
#' @importMethodsFrom BiocGenerics design 
#' @exportMethod "design"
setMethod("design", signature(object = "PBMDesign"),
          function(object) {
              object@design
          })

#' @rdname PBMDesign-replace
#' @exportMethod "PBMDesign<-"
setReplaceMethod("PBMDesign",
                 signature(object = "PBMExperiment", value = "PBMDesign"),
                 function(object, value) {
                     se <- as(object, "SummarizedExperiment")
                     rowData(se)$probeID <- NULL
                     rowData(se)$Sequence <- NULL
                     PBMExperiment(se, pbmDesign = value)
                 })
