#' PBMExperiment class
#'
#' @description
#' Extension of the \code{\link{SummarizedExperiment}} class to store
#' universal protein binding microarray (uPBM) data for use with functions
#' in the \pkg{upbm} package. 
#'
#' @slot probeFilter an optional named list of probe filters to be used to subset
#'       probes during data analysis steps. List names must correspond to columns in
#'       RowData. List entries must be single-parameter functions to be called on the
#'       corresponding column to return a logical vector of probes to keep (TRUE) and
#'       drop (FALSE) during analysis.
#' @slot probeTrim an optional integer vector of length 2 specifying start and end
#'       positions in probe `Sequence' to use in analysis steps.
#' @slot probeCols an optional character vector of rowData column names corresponding
#'       to probe design information. 
#' 
#' @aliases PBMExperiment-class
#' @import SummarizedExperiment
#' @export
#' @exportClass PBMExperiment
#' @name PBMExperiment-class
#' @author Patrick Kimes
setClass("PBMExperiment", contains = "SummarizedExperiment",
         slots = c(probeFilter = "list", probeTrim = "numeric",
                   probeCols = "character"))

setValidity("PBMExperiment",
            function(object) {
                rd <- rowData(object)
                if (!all(c("Sequence", "ID") %in% colnames(rd))) {
                    stop("PBMExperiment must include probe sequence information ",
                         "as rowData columns: 'Sequence' and 'ID'.")
                }
                if (!is.vector(rd$ID, mode = "character")) {
                    stop("PBMExperiment probe IDs must be character strings.\n",
                         "Please check 'ID' rowData values.")
                }
                if (!is.vector(rd$Sequence, mode = "character")) {
                    stop("PBMExperiment probe sequences must be character strings.\n",
                         "Please check 'Sequence' rowData values.")
                }
                ## if (any(nchar(rd$Sequence) != nchar(rd$Sequence[1]))) {
                ##     stop("PBMExperiment probe sequences must all be same length.\n",
                ##          "Please check 'Sequence' rowData values.")
                ## }
                if (any(duplicated(rd$ID))) {
                    stop("PBMExperiment probe IDs must be unique.\n",
                         "Please check 'ID' rowData values.")
                }
                ## check probeFilter specification
                if (length(object@probeFilter) > 1) {
                    if (is.null(names(object@probeFilter)) | any(names(object@probeFilter) == "")) {
                        stop("PBMExperiment 'probeFilter' must be a named list.")
                    }
                    if (!all(names(object@probeFilter) %in% colnames(rd))) {
                        stop("PBMExperiment 'probeFilter' must be a named list matching ",
                             "rowData columns.")
                    }
                }
                ## check probeTrim specification
                if (length(object@probeTrim) != 0L && length(object@probeTrim) != 2L) {
                    stop("PBMExperiment 'probeTrim' is too long. \n",
                         "If specified, 'probeTrim' must be a vector of length 2 corresponding ",
                         "to the start and end of the probe sequence to keep for analysis.")
                }
                TRUE
            })


#' PBMDesign class
#'
#' @description
#' Simple class to store array design information for universal protein binding
#' microarray (uPBM) data for use with functions in the \pkg{upbm} package.
#' Array information is comprised of three elements.
#'
#' @slot design a data.frame with each row corresponding to a probe on the array.
#'       Must include `Sequence' and (unique) `probeID' columns, along with any
#'       other metadata for probes, e.g. array `Row' or `Column' spatial coordinates. 
#' @slot probeFilter an optional named list of probe filters to be used to subset
#'       probes during data analysis steps. List names must correspond to columns in `design'.
#'       List entries must be single-parameter functions to be called on the
#'       corresponding column to return a logical vector of probes to keep (TRUE) and
#'       drop (FALSE) during analysis.
#' @slot probeTrim an optional integer vector of length 2 specifying start and end
#'       positions in probe `Sequence' to use in analysis steps.
#' 
#' @aliases PBMDesign-class
#' @export
#' @exportClass PBMDesign
#' @name PBMDesign-class
#' @author Patrick Kimes
setClass("PBMDesign",
         slots = c(design = "data.frame", probeFilter = "list", probeTrim = "numeric"))

setValidity("PBMDesign",
            function(object) {
                ## check necessary columns are specified
                cn <- colnames(object@design)
                if (!all(c("Sequence", "probeID") %in% cn)) {
                    stop("PBMDesign 'design' data.frame must contain 'Sequence' and 'probeID' columns.")
                }
                ## check probeFilter specification
                if (length(object@probeFilter) > 1) {
                    if (is.null(names(object@probeFilter)) | any(names(object@probeFilter) == "")) {
                        stop("PBMDesign 'probeFilter' must be a named list.")
                    }
                    if (!all(names(object@probeFilter) %in% colnames(rd))) {
                        stop("PBMDesign 'probeFilter' must be a named list matching ",
                             "'design' data.frame columns.")
                    }
                }
                ## check probeTrim specification
                if (length(object@probeTrim) != 0L && length(object@probeTrim) != 2L) {
                    stop("PBMDesign 'probeTrim' is too long. \n",
                         "If specified, 'probeTrim' must be a vector of length 2 corresponding ",
                         "to the start and end of the probe sequence to keep for analysis.")
                }
                TRUE
            })
