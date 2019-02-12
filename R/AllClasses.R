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
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @exportClass PBMExperiment
#' @name PBMExperiment-class
#' @author Patrick Kimes
.PBMExperiment <- setClass("PBMExperiment", contains = "SummarizedExperiment",
                           slots = c(probeFilter = "list", probeTrim = "numeric",
                                     probeCols = "character"))

setValidity2("PBMExperiment",
             function(object) {
                 ## check rowData necessary columns (Sequence, probeID)
                 rd <- rowData(object)
                 if (!all(c("Sequence", "probeID") %in% colnames(rd))) {
                     stop("PBMExperiment must include probe sequence information ",
                          "as rowData columns: 'Sequence' and 'probeID'.")
                 }
                 if (!is.vector(rd$probeID, mode = "character")) {
                     stop("PBMExperiment probe IDs must be character strings.\n",
                          "Please check 'probeID' rowData values.")
                 }
                 if (!is.vector(rd$Sequence, mode = "character")) {
                     stop("PBMExperiment probe sequences must be character strings.\n",
                          "Please check 'Sequence' rowData values.")
                 }
                 ## check probeFilter specification
                 if (length(object@probeFilter) > 1) {
                     if (is.null(names(object@probeFilter)) | any(names(object@probeFilter) == "")) {
                         stop("PBMExperiment 'probeFilter' must be a named list.")
                     }
                     if (!all(names(object@probeFilter) %in% intersect(colnames(rd), object@probeCols))) {
                         stop("PBMExperiment 'probeFilter' must be a named list matching ",
                              "rowData columns also in 'probeCols'.")
                     }
                 }
                 ## check probeTrim specification
                 if (length(object@probeTrim) != 0L && length(object@probeTrim) != 2L) {
                     stop("PBMExperiment 'probeTrim' specification is invalid. \n",
                          "If specified, 'probeTrim' must be a vector of length 2 corresponding ",
                          "to the start and end of the probe sequence to keep for analysis.")
                 }
                 if (length(object@probeTrim) == 2L &&
                     (object@probeTrim[1] > object@probeTrim[2] || object@probeTrim[1] < 0)) {
                     stop("PBMExperiment 'probeTrim' specification is invalid. \n",
                          "If specified, 'probeTrim' must be a vector of length 2 with ",
                          "0 <= probeTrim[1] <= probeTrim[2].")
                 }
                 ## check probeCols specification
                 if (!all(c("Sequence", "probeID") %in% object@probeCols)) {
                     stop("PBMExperiment 'probeCols' must include names of columns in rowData ",
                          "corresponding to probe design information.\n",
                          "At a minimum, this should include: 'Sequence' and 'probeID'.")
                 }
                 if (!all(object@probeCols %in% colnames(rd))) {
                     stop("PBMExperiment 'probeCols' should only include names of columns in rowData ",
                          "corresponding to probe design information.\n",
                          "The following columns listed in 'probeCols' are not found in rowData: \n",
                          paste0(setdiff(object@probeCols, colnames(rd)), collapse = ", "), ".")
                 }
                 ## check properties of filtered object
                 objf <- pbmFilterProbes(object)
                 ## check filtered probeIDs are unique
                  if (any(duplicated(rowData(objf)$probeID))) {
                     stop("PBMExperiment probe IDs must be unique after filtering.\n",
                          "Please check 'probeID' rowData values.")
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
.PBMDesign <- setClass("PBMDesign",
                       slots = c(design = "data.frame", probeFilter = "list",
                                 probeTrim = "numeric"))

setValidity2("PBMDesign",
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
                     stop("PBMDesign 'probeTrim' specification is invalid. \n",
                          "If specified, 'probeTrim' must be a vector of length 2 corresponding ",
                          "to the start and end of the probe sequence to keep for analysis.")
                 }
                 if (length(object@probeTrim) == 2L &&
                     (object@probeTrim[1] > object@probeTrim[2] || object@probeTrim[1] < 0)) {
                     stop("PBMDesign 'probeTrim' specification is invalid. \n",
                          "If specified, 'probeTrim' must be a vector of length 2 with ",
                          "0 <= probeTrim[1] <= probeTrim[2].")
                 }
                 ## check properties of filtered object
                 objf <- pbmFilterProbes(object)
                 ## check filtered probeIDs are unique
                  if (any(duplicated(objf@design$probeID))) {
                     stop("PBMDesign probe IDs must be unique after filtering.\n",
                          "Please check 'probeID' column values of design.")
                 }
                 TRUE
             })
