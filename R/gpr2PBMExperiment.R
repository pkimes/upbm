#' @title Load raw PBM scan data as PBMExperiment
#'
#' @description
#' This function reads PBM scan data from GenePix Results (GPR) files and returns
#' a single \code{\link[=PBMExperiment-class]{PBMExperiment}} object containing
#' foreground and optionally background probe intensity data in assays. Each
#' GPR is stored as a column in the object, with metadata specified in the
#' input \code{scans} table stored in the colData of the object.
#'
#' If probe design information is also specified as a \code{\link[=PBMDesign-class]{PBMDesign}}
#' object to the optional \code{probes=} argument, this information is also stored
#' in the rowData of the returned PBMExperiment.
#'
#' @param scans a data.frame of GPR scan information with at least a single `gpr'
#'        column, and preferrably a second `type' column. See details for more
#'        information on the content of these columns.
#' @param useMedian a logical value whether to use median fluorescent intensity
#'        for each probe rather than mean fluorescent intensity.
#'        (default = TRUE)
#' @param filterFlags a logical value whether to replace intensity values at probes
#'        flagged manually or automatically as being low quality
#'         (`Bad': -100, `Absent': -75, `Not Found': -50) with NA. (default = TRUE)
#' @param readBackground a logical value whether to also read in probe background
#'        intensities. (default = TRUE)
#' @param probes a \code{\link[=PBMDesign-class]{PBMDesign}} object containing
#'        the probe design used for the data. (default = NULL)
#' @param ... optionally, probe design information can be directly specified using
#'        \code{design=}, \code{probeFilter=}, \code{probeTrim=} parameters. These
#'        will be ignored if \code{probes=} is specified.
#'
#' @return
#' \code{\link[=PBMExperiment-class]{PBMExperiment}} object with assays containing
#' intensities from the GPR files specified in the table of scans. The object
#' includes the following assays:
#' * \code{fore} - foreground probe intensities,
#' * \code{back} - background probe intensities (unless \code{readBackground = FALSE}).
#'
#' Scan metadata included in the data.frame are stored as colData.
#' Probe design information, if specified, is included as rowData.
#' 
#' @details
#' The primary argument, \code{scans}, must be a data.frame of GPR scan metadata.
#' Each row of the data.frame should correspond to a single GPR scan. At a minimum
#' the data.frame must include a single \code{"gpr"} column specifying the paths
#' to the corresponding GPR files to be read in. 
#'
#' Additionally, GPR files are parsed differently depending on the file type which must be
#' one of ``Alexa", ``Cy3", or ``Masliner". Both ``Alexa" and ``Cy3" scans
#' are treated as untouched GPR files, while ``Masliner" files treated as Masliner-processed
#' GPR files. If possible the scan type of each file should be
#' specified in a \code{"type"} column of the \code{scans} data.frame.
#' If a \code{"type"} column is not included, all GPR files are assumed to be Alexa488
#' scans unless the file name in the \code{"gpr"} column contains ``Masliner".
#'
#' @seealso \code{\link{PBMExperiment-class}}, \code{\link{PBMExperiment}}, \code{\link{PBMDesign-class}}
#' @importFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#' @importFrom purrr reduce
#' @importFrom dplyr select left_join
#' @importFrom tibble as_tibble
#' @name gpr2PBMExperiment
#' @export
#' @author Patrick Kimes
gpr2PBMExperiment <- function(scans, useMedian = TRUE,
                              filterFlags = TRUE, readBackground = TRUE,
                              probes = NULL, ...) {

    ## check validity of inputs
    stopifnot(is.data.frame(scans))
    stopifnot("gpr" %in% names(scans))
    stopifnot(is.logical(useMedian))
    stopifnot(is.logical(filterFlags))

    ## guess scan type if not specified - only care if Masliner or not
    if ("type" %in% names(scans)) {
        tab_type <- scans$type
        if (!all(tab_type %in% c("Alexa", "Cy3", "Masliner"))) {
            stop("\"type\" column entries can only take the following values: \n",
                 "\"Alexa\", \"Cy3\", \"Masliner\".")
        }
    } else {
        ## don't add to 'types' since don't want in colData of output
        tab_type <- rep("Alexa", nrow(scans))
        tab_type[grep("MaslinerOutput", scans$gpr,
                      ignore.case = TRUE)] <- "Masliner"
    }
    
    ## parse probe parameters
    stopifnot(is.null(probes) || is(probes, "PBMDesign"))
    if (is.null(probes)) {
        pparam <- list(...)
        if (length(pparam) > 0L) {
            probes <- new("PBMDesign", ...)
        }
    }
    
    ## read in all data files
    assay_table <- mapply(readGPR, gpr_path = scans$gpr, gpr_type = tab_type,
                          useMedian = useMedian, filterFlags = filterFlags,
                          readBackground = readBackground, SIMPLIFY = FALSE)

    ## merge GPR data across samples
    if (readBackground) {
        assay_btable <- lapply(assay_table, `[`, c("Column", "Row", "background"))
        assay_btable <- purrr::reduce(assay_btable, left_join, by = c("Column", "Row"))
    }
    assay_table <- lapply(assay_table, `[`, c("Column", "Row", "foreground"))
    assay_table <- purrr::reduce(assay_table, left_join, by = c("Column", "Row"))

    ## convert GPR data to list of DataFrames for SummarizedExperiment assay slot
    assaydat <- list(fore = S4Vectors::DataFrame(dplyr::select(assay_table, -Column, -Row)))
    names(assaydat[["fore"]]) <- paste0("s", 1:ncol(assaydat[["fore"]]))
    if (readBackground) {
        assaydat[["back"]] <- S4Vectors::DataFrame(dplyr::select(assay_btable, -Column, -Row))
        names(assaydat[["back"]]) <- names(assaydat[["fore"]])
    }
    
    ## row/probe-level metadata from GPR files
    rowdat <- dplyr::select(assay_table, Column, Row)

    ## check if probe sequence dimension matches data
    if (!is.null(probes) && nrow(probes@design) != nrow(rowdat)) {
        warning("Dimension of specified 'probes' does not match GPR files.\n",
                "Ignoring specified 'probes' input.")
        probes <- NULL
    }

    ## add probes to metadata if provided
    if (!is.null(probes)) {
        pdesign <- probes@design
        ovnames <- intersect(names(rowdat), names(pdesign))
        if (length(ovnames) == 0) {
            warning("Specified 'probes' does not have columns matching GPR row metadata.\n",
                    "The probes will be added to the table in the order provided.\n",
                    "WARNING: This may lead to incorrect joining.\n",
                    "  Consider adding probe 'Row' and 'Column' information to the\n",
                    "  probe design table for safer merging.")
            rowdat <- cbind(rowdat, pdesign)
        } else {
            if (!all(c("Column", "Row") %in% ovnames)) {
                warning("Specified 'probes' does not have 'Column' and 'Row' columns.\n",
                        "The probes will be matched using other columns.\n",
                        "WARNING: This may lead to incorrect joining.",
                        "  Consider adding probe 'Row' and 'Column' information to the\n",
                        "  probe design table for safer merging.")
            }
            rowdat <- dplyr::left_join(rowdat, tibble::as_tibble(as.data.frame(pdesign)),
                                       by = ovnames)
        }
    }

    ## column/condition-level metadata
    coldat <- S4Vectors::DataFrame(dplyr::select(scans, -gpr))
    rownames(coldat) <- names(assaydat[["fore"]])

    if (is.null(probes)) {
        pe <- PBMExperiment(assays = assaydat, rowData = rowdat, colData = coldat)
    } else {
        pe <- PBMExperiment(assays = assaydat, rowData = rowdat, colData = coldat,
                            probeFilter = probes@probeFilter,
                            probeTrim = probes@probeTrim,
                            probeCols = names(probes@design))
    }
    pe
}


#' Read GPR File
#' 
#' Helper function for reading in a single GPR file.
#'
#' @param gpr_path path to GPR file.
#' @param gpr_type scan type; one of "Alexa", "Cy3", "Masliner".
#' @param useMedian logical whether to use median fluorescent intensity
#'        for each probe rather than mean fluorescent intensity.
#'        (default = TRUE)
#' @param filterFlags logical whether to replace intensity values at probes
#'        flagged manually or automatically as being low quality.
#'         ('Bad': -100, 'Absent': -75, 'Not Found': -50) with NA. (default = TRUE)
#' @param readBackground logical whether to also read in probe background
#'        intensities. (default = TRUE)
#'
#' @return
#' tibble (data.frame-like) object of a single GPR file with three or four
#' columns: 'Column', 'Row', 'foreground', (optionally) 'background'.
#' ('ID' and 'Name' columns are ignored as these may be incorrect in the GPR file.) 
#' 
#' @details
#' Since the names of foreground and background columns can vary across samples,
#' they are identified using regular expressions. Therefore, if the column names
#' deviate from the expected format ("F.* Mean", "F.* Median". "B.* Mean", "B.* Median"),
#' the function may fail to read the foreground and background intensities.
#'
#' @keywords internal
#' @importFrom readr read_tsv read_lines
#' @importFrom dplyr select
#' @author Patrick Kimes
readGPR <- function(gpr_path, gpr_type, useMedian = TRUE, filterFlags = TRUE,
                    readBackground = TRUE) {

    if (! gpr_type %in% c("Alexa", "Cy3", "Masliner")) {
        stop("gpr_type must be one of: \"Alexa\", \"Cy3\", \"Masliner\".")
    }

    ## number of rows to skip appears to be variable - determine from reading raw
    header <- readr::read_lines(gpr_path, n_max = 50)
    header_ln <- grep(".*Column.*Row.*Name.*", header)
    header <- strsplit(header[header_ln], "\\t")[[1]]
    if (length(header_ln) != 1) {
        stop("After checking first 50 lines, cannot determine header row ",
             "in file: \n", basename(gpr_path))
    }

    ## determine number of columns
    p <- length(header)

    ## specify which columns to read in - Column, Row index
    colt <- rep("-", p)
    colt[c(2, 3)] <- 'i'
    if (filterFlags) {
        colt[grep("\"Flags\"", header)] <- 'i'
    }
    
    ## column corresponding to intensities
    if (useMedian) {
        fore_idx <- grep("^\"F.* Median\"$", header)
        back_idx <- grep("^\"B.* Median\"$", header)
    } else {
        fore_idx <- grep("^\"F.* Mean\"$", header)
        back_idx <- grep("^\"B.* Mean\"$", header)
    }
    if (length(fore_idx) == 0) {
        stop("Foreground column could not be identified for the following file!\n",
             gpr_path, "\n",
             "Please verify that this is a valid GPR or Masliner GPR file.")
    }
    colt[c(fore_idx)] <- 'd'
    if (readBackground) {
        if (length(back_idx) == 0) {
            stop("Background column could not be identified for the following file!\n",
                 gpr_path, "\n",
                 "Please verify that this is a valid GPR or Masliner GPR file.\n",
                 "If background intensities are not needed, try rerunning ",
                 "with readBackground = FALSE.")
        }
        colt[c(back_idx)] <- 'd'
    }

    colts <- paste(colt, collapse = "")
    vals <- readr::read_tsv(gpr_path, skip = header_ln - 1, col_types = colts,
                            progress = FALSE)
    names(vals)[which(which(colt != "-") == fore_idx)] <- 'foreground'
    if (readBackground) {
        names(vals)[which(which(colt != "-") == back_idx)] <- 'background'
    }

    ## remove negative (low quality) flagged probes
    if (filterFlags) {
        vals$foreground[vals$Flags < 0] <- NA_real_
        if (readBackground) {
            vals$background[vals$Flags < 0] <- NA_real_
        }
        vals <- dplyr::select(vals, -Flags)
    }
    vals
}


