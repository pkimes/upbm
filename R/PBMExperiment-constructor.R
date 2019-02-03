#' @title Load raw PBM scan data as PBMExperiment
#'
#' @description
#' Read PBM scan data from GenePix Results (GPR) files and combine into
#' a single \code{\link[PBMExperiment-class]{PBMExperiment}} object.
#' Paths to GPR file must be specified as part of a larger 
#' 
#' @param x a data.frame of scans with at least a single `gpr' column
#'        corresponding to GPR file paths.
#' @param useMedian a logical value whether to use median fluorescent intensity
#'        for each probe rather than mean fluorescent intensity.
#'        (default = TRUE)
#' @param filterFlags a logical value whether to replace intensity values at probes
#'        flagged manually or automatically as being low quality
#'         (`Bad': -100, `Absent': -75, `Not Found': -50) with NA. (default = TRUE)
#' @param readBackground a logical value whether to also read in probe background
#'        intensities. (default = TRUE)
#' @param probes a \code{\link[PBMDesign-class]{PBMDesign}} object containing
#'        the probe design used for the data. (default = NULL)
#' @param ... optionally, probe design information can be directly specified using
#'        \code{design=}, \code{probeFilter=}, \code{probeTrim=} parameters. These
#'        will be ignored if \code{probes=} is specified.
#'
#' @return
#' \code{\link[PBMExperiment-class]{PBMExperiment}} object with assays containing
#' intensities from the GPR files specified in the table of scans. The object
#' includes the following assays:
#' * \code{fore} - foreground probe intensities,
#' * \code{back} - background probe intensities (unless \code{readBackground = FALSE}).
#'
#' Scan metadata included in the data.frame are stored as column data.
#' Probe design information, if specified, is included as row data.
#' 
#' @details
#' GPR files are parsed differently depending on the file type which can take values
#' of ``Alexa", ``Cy3", or ``Masliner". Both ``Alexa" and ``Cy3" scans
#' are treated as untouched GenePix Results (GPR) files, while ``Masliner" files are
#' treated as Masliner-processed GPR files.
#' 
#' If a `scan' column is not provided in the input table, GPR files are assumed to
#' be Alexa488 scans unless the file name in the `gpr' column contains "Masliner".
#' 
#' If ``raw data" files (processed GPR data files, e.g. uploaded to UniPROBE)
#' include probe ID, Name, or Sequence information, this is compared across samples and
#' an error is thrown if the probe-to-sequence mapping is not common across all
#' samples.
#'
#' @seealso \code{\link{PBMExperiment-class}}
#' @importFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#' @importFrom purrr reduce
#' @importFrom dplyr select left_join
#' @importFrom tibble as_tibble
#' @name PBMExperiment
#' @export
#' @author Patrick Kimes
PBMExperiment <- function(x, useMedian = TRUE,
                          filterFlags = TRUE, readBackground = TRUE,
                          probes = NULL, ...) {

    ## check validity of inputs
    stopifnot(is.data.frame(tab))
    stopifnot("gpr" %in% names(tab))
    stopifnot(is.logical(useMedian))
    stopifnot(is.logical(filterFlags))

    ## guess scan type if not specified - only care if Masliner or not
    if ("scan" %in% names(tab)) {
        tab_scan <- tab$scan
        if (!all(tab_scan %in% c("Alexa", "Cy3", "Masliner"))) {
            stop("\"scan\" column entries can only take the following values: \n",
                 "\"Alexa\", \"Cy3\", \"Masliner\".")
        }
    } else {
        ## don't add to 'tab' since don't want in colData of output
        tab_scan <- rep("Alexa", nrow(tab))
        tab_scan[grep("MaslinerOutput", tab$gpr,
                      ignore.case = TRUE)] <- "Masliner"
    }
    
    ## parse probe parameters
    stopifnot(is.null(probes) || is(probes, "PBMDesign"))
    if (is.null(probes)) {
        pparam <- list(...)
        if (length(pparam) > 0L) {
            probes <- PBMDesign(...)
        }
    }
    
    ## read in all data files
    assay_table <- mapply(readGPR, gpr_path = tab$gpr, gpr_type = tab_scan,
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
    coldat <- S4Vectors::DataFrame(dplyr::select(tab, -gpr))
    rownames(coldat) <- names(assaydat[["fore"]])
    
    new("PBMExperiment", assays = assaydat,
        rowData = rowdat, colData = coldat,
        metadata = list(steps = list(),
                        spatialAdjustment = NULL,
                        backgroundCorrection = NULL,
                        betweenArrayNormalization = NULL),
        probeFilter = ifelse(is.null(probes), list(), probes@probeFilter),
        probeTrim = ifelse(is.null(probes), numeric(), probes@probeTrim),
        probeCols = ifelse(is.null(probes), character(), names(probeCols)))
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
        stop("Foreground column could not be identified!\n",
             "Try setting 'scan' column in sample table to 'RawData'.")
    }
    colt[c(fore_idx)] <- 'd'
    if (readBackground) {
        if (length(back_idx) == 0) {
            stop("Background column could not be identified! \n",
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


