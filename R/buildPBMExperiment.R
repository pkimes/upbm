#' Read PBMExperiment Data
#'
#' Read PBM data from a table containing paths to GPR files. 
#' 
#' @param tab table of samples with at least a single `gpr` column
#'        corresponding to GPR file paths.
#' @param useMean logical whether to use mean fluorescent intensity
#'        for each probe rather than median fluorescent intensity.
#'        (default = FALSE)
#' @param filterFlags logical whether to replace intensity values at probes
#'        flagged manually or automatically as being low quality
#'         ('Bad': -100, 'Absent': -75, 'Not Found': -50) with NA. (default = TRUE)
#' @param readBackground logical whether to also read in probe background
#'        intensities. (default = TRUE)
#' @param probes data.frame containing probe sequences in a column, 'Sequence',
#'        or a character vector specifying the probe sequences. If specified,
#'        these values will be added to the rowData of the returned
#'        SummarizedExperiment object. (default = NULL)
#'
#' @return
#' SummarizedExperiment object with assays containing intensities for the
#' samples specified in the input 'tab'. The object currently includes the 
#' following assays:
#' * 'fore' - foreground probe intensities,
#' * 'back' - background probe intensities (unless `readBackground = FALSE`).
#'
#' @details
#' If raw data files (cleaned GPR data files, e.g. uploaded to UniPROBE) include
#' probe ID, Name, or Sequence information, this is compared across samples and
#' an error is thrown if the probe-to-sequence mapping is not common across all
#' samples. 
#' 
#' @import SummarizedExperiment
#' @importFrom purrr reduce
#' @importFrom dplyr select left_join
#' @importFrom tibble as_tibble
#' @export
#' @author Patrick Kimes
buildPBMExperiment <- function(tab, useMean = FALSE, filterFlags = TRUE,
                               readBackground = TRUE, probes = NULL) {
    ## check validity of inputs
    stopifnot(is.data.frame(tab))
    stopifnot("gpr" %in% names(tab))
    stopifnot(is.logical(useMean))
    stopifnot(is.logical(filterFlags))
    
    ## guess scan type if not specified - only care if Masliner or not
    if ("scan" %in% names(tab)) {
        tab_scan <- tab$scan
        if (!all(tab_scan %in% c("Alexa", "Cy3", "Masliner", "RawData"))) {
            stop("\"scan\" column entries can only take the following values: \n",
                 "\"Alexa\", \"Cy3\", \"Masliner\", \"RawData\".")
        }
    } else {
        ## don't add to 'tab' since don't want in colData of output
        tab_scan <- rep("Alexa", nrow(tab))
        tab_scan[grep("MaslinerOutput", tab$gpr,
                      ignore.case = TRUE)] <- "Masliner"
    }
    
    ## read in all data files
    assay_table <- mapply(readPBM, gpr_path = tab$gpr, gpr_type = tab_scan,
                          useMean = useMean, filterFlags = filterFlags,
                          readBackground = readBackground, SIMPLIFY = FALSE)

    ## handle rawdata files w/ ID, Name, Sequence information
    rowdat_extra <- NULL
    if (any(tab_scan == "RawData")) {
        gpr_cols <- lapply(assay_table, names)
        extra_cols <- lapply(gpr_cols, setdiff, c("Column", "Row", "foreground", "background"))
        extra_cols <- unique(extra_cols)
        
        ## check if same across samples
        if (length(extra_cols) > 1) {
            stop("If files of type RawData are read, they must all have the same ",
                 "set of 'ID' and 'Name' columns. \nConsider reading in ",
                 "RawData files separately if necessary.")
        }
        extra_cols <- extra_cols[[1]]
        
        ## make sure that assay meta data are equal across samples
        if (length(extra_cols) > 0) {
            rowdat_extra <- lapply(assay_table, `[`, c("Column", "Row", extra_cols))
            if (length(tab_scan) > 1) {
                extra_eq <- lapply(rowdat_extra[-1], dplyr::all_equal, rowdat_extra[[1]])
                extra_eq <- sapply(extra_eq, isTRUE) 
                if (!all(extra_eq)) {
                    stop("Not all file of type RawData have same '",
                         paste(extra_cols[[1]], collapse = "' and '"), "' column values.")
                }
            }
            rowdat_extra <- rowdat_extra[[1]]
        }
    }
    
    
    ## merge GPR data across samples
    if (readBackground) {
        assay_btable <- lapply(assay_table, `[`, c("Column", "Row", "background"))
        assay_btable <- purrr::reduce(assay_btable, left_join, by = c("Column", "Row"))
    }
    assay_table <- lapply(assay_table, `[`, c("Column", "Row", "foreground"))
    assay_table <- purrr::reduce(assay_table, left_join, by = c("Column", "Row"))

    ## convert GPR data to list of DataFrames for SummarizedExperiment assay slot
    assaydat <- list(fore = DataFrame(dplyr::select(assay_table, -Column, -Row)))
    names(assaydat[["fore"]]) <- paste0("s", 1:ncol(assaydat[["fore"]]))
    if (readBackground & !all(tab_scan == "RawData")) {
        assaydat[["back"]] <- DataFrame(dplyr::select(assay_btable, -Column, -Row))
        names(assaydat[["back"]]) <- names(assaydat[["fore"]])
    }
    
    ## row/probe-level metadata from GPR files
    rowdat <- dplyr::select(assay_table, Column, Row)
    if (!is.null(rowdat_extra)) {
        rowdat <- dplyr::left_join(rowdat, rowdat_extra, by = c("Column", "Row"))
    }

    ## check if probe sequences read in from rawdata
    if (!is.null(probes) & !is.null(rowdat_extra)) {
        if ("Sequence" %in% names(rowdat)) {
            stop("'probes' can not be specified because RawData files already include ",
                 "Sequence information. Consider reading in RawData files separately from ",
                 "other raw GPR files.")
        }
    }
    ## check if probe sequences are valid
    if (!is.null(probes)) {
        if (!is.null(rowdat_extra) & "Sequence" %in% names(rowdat)) {
            stop("'probes' can not be specified because RawData files already include ",
                 "Sequence information. Consider reading in RawData files separately from ",
                 "other raw GPR files.")
        }
        if (is.vector(probes, mode = "character")) {
            probes <- DataFrame(Sequence = probes)
        }
        if (!is(probes, "DataFrame") & !is(probes, "data.frame")) {
            warning("Specified 'probes' must be a DataFrame, data.frame, or character vector ",
                    "of probe sequences equal to the size of each GPR array result.\n",
                    "Ignoring specified 'probes' input.")
            probes <- NULL
        } else if (nrow(probes) != nrow(rowdat)) {
            warning("Dimension of specified 'probes' does not match GPR files.\n",
                    "Ignoring specified 'probes' input.")
            probes <- NULL
        } else if (! "Sequence" %in% names(probes)) {
            warning("Specified 'probes' must have a column named 'Sequence' with the probe sequences.\n",
                    "Ignoring specified 'probes' input.")
            probes <- NULL
        }
    }

    ## add probes to metadata if provided
    if (!is.null(probes)) {
        ovnames <- intersect(names(rowdat), names(probes))
        if (length(ovnames) == 0) {
            warning("Specified 'probes' does not have columns matching GPR row metadata.\n",
                    "The probes will be added to the table in the order provided.\n",
                    "WARNING: This may lead to incorrect joining.")
            rowdat <- cbind(rowdat, probes)
        } else {
            if (!all(c("Column", "Row") %in% ovnames)) {
                warning("Specified 'probes' does not have 'Column' and 'Row' columns.\n",
                        "The probes will be matched using other columns.\n",
                        "WARNING: This may lead to incorrect joining.")
            }
            rowdat <- dplyr::left_join(rowdat, tibble::as_tibble(as.data.frame(probes)),
                                       by = ovnames)
        }
    }
    
    ## column/condition-level metadata
    coldat <- DataFrame(dplyr::select(tab, -gpr))
    rownames(coldat) <- names(assaydat[["fore"]])
    
    ## SummarizedExperiment
    SummarizedExperiment(assays = assaydat, rowData = rowdat, colData = coldat,
                         metadata = list(steps = list(),
                                         spatialAdjustment = NULL,
                                         backgroundCorrection = NULL,
                                         betweenArrayNormalization = NULL))
}


#' Read PBM Data File
#' 
#' Helper function for reading in a single PBM raw data file either
#' in GPR format or processed format, e.g. such as files uploaded
#' to the UniPROBE database.
#' If \code{gpr_type = "RawData"} is specified, \code{\link{readRawData}}
#' is used to read in the specified file. Otherwise, \code{\link{readGPR}}
#' is used. Parameters are passed directly to the appropriate function.
#'
#' @param gpr_path path to GPR file.
#' @param gpr_type scan type; one of "Alexa", "Cy3", "Masliner", "RawData".
#' @param useMean logical whether to use mean fluorescent intensity
#'        for each probe rather than median fluorescent intensity.
#'        Ignored if \code{gpr_type = "RawData"}. (default = FALSE)
#' @param filterFlags logical whether to replace intensity values at probes
#'        flagged manually or automatically as being low quality.
#'         ('Bad': -100, 'Absent': -75, 'Not Found': -50) with NA. (default = TRUE)
#' @param readBackground logical whether to also read in probe background
#'        intensities. Ignored if \code{gpr_type = "RawData"}. (default = TRUE)
#'
#' @return
#' See documentation for the appropriate function for more details.
#'
#' @seealso readGPR readRawData
#' @author Patrick Kimes
readPBM <- function(gpr_path, gpr_type, useMean = FALSE, filterFlags = TRUE,
                    readBackground = TRUE) {
    if (gpr_type == "RawData") {
        dat <- readRawData(gpr_path, filterFlags)
    } else if (gpr_type %in% c("Alexa", "Cy3", "Masliner")) {
        dat <- readGPR(gpr_path, gpr_type, useMean = FALSE, filterFlags = TRUE,
                       readBackground = TRUE)
    } else {
        stop("gpr_type must be one of: \"Alexa\", \"Cy3\", \"Masliner\", \"RawData\".")
    }
    return(dat)
}


#' Read GPR File
#' 
#' Helper function for reading in a single GPR file.
#'
#' @param gpr_path path to GPR file.
#' @param gpr_type scan type; one of "Alexa", "Cy3", "Masliner".
#' @param useMean logical whether to use mean fluorescent intensity
#'        for each probe rather than median fluorescent intensity.
#'        (default = FALSE)
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
#' @importFrom readr read_tsv read_lines
#' @importFrom dplyr select
#' @author Patrick Kimes
readGPR <- function(gpr_path, gpr_type, useMean = FALSE, filterFlags = TRUE,
                    readBackground = TRUE) {
    
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

    ## specify which columns to read in
    colt <- rep("-", p)
    colt[c(2, 3)] <- 'i'
    if (filterFlags) {
        colt[grep("\"Flags\"", header)] <- 'i'
    }
    
    ## column corresponding to intensities
    if (useMean) {
        fore_idx <- grep("^\"F.* Mean\"$", header)
        back_idx <- grep("^\"B.* Mean\"$", header)
    } else {
        fore_idx <- grep("^\"F.* Median\"$", header)
        back_idx <- grep("^\"B.* Median\"$", header)
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


#' Read Raw Data File
#' 
#' Helper function for reading in a partially "cleaned" 'rawdata' file, e.g.
#' available on the UniPROBE database.
#' Raw GenePix scan output files typically include a 30 to 40 line header, as well
#' as a large number of columns with various probe-level metrics. However, data uploaded
#' to the UniPROBE database typically include a small subset of the columns, and no
#' header. This function is written to read in raw GPR data stored in the latter format.
#' Files in the first format can be read in using the \code{readGPR} function.
#'
#' @param gpr_path path to text file containing GPR data. See above description
#'        for more information on how this differs from \code{readGPR}. 
#' @param filterFlags logical whether to replace intensity values at probes
#'        flagged manually or automatically as being low quality.
#'        ('Bad': -100, 'Absent': -75, 'Not Found': -50) with NA.
#'        (default = TRUE)
#' @param readBackground logical whether output should include a column of NAs
#'        for background intensities. Only useful for matching column names with
#'        raw GPR data. (default = TRUE)
#' 
#' @return
#' tibble (data.frame-like) object of a single raw data file with between three to six
#' columns: 'Column', 'Row', 'foreground', and if available, 'ID', 'Name', 'Sequence'.
#' Unlike with raw GPR files, the ID and Name columns are assumed to be correct in
#' raw data files.
#' 
#' @details
#' Columns included in \cpde{rawdata} files on UniPROBE can vary substantially.
#' Based on a scan of data sets available on the database, we assume that the
#' Alexa intensity columns correspond to median background subtracted foreground
#' intensities. The first column containing the string "flag" following the
#' Alexa intensity column is assumed to be the corresponding flags for the
#' Alexa scan.
#'
#' @importFrom readr read_tsv read_lines
#' @importFrom dplyr select
#' @author Patrick Kimes
readRawData <- function(gpr_path, filterFlags = TRUE, readBackground = TRUE) {
    
    ## number of rows to skip appears to be variable - determine from reading raw
    header <- readr::read_lines(gpr_path, n_max = 1)

    ## verify that necessary columns are labeled
    if (!grepl(".*Column.*Row.*ID.*", header, ignore.case = TRUE)) {
        stop("header line does not include necessary ",
             "'Column', 'Row', and 'ID' columns.")
    }
    header <- strsplit(header, "\\t")[[1]]

    ## determine number of columns
    p <- length(header)

    ## determine column indicies
    icol <- grep("column", header, ignore.case = TRUE) 
    irow <- grep("row", header, ignore.case = TRUE)
    ival <- grep("alexa", ifelse(grepl("flag", header, ignore.case = TRUE),
                                 NA, header), ignore.case = TRUE)
    iid <- grep("^id$", header, ignore.case = TRUE)
    iname <- grep("name", header, ignore.case = TRUE)
    iseq <- grep("sequence", header, ignore.case = TRUE)
    
    ## make sure column, row, value indicies only occur once
    stopifnot(sapply(list(icol, irow, ival), length) == 1)
    ## make sure id, name, sequence indicies at most once
    stopifnot(sapply(list(iid, iname, iseq), length) <= 1)
    
    ## subset flag column to be after alexa value column
    if (filterFlags) {
        iflag <- grep("flag", header, ignore.case = TRUE)
        iflag <- iflag[iflag > ival]
        if (length(iflag) == 0) {
            warning("No valid 'Flag' column was found. Skipping filtering.")
            filterFlags <- FALSE
        } else if (length(iflag) > 1) {
            warning("More than one 'Flag' column was found. Skipping filtering.")
            filterFlags <- FALSE
        }
    }

    ## specify which columns to read in
    colt <- rep("-", p)
    colt[c(icol, irow)] <- 'i'
    colt[ival] <- 'd'
    if (filterFlags) {
        colt[iflag] <- 'i'
    }
    if (length(iid) > 0) {
        colt[iid] <- 'c'
    }
    if (length(iname) > 0) {
        colt[iname] <- 'c'
    }
    if (length(iseq) > 0) {
        colt[iseq] <- 'c'
    }

    colts <- paste(colt, collapse = "")
    suppressWarnings(
        vals <- readr::read_tsv(gpr_path, col_names = TRUE, col_types = colts,
                                progress = FALSE)
    )
    names(vals)[which(which(colt != "-") == ival)] <- 'foreground'
    names(vals)[which(which(colt != "-") == icol)] <- 'Column'
    names(vals)[which(which(colt != "-") == irow)] <- 'Row'
    if (length(iid) > 0) {
        names(vals)[which(which(colt != "-") == iid)] <- 'ID'
    }
    if (length(iname) > 0) {
        names(vals)[which(which(colt != "-") == iname)] <- 'Name'
    }
    if (length(iseq) > 0) {
        names(vals)[which(which(colt != "-") == iseq)] <- 'Sequence'
    }
    if (filterFlags) {
        names(vals)[which(which(colt != "-") == iflag)] <- 'Flags'
    }

    ## add NA column of 'background' intensities
    if (readBackground) {
        vals$background <- NA_real_
    }
    
    ## remove negative (low quality) flagged probes
    if (filterFlags) {
        vals$foreground[vals$Flags < 0] <- NA_real_
        vals <- dplyr::select(vals, -Flags)
    }
    vals
}

