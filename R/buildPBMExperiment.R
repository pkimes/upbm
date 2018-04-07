#' Read PBMExperiment Data
#'
#' Read PBM data from a table containing paths to GPR files. 
#' 
#' @param tab table of samples with at least `gpr` and `vers` columns
#'        corresponding to GPR file paths and the corresponding PBM
#'        design version.
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
    stopifnot(c("vers", "gpr") %in% names(tab))
    stopifnot(is.logical(useMean))
    stopifnot(is.logical(filterFlags))
    
    ## currently only support all scans with same design
    if (length(unique(tab$vers)) > 1) {
        stop("All samples/scans must have the same assay design version")
    }

    ## guess scan type if not specified - only care if Masliner or not
    if ("scan" %in% names(tab)) {
        tab_scan <- tab$scan
    } else {
        ## don't add to 'tab' since don't want in colData of output
        tab_scan <- rep("Alexa488", nrow(tab))
        tab_scan[grep("MaslinerOutput", tab$gpr,
                      ignore.case = TRUE)] <- "Masliner"
    }
    
    ## read in all GPR scans
    assay_table <- mapply(readGPR, gpr_path = tab$gpr, gpr_type = tab_scan,
                          useMean = useMean, filterFlags = filterFlags,
                          readBackground = readBackground, SIMPLIFY = FALSE)
    
    ## merge GPR data across samples
    if (readBackground) {
        assay_btable <- lapply(assay_table, `[`, c("Column", "Row", "background"))
        assay_btable <- purrr::reduce(assay_btable, left_join, by = c("Column", "Row"))
        assay_table <- lapply(assay_table, `[`, c("Column", "Row", "foreground"))
    }
    assay_table <- purrr::reduce(assay_table, left_join, by = c("Column", "Row"))

    ## convert GPR data to list of DataFrames for SummarizedExperiment assay slot
    assaydat <- list(fore = DataFrame(dplyr::select(assay_table, -Column, -Row)))
    names(assaydat[["fore"]]) <- paste0("s", 1:ncol(assaydat[["fore"]]))
    if (readBackground) {
        assaydat[["back"]] <- DataFrame(dplyr::select(assay_btable, -Column, -Row))
        names(assaydat[["back"]]) <- names(assaydat[["fore"]])
    }

    ## row/probe-level metadata from GPR files
    rowdat <- dplyr::select(assay_table, Column, Row)

    ## check probes
    if (!is.null(probes)) {
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


#' Read GPR File as Assay
#' 
#' Helper function for reading in a single GPR file.
#'
#' @param gpr_path path to GPR file.
#' @param gpr_type scan type; one of "Alexa488", "Cy3", "Masliner".
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
        stop("Foreground column could not be identified!")
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

