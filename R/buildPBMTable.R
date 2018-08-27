#' Build PBM Sample Table
#'
#' Helper function to create table of GPR sample files from specified directory path or
#' manually specified list of GPR file names. One of \code{gpr_dir} or \code{gpr_files}
#' must be specified. Parsing of file names is based on previously observed patterns and may not be exact.
#' Expected naming conventions of GPR files are described in the Details below. 
#' 
#' @param gpr_dir path to directory containing GPR files. (default = NULL)
#' @param gpr_type string value specifying type of GPR files to find; must
#'        be one of: "Alexa", "Cy3", "Masliner", "all". (default = "Alexa")
#' @param gpr_files alternatively, rather than specifying a directory containing
#'        GPR files, can specify a vector of GPR file names. (default = NULL)
#' 
#' @return
#' tibble with each row corresponding to a sample GPR file.
#'
#' @details
#' For the three types of GPR files currently supported, files are assumed to be named
#' according to the following conventions:  
#' * **Alexa**: `{date}_{version}_{runID}_{condition1}_...` \cr `_{condition8}_Alexa[0-9]+_lp{laser}pg{gain}_{index}_8.gpr`
#' * **Masliner**: `{date}_{version}_{runID}_{condition1}_...` \cr `_{condition8}_Alexa[0-9]+_MaslinerOutput_{index}_8.gpr`
#' * **Cy3**: `{date}_{version}_{runID}_{miscellaneous}_Cy3_lp{laser}pg{gain}_{index}_8.gpr`
#'
#' @md
#' @importFrom tibble tibble
#' @importFrom dplyr left_join filter select
#' @export
#' @author Patrick Kimes
buildPBMTable <- function(gpr_dir = NULL, gpr_type = "Alexa", gpr_files = NULL) {
    stopifnot(gpr_type %in% c("Alexa", "Masliner", "Cy3", "all"))
    stopifnot(!is.null(gpr_dir) || !is.null(gpr_files))
    if (!is.null(gpr_dir) && !is.null(gpr_files)) {
        stop("Only one of 'gpr_dir' or 'gpr_files' should be specified.")
    }
    
    ## find all GPR files
    if (!is.null(gpr_dir)) {
        flfull <- list.files(gpr_dir, pattern = "\\.gpr$", full.names = TRUE, ignore.case = TRUE)
    } else {
        if (!all(grepl("\\.gpr$", gpr_files, ignore.case = TRUE))) {
            stop("Files passed to 'gpr_files' must end in '.gpr'.")
        }
        flfull <- gpr_files
    }
    
    ## determine Cy3, Alexa488 GPRs
    flcy3 <- grep("_Cy3_", basename(flfull), ignore.case = TRUE)
    fl488 <- grep("_Alexa[0-9]*_(?!.*MaslinerOutput)", basename(flfull), ignore.case = TRUE, perl = TRUE)
    flmas <- grep("MaslinerOutput", basename(flfull), ignore.case = TRUE)
    
    ## check for any ambiguities if masliner or Cy3 desired
    if (length(intersect(flcy3, fl488)) > 0) { 
        stop("Can't tell if some files are Alexa or Cy3 scans: \n",
             paste(basename(flfull)[intersect(flcy3, fl488)], collapse = "\n"))
    }
   if (length(intersect(flmas, fl488)) > 0) { 
        stop("Can't tell if some files are Alexa or Masliner scans: \n",
             paste(basename(flfull)[intersect(flmas, fl488)], collapse = "\n"))
    }
   if (length(intersect(flcy3, flmas)) > 0) { 
        stop("Can't tell if some files are Cy3 or Masliner scans: \n",
             paste(basename(flfull)[intersect(flcy3, fl488)], collapse = "\n"))
    }
    if (length(c(flcy3, fl488, flmas)) < length(flfull)) {
        fldrop <- setdiff(1:length(flfull), c(flcy3, fl488, flmas))
        warning("Can't tell if following files were Alexa, Cy3, or Masliner scans: \n",
                paste(basename(flfull)[fldrop], collapse = "\n"))
    }
        
    fl_type <- rep("", length(flfull))
    fl_type[flcy3] <- "Cy3"
    fl_type[fl488] <- "Alexa"
    fl_type[flmas] <- "Masliner"
    
    if (gpr_type == "all") {
        flfull <- flfull[fl_type != ""]
        fl_type <- fl_type[fl_type != ""]
    } else {
        flfull <- flfull[fl_type == gpr_type]
        fl_type <- rep(gpr_type, length(flfull))
    }
    
    ## trim off common tail
    fl <- gsub("\\.gpr", "", basename(flfull), ignore.case = TRUE)

    ## peel off date prefix (if available)
    fl_date <- gsub("(^.*?)_?v[[:digit:]]+.*", "\\1", fl, ignore.case = TRUE)
    fl <- gsub("^.*?_?(v[[:digit:]]+.*)", "\\1", fl, ignore.case = TRUE)

    ## peel off assay version (if available)
    fl_vers <- tolower(gsub("(^v?[[:digit:]]+)_?.*", "\\1", fl))
    fl <- gsub("^v?[[:digit:]]+_?(.*)", "\\1", fl)

    ## peel off assay id (if available)
    fl_id <- gsub("(^[[:digit:]_]*)[[:alpha:]]+.*", "\\1", fl)
    fl_id <- gsub("_$", "", fl_id)
    fl <- gsub("^[[:digit:]_]*([[:alpha:]]+.*)", "\\1", fl)

    ## peel off index from tail (must be available)
    if (!all(grepl("[[:digit:]]+-[[:digit:]]+$", fl))) {
        stop("At least one sample file does not have expected index suffix.")
    }
    fl_tot <- as(gsub("^.*?-([[:digit:]]+)$", "\\1", fl), "integer")
    fl_idx <- as(gsub("^.*?([[:digit:]]+)-[[:digit:]]+$", "\\1", fl), "integer")
    fl <- gsub("(^.*?)_?[[:digit:]]+-[[:digit:]]+$", "\\1", fl)

    ## peal off laser power and power gain from tail (if available)
    fl_lp <- rep(NA_real_, length(flfull))
    fl_pg <- rep(NA_real_, length(flfull))
    if (any(fl_type != "Masliner")) {
        fl_lp[fl_type != "Masliner"] <- gsub(".*lp([[:digit:]]+)pg([[:digit:]]+)$",
                                             "\\1", fl[fl_type != "Masliner"])
        fl_pg[fl_type != "Masliner"] <- gsub(".*lp([[:digit:]]+)pg([[:digit:]]+)$",
                                             "\\2", fl[fl_type != "Masliner"])
        fl_lp <- as(fl_lp, "numeric")
        fl_pg <- as(fl_pg, "numeric")

        fl[fl_type != "Masliner"] <- gsub("(.*?)_?lp[[:digit:]]+pg[[:digit:]]+$", "\\1",
                                          fl[fl_type != "Masliner"])
    }
    fl[fl_type == "Masliner"] <- gsub("(.*?)_?Maslin.*", "\\1", fl[fl_type == "Masliner"],
                                      ignore.case = TRUE)

    ## peal off last part of file (should be Alexa or Cy3)
    fl <- gsub("(.*?)_?[[:alnum:]]+$", "\\1", fl)

    ## parse remaining parts (should contain conditions)
    fl <- strsplit(fl, "_")
    fl_nparts <- sapply(fl, length)
    
    ## extract condition names from file names (if available)
    fl_name <- rep("unknown", length(flfull))
    if (any(fl_nparts > 0)) {
        fl_name[fl_nparts > 0] <- mapply(`[`, fl[fl_nparts > 0], fl_idx[fl_nparts > 0])
    }

    ## construct primary sample table 
    tab <- tibble::tibble(date = fl_date, vers = fl_vers, id = fl_id, idx = fl_idx,
                          condition = fl_name, lp = fl_lp, pg = fl_pg, scan = fl_type,
                          gpr = flfull)

    ## verify that version/id of experiment gives unique name:idx mapping 
    idx2name <- dplyr::distinct(dplyr::filter(tab, scan == "Alexa", condition != "unknown"),
                                vers, id, idx, condition)
    if (nrow(idx2name) != nrow(dplyr::distinct(idx2name, vers, id, idx))) {
        stop("Some version/id/idx values map to multiple condition names")
    }

    ## fill in Cy3 scan condition names
    if (any(fl_type == "Cy3")) {
        tab <- dplyr::bind_rows(
                          dplyr::left_join(dplyr::select(dplyr::filter(tab, scan == "Cy3"), -condition),
                                           idx2name, by = c("vers", "id", "idx")),
                          dplyr::filter(tab, scan != "Cy3")
                      )
    }

    return(tab)
}
