#' Make PBM Sample Table
#'
#' Helper function to create table of GPR sample files from specified directory path.
#' Parsing of file names is based on previously observed patterns and may not be exact.
#' 
#' @param gpr_dir path to directory containing GPR files.
#' @param gpr_type string value specifying type of GPR files to find; must
#'        be one of: "Alexa488", "Cy3", "Masliner", "all". (default = "Alexa488")
#' 
#' @return
#' tibble with each row corresponding to a sample GPR file.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr left_join filter select
#' @export
#' @author Patrick Kimes
makePBMTable <- function(gpr_dir, gpr_type = "Alexa488") {
    stopifnot(gpr_type %in% c("Alexa488", "Masliner", "Cy3", "all"))
    
    ## find all GPR file types
    flfull <- list.files(gpr_dir, pattern = "\\.gpr$", full.names = TRUE, ignore.case = TRUE)

    ## determine Cy3, Alexa488 GPRs
    flcy3 <- grep("_Cy3_", basename(flfull), ignore.case = TRUE)
    fl488 <- grep("_Alexa488_(?!.*MaslinerOutput)", basename(flfull), ignore.case = TRUE, perl = TRUE)
    flmas <- grep("MaslinerOutput", basename(flfull), ignore.case = TRUE)
    
    if (length(intersect(flcy3, fl488)) > 0) { 
        stop("Can't tell if some files are Alexa488 or Cy3 scans: \n",
             paste(basename(flfull)[intersect(flcy3, fl488)], collapse = "\n"))
    }
   if (length(intersect(flmas, fl488)) > 0) { 
        stop("Can't tell if some files are Alexa488 or Masliner scans: \n",
             paste(basename(flfull)[intersect(flmas, fl488)], collapse = "\n"))
    }
   if (length(intersect(flcy3, flmas)) > 0) { 
        stop("Can't tell if some files are Cy3 or Masliner scans: \n",
             paste(basename(flfull)[intersect(flcy3, fl488)], collapse = "\n"))
    }
    if (length(c(flcy3, fl488, flmas)) < length(flfull)) {
        fldrop <- setdiff(1:length(flfull), c(flcy3, fl488, flmas))
        warning("Can't tell if following files were Alexa488, Cy3, or Masliner scans: \n",
                paste(basename(flfull)[fldrop], collapse = "\n"))
    }
        
    fl_type <- rep("", length(flfull))
    fl_type[flcy3] <- "Cy3"
    fl_type[fl488] <- "Alexa488"
    fl_type[flmas] <- "Masliner"
    
    ## check for any ambiguities if masliner or cy3 desired
    if (gpr_type == "all") {
        flfull <- flfull[fl_type != ""]
        fl_type <- fl_type[fl_type != ""]
    } else {
        flfull <- flfull[fl_type == gpr_type]
        fl_type <- rep(gpr_type, length(flfull))
    }
    
    ## trim off common tail
    fl <- gsub("\\.gpr", "", basename(flfull), ignore.case = TRUE)
    ## parse assuming "_" delimiters
    fl <- strsplit(fl, "_")
    ## throw error if any file doesn't match expected format
    if (any(sapply(fl[fl_type != "Cy3"], length) != 14)) {
        stop("file name doesn't match expectation")
    }

    ## cy3 GPR file names can be of varying length
    fl_nparts <- sapply(fl, length)
    
    ## metadata for experiment run 
    fl_date <- sapply(fl, `[`, 1)
    fl_vers <- tolower(sapply(fl, `[`, 2))
    fl_id <- sapply(fl, `[`, 3)

    ## index of sample/experiment on 8-plexed chip
    fl_idx <- mapply(`[`, fl, fl_nparts)
    fl_idx <-  as(gsub("(.*)-8", "\\1", fl_idx), "numeric")

    ## laser power and power gain
    fl_lp <- rep(NA_real_, length(flfull))
    fl_pg <- rep(NA_real_, length(flfull))
    if (any(fl_type != "Masliner")) {
        fl_lppg <- mapply(`[`, fl[fl_type != "Masliner"], fl_nparts[fl_type != "Masliner"] - 1)
        fl_lp[fl_type != "Masliner"] <- as(gsub("lp([0-9]+)pg([0-9]+)", "\\1", fl_lppg), "numeric")
        fl_pg[fl_type != "Masliner"] <- as(gsub("lp([0-9]+)pg([0-9]+)", "\\2", fl_lppg), "numeric")
    }
    
    ## extract condition names from file names (shift +3 for: date,vers,id)
    fl_name <- rep("Cy3", length(flfull))
    if (any(fl_type != "Cy3")) {
        fl_name[fl_type != "Cy3"] <- mapply(`[`, fl[fl_type != "Cy3"], fl_idx[fl_type != "Cy3"] + 3)
    }

    ## construct primary sample table 
    tab <- tibble::tibble(date = fl_date, vers = fl_vers, id = fl_id, idx = fl_idx,
                          condition = fl_name, lp = fl_lp, pg = fl_pg, scan = fl_type,
                          gpr = flfull)

    ## verify that version/id of experiment gives unique name:idx mapping 
    idx2name <- dplyr::distinct(dplyr::filter(tab, scan == "Alexa488"), vers, id, idx, condition)
    if (nrow(idx2name) != nrow(dplyr::distinct(idx2name, vers, id, idx))) {
        stop("Some version/id/idx values map to multiple condition names")
    }
    if (nrow(idx2name) != nrow(dplyr::distinct(idx2name, vers, id, condition))) {
        stop("Some version/id/condition values map to multiple array indices")
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
