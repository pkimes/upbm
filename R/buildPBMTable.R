#' Build PBM Sample Table
#'
#' Helper function to create table of GPR sample files from specified directory path.
#' Parsing of file names is based on previously observed patterns and may not be exact.
#' 
#' @param gpr_dir path to directory containing 'gpr' 
#'
#' @return
#' tibble with each row corresponding to a sample GPR file.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr left_join filter select
#' @export
#' @author Patrick Kimes
buildPBMTable <- function(gpr_dir) {
    ## find all gpr files
    flfull <- list.files(gpr_dir, pattern = "\\.gpr$", full.names = TRUE, ignore.case = TRUE)

    ## determine Cy3, Alexa488 GPRs
    flcy3 <- grep("_Cy3_", flfull, ignore.case = TRUE)
    fl488 <- grep("_Alexa488_", flfull, ignore.case = TRUE)
    if (length(intersect(flcy3, fl488)) > 0) { 
        stop("Can't tell if some files are Cy3 or Alexa488 scans: \n",
             paste(basename(flfull)[intersect(flcy3, fl488)], collapse = "\n"))
    }
    ## only keep Cy3, Alexa488 GPRs
    flfull <- flfull[sort(c(flcy3, fl488))]
    ## re-determine Cy3, Alexa488 GPRs
    flcy3 <- grep("_Cy3_", flfull, ignore.case = FALSE)
    fl488 <- setdiff(1:length(flfull), flcy3)
    
    ## trim off common tail
    fl <- gsub("\\.gpr", "", basename(flfull))
    ## parse assuming "_" delimiters
    fl <- strsplit(fl, "_")
    ## throw error if any file doesn't match expected format
    if (any(sapply(fl[fl488], length) != 14)) {
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
    fl_lppg <- mapply(`[`, fl, fl_nparts-1)
    fl_lp <- as(gsub("lp([0-9]+)pg([0-9]+)", "\\1", fl_lppg), "numeric")
    fl_pg <- as(gsub("lp([0-9]+)pg([0-9]+)", "\\2", fl_lppg), "numeric")

    ## extract condition names from file names (shift +3 for: date,vers,id)
    fl_name <- rep("Cy3", length(flfull))
    fl_name[fl488] <- mapply(`[`, fl[fl488], fl_idx[fl488] + 3)

    fl_scan <- rep("Cy3", length(flfull))
    fl_scan[fl488] <- "Alexa488"
    
    tab <- tibble::tibble(date = fl_date, vers = fl_vers, id = fl_id, idx = fl_idx,
                          condition = fl_name, lp = fl_lp, pg = fl_pg, scan = fl_scan,
                          gpr = flfull)

    ## verify that version/id of experiment gives unique name:idx mapping 
    idx2name <- distinct(dplyr::filter(tab, scan == "Alexa488"), vers, id, idx, condition)
    if (nrow(idx2name) != nrow(distinct(idx2name, vers, id, idx))) {
        stop("Some version/id/idx values map to multiple condition names")
    }
    if (nrow(idx2name) != nrow(distinct(idx2name, vers, id, condition))) {
        stop("Some version/id/condition values map to multiple array indices")
    }

    ## fill in Cy3 scan condition names
    bind_rows(
        dplyr::left_join(dplyr::select(filter::filter(tab, scan == "Cy3"), -condition),
                         idx2name, by = c("vers", "id", "idx")),
        dplyr::filter(tab, scan == "Alexa488")
    )
}
