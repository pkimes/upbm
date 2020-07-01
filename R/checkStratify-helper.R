.pbmCheckStratify <- function(s, strat, bl, gp = NULL, needbl = FALSE, verb = FALSE) {
    
    ## check validity of stratifying colData column
    stopifnot(strat %in% names(colData(s)))
    stopifnot(is.null(gp) || gp %in% names(colData(s)))

    if (is.null(gp)) {
        sg <- dplyr::tibble(sample = colnames(s),
                            Stratify = colData(s)[[strat]],
                            Group = "noGroups")
    } else {
        sg <- dplyr::tibble(sample = colnames(s),
                            Stratify = colData(s)[[strat]],
                            Group = colData(s)[[gp]])
    }

    ## determine baseline if necessary; check validity
    if (is.null(bl)) {
        bl <- grep("ref$", unique(sg$Stratify), value = TRUE, ignore.case = TRUE)
        if (verb && length(bl) > 1) {
            warning("Too many candidate baseline states in '", strat, "' column: ",
                    paste0(bl, collapse = ", "), ".\n",
                    "Using first match: ", bl[1], ".\n",
                    "If other value should be used, specify correct baseline condition w/ 'baseline='.")
            bl <- bl[1]
        }
    } else if (! bl %in% sg$Stratify) {
        stop(bl, " is not a value in '", strat, "' column.\n",
             "Specify correct baseline condition w/ 'baseline='.")
    }
    
    sg_tab <- dplyr::count(sg, Stratify, Group)
    sg_tab <- tidyr::pivot_wider(sg_tab, names_from = Stratify, values_from = n)
    
    if (needbl && any(is.na(sg_tab[[bl]]))) {
        stop("Baseline condition [", bl, "] is not present in all groups.\n",
             "Missing from groups: ",
             paste(sg_tab$Group[is.na(sg_tab[[bl]])], collapse = ", "))
    }
    
    ## determine baseline sample for each group
    sg_bl <- dplyr::filter(sg, Stratify == !! bl)
    sg_bl <- dplyr::group_by(sg_bl, Group)
    if (verb && min(sg_tab[[bl]], na.rm = TRUE) > 1L) {
        warning("Too many samples with baseline states in some groups.\n",
                "Multiple samples in groups: ",
                paste(sg_tab$Group[na.omit(sg_tab[[bl]] > 1L)], collapse = ", "), "\n",
                "Using first sample (by column name) in each group as baseline sample.")
    }
    sg_bl <- dplyr::top_n(sg_bl, 1L, sample)
    sg_bl <- dplyr::ungroup(sg_bl)
    sg_bl <- dplyr::mutate(sg_bl, isBaseline = TRUE)

    ## add baseline status to table of samples
    sg <- dplyr::left_join(sg, sg_bl, by = c("sample", "Stratify", "Group"))
    sg <- tidyr::replace_na(sg, list(isBaseline = FALSE))
    
    ## remove group labels if just using dummy variables
    if (is.null(gp)) {
        sg <- dplyr::select(sg, -Group)
    }

    return(list(coldat = sg, baseline = bl))
}
