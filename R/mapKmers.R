#' @title Map K-mer Motifs to PBM Probe Sequences
#'
#' @description
#' This helper function takes a PBMDesign object and set of k-mer
#' sequences of interest, and returns a mapping between probes and 
#' k-mer sequences.
#' Since no probe filtering or trimming is performed as part of this
#' function, \code{pbmTrimProbes} and \code{pbmFilterProbes} should
#' be called before passing the PBMDesign object to this function.
#'
#' @param pd a PBMDesign object containing probe sequence information.
#' @param kmers character vector of K-mers (of equal length) to which
#'        probes should be mapped.
#'
#' @return
#' A table of k-mer sequences and probes, with each row
#' corresponding to an occurrence of a k-mer sequence on a probe.
#' The table inclues all columns of the input \code{pd} design
#' as well as columns indicating the orientation (\code{orient}) and
#' position (\code{pos}) of the k-mer in the probe. 
#'
#' @importFrom tibble tibble rownames_to_column
#' @importFrom tidyr unnest
#' @importFrom dplyr mutate filter bind_rows select left_join
#' @importFrom Biostrings DNAStringSet reverseComplement
#' @importFrom tidyselect everything
#' @export
#' @author Patrick Kimes
mapKmers <- function(pd, kmers) {
    stopifnot(is(pd, "PBMDesign"))
    stopifnot(is.vector(kmers, mode = "character"))
    stopifnot("probeID" %in% names(pd@design))
    
    pdd <- pd@design
    
    ## check if kmers/motfs are of uniform length
    k <- unique(nchar(kmers))
    if (length(k) > 1) {
        stop("Specified 'kmers' includes motifs of varying length.\n",
             "Method currently only supports kmers of uniform length.")
    }
    
    ## check if sequences are of uniform length
    seql <- unique(nchar(pdd$Sequence))
    if (length(seql) > 1) {
        stop("Probe sequences should all be of equal length.") 
    }
    
    ## count up occurrences of kmers - note: will do ALL kmers
    rolls <- lapply(pdd$Sequence, substring, 1:(seql - k + 1), k:seql)
    rolls <- tibble::tibble(probeID = pdd$probeID, fwd_seq = rolls)

    ## add start position (1-indexed)
    rolls <- dplyr::mutate(rolls, pos = list(1:(seql - k + 1)))

    ## unnest and add reverse complement sequence
    rolls <- tidyr::unnest(rolls, tidyselect::everything())
    rolls <- dplyr::mutate(rolls, rev_seq = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(fwd_seq))))
    
    rolls_fwd <- dplyr::filter(rolls, fwd_seq %in% kmers)
    rolls_fwd <- dplyr::mutate(rolls_fwd, seq = fwd_seq, orient = 'fwd') 
    
    rolls_rev <- dplyr::filter(rolls, rev_seq %in% kmers, ! rev_seq == fwd_seq)
    rolls_rev <- dplyr::mutate(rolls_rev, seq = rev_seq, orient = 'rev') 
    
    rolls <- dplyr::bind_rows(rolls_fwd, rolls_rev)
    rolls <- dplyr::select(rolls, probeID, seq, pos, orient)
    
    ## add all probe identifiers except for sequence
    if (length(pdd) > 1) {
        rolls <- dplyr::left_join(rolls, pdd, by = "probeID")
    }
    
    return(rolls)
}
