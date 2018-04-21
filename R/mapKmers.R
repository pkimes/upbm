#' Map K-mer Motifs to PBM Probe Sequences
#'
#' This helper function takes a standardized PBM sequence design file,
#' e.g. "8x60k_v1_sequences.txt", and k-mer sequences of interest,
#' and identifies the subset of probes containing each k-mer sequence.
#' The output table can be used to compute k-mer-level affinities from
#' the probe-level intensity values.
#'
#' @param probes table containing sequence information for probes,
#'        e.g. read in from "8x60k_v1_sequences.txt". At a minimum,
#'        must have a column called "Sequence" with the array probe
#'        sequences.
#' @param kmers character vector of K-mers (of equal length) to which
#'        probes should be mapped.
#'
#' @return
#' a table with 8mer sequences (seq), and probes, where each row
#' corresponds to a unique occurrence of the sequence on a probe and on
#' the array. The table inclues all columns of the input \code{probes}
#' table except for the \code{Sequences} column, as well as columns
#' indicating the orientation (orient) and position (pos) of the k-mer
#' in the probe. 
#'
#' @importFrom tibble tibble rownames_to_column
#' @importFrom tidyr unnest
#' @importFrom dplyr mutate filter bind_rows select left_join
#' @importFrom Biostrings DNAStringSet reverseComplement
#' @export
#' @author Patrick Kimes
mapKmers <- function(probes, kmers) {
    ## check validity of inputs
    if (is(probes, "DataFrame")) {
        probes <- as.data.frame(probes)
    } else if (!is(probes, "data.frame")) {
        stop("Specified 'probes' must be a DataFrame or data.frame ",
             "of probe sequences.")
    }
    stopifnot("Sequence" %in% names(probes))
    stopifnot(is.vector(kmers, mode = "character"))

    ## check if kmers/motfs are of uniform length
    k <- unique(nchar(kmers))
    if (length(k) > 1) {
        stop("Specified 'kmers' includes motifs of varying length.\n",
             "Method currently only supports kmers of uniform length.")
    }

    if (! "probe_idx" %in% colnames(probes)) {
        probes <- dplyr::mutate(probes, probe_idx = 1:n())
    }  else {
        if (any(duplicated(probes$probe_idx))) {
            stop("Probe DataFrame has 'probe_idx' column with non-unique entries.\n",
                 "If 'probe_idx' is specified, entries must be unique, ",
                 "or column should be removed.")
        }
    }
    
    ## check if sequences are of uniform length
    seql <- unique(nchar(probes$Sequence))
    if (length(seql) > 1) {
        stop("Probe sequences should all be of equal length.") 
    }
    
    ## count up occurrences of kmers - note: will do ALL kmers
    rolls <- lapply(probes$Sequence, substring, 1:(seql - k + 1), k:seql)
    rolls <- tibble::tibble(probe_idx = probes$probe_idx, fwd_seq = rolls)

    ## add start position (1-indexed)
    rolls <- dplyr::mutate(rolls, pos = list(1:(seql - k + 1)))

    ## unnest and add reverse complement sequence
    rolls <- tidyr::unnest(rolls)
    rolls <- dplyr::mutate(rolls, rev_seq = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(fwd_seq))))
    
    rolls_fwd <- dplyr::filter(rolls, fwd_seq %in% kmers)
    rolls_fwd <- dplyr::mutate(rolls_fwd, seq = fwd_seq, orient = 'fwd') 
    
    rolls_rev <- dplyr::filter(rolls, rev_seq %in% kmers, ! rev_seq == fwd_seq)
    rolls_rev <- dplyr::mutate(rolls_rev, seq = rev_seq, orient = 'rev') 
    
    rolls <- dplyr::bind_rows(rolls_fwd, rolls_rev)
    rolls <- dplyr::select(rolls, probe_idx, seq, pos, orient)
    
    ## add all probe identifiers except for sequence
    if (length(probes) > 1) {
        rolls <- dplyr::left_join(rolls, dplyr::select(probes, -Sequence), by = "probe_idx")
    }
    
    return(rolls)
}
