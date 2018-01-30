#' Map 8-mer Motifs to PBM Probe Sequences
#'
#' This helper function takes a standardized PBM sequence design file,
#' e.g. "8x60k_v1_sequences.txt", and 8-mer sequences of interest,
#' and identifies the subset of probes containing each 8-mer sequence.
#' The output table can be used to compute 8-mer-level affinities from
#' the probe-level intensity values.
#'
#' @param seqfile path to table containing sequence information for probes,
#'        e.g. "8x60k_v1_sequences.txt". At a minimum, must have a column
#'        called "Sequence" with the array probe sequences.
#' @param kmers character vector of 8mers to map probe IDs to.
#'
#' @return
#' a table with 8mer sequences (seq), and probes, where each row
#' corresponds to a unique occurrence of the sequence on a probe and on
#' the array.
#'
#' @author Patrick Kimes
map8mers <- function(seqfile, kmers) {
    
    ## map from chip coordinates (Column, Row) to 60nt probe sequence
    probes <- read_tsv(seqfile)

    ## only keep de Bruijn probes and trim out primer sequence
    probes <- probes %>%
        filter(grepl("^dBr_", NAME)) %>%
        mutate(Sequence = subseq(Sequence, 1, 36),
               probe_idx = row_number())

    ## count up occurrences of 8mers
    rolls <- lapply(probes$Sequence, substring, 1:29, 8:36)
    rolls <- tibble(probe_idx = probes$probe_idx, fwd_seq = rolls)
    rolls <- unnest(rolls)
    rolls <- mutate(rolls, rev_seq = as.character(reverseComplement(DNAStringSet(fwd_seq))))
    
    rolls_fwd <- filter(rolls, fwd_seq %in% kmers)
    rolls_fwd <- mutate(rolls_fwd, seq = fwd_seq, orient = 'fwd') 
    
    rolls_rev <- filter(rolls, rev_seq %in% kmers, ! rev_seq == fwd_seq)
    rolls_rev <- mutate(rolls_rev, seq = rev_seq, orient = 'rev') 
    
    rolls <- bind_rows(rolls_fwd, rolls_rev)
    rolls <- select(rolls, probe_idx, seq, orient)

    ## add all probe identifiers except for sequence
    rolls <- left_join(rolls, select(probes, -Sequence), by = "probe_idx")

    return(rolls)
}
