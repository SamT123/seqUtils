#' Wrapper for Biostrings::translate
#'
#' Accepts character vector, and handles deletions
#'
#' @param sequences Character vector of DNA/RNA sequences
#' @param reference (optional) Sequences are aligned to `reference` if provided. Useful if, for example, there are deletions in sequences
#'
#' @return character vector of AA sequences, with names preserved
#' @importFrom Biostrings translate DNAStringSet
#' @importFrom stringr str_detect str_remove_all
#'
#' @export
translate = function(sequences, reference_aas = NULL){

  any_deletions = any(stringr::str_detect(sequences, "-"))

  if (any_deletions & is.null(reference_aas)){
    warning("Deletions detected & no reference for alignment provided. Output will be unaligned.")
  }

  sequences = stringr::str_remove_all(sequences, "-")

  aa_sequences = Biostrings::translate(
    Biostrings::DNAStringSet(sequences),
    if.fuzzy.codon = "X"
  )

  aa_sequences = as.character(aa_sequences)
  names(aa_sequences) = names(sequences)

  if (!is.null(reference_aas)){
    aa_sequences = mafft_align(aa_sequences, reference_aas)
  }

  aa_sequences
}
