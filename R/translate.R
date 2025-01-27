#' Wrapper for Biostrings::translate so a character vector can be provided
#'
#' @param sequences character vector of DNA/RNA sequences
#'
#' @return character vector of AA sequences, with names preserved
#' @importFrom Biostrings translate DNAStringSet
#'
#' @export
translate = function(sequences){
  aa_sequences = Biostrings::translate(
    Biostrings::DNAStringSet(sequences),
    if.fuzzy.codon = "X"
  )

  aa_sequences = as.character(aa_sequences)
  names(aa_sequences) = names(sequences)

  aa_sequences
}
