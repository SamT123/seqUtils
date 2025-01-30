#' Wrapper for Biostrings::translate
#'
#' Accepts character vector, and handles deletions
#'
#' @param sequences Character vector of DNA/RNA sequences
#' @param reference_aas (optional) Sequences are aligned to `reference` if provided. Useful if, for example, there are deletions in sequences
#'
#' @return character vector of AA sequences, with names preserved
#' @importFrom Biostrings translate DNAStringSet
#' @importFrom stringr str_detect str_remove_all
#'
#' @export
translate = function(sequences, reference_aas = NULL){

  sequences[substr(sequences, 1, 1) == "-"] = map_chr(
    sequences[substr(sequences, 1, 1) == "-"],
    align_dels_to_codons_for_translation
  )

  if (any(stringr::str_detect(sequences, "-")) & is.null(reference_aas)){
    warning("Deletions present but no reference for alignment provided. Output will be unaligned.")
  }

  sequences = stringr::str_remove_all(sequences, "-")
  sequences = align_end_to_codon(sequences)

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


count_flanking_char = function(string, char, leading = T){
  string = strsplit(string, "")[[1]]
  if (!leading) string = rev(string)
  n = 0
  while(string[n+1] == char){
    n = n+1
  }
  as.integer(n)
}

align_end_to_codon = function(sequences){
  stringr::str_sub(
    sequences,
    1,
    floor(nchar(sequences)/3)*3
  )

}

align_dels_to_codons_for_translation = function(sequence){
  if (is.na(sequence)) return(sequence)
  leading_dels = count_flanking_char(sequence, "-", leading = T)
  substr(sequence, 1, 3*ceiling(leading_dels/3)) = paste0(
    rep("-", 3*ceiling(leading_dels/3)),
    collapse = ""
  )

  sequence
}
