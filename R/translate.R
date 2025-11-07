#' Wrapper for Biostrings::translate with deletion handling
#'
#' Translates DNA/RNA sequences to amino acids, with special handling for deletions
#' (gap characters "-"). The function processes deletions as follows:
#'
#' 1. Leading deletions are padded to codon boundaries (multiples of 3)
#' 2. All deletions are removed before translation
#' 3. If `reference_aas` is provided, translated sequences are realigned to the reference
#' 4. If deletions exist without a reference, output will be unaligned (warning issued)
#'
#' @param sequences Character vector of DNA/RNA sequences. May contain gap characters ("-").
#' @param reference_aas (optional) Reference amino acid sequence. If provided, translated
#'   sequences will be aligned to this reference using MAFFT. This is essential when
#'   input sequences contain deletions and you want to preserve alignment in the output.
#'
#' @return Character vector of amino acid sequences, with names preserved. If `reference_aas`
#'   is provided, sequences will be aligned. Otherwise, deletions are simply removed and
#'   sequences may have different lengths.
#'
#' @importFrom Biostrings translate DNAStringSet
#' @importFrom stringr str_detect str_remove_all
#'
#' @export
translate = function(sequences, reference_aas = NULL) {
  if (any(is.na(sequences))) {
    warning("Some sequences are NA!")
  }

  do_translation = function(sequences, reference_aas) {
    sequences[substr(sequences, 1, 1) == "-"] = map_chr(
      sequences[substr(sequences, 1, 1) == "-"],
      align_dels_to_codons_for_translation
    )

    if (any(stringr::str_detect(sequences, "-")) & is.null(reference_aas)) {
      stop(
        "Deletions present but no reference for alignment provided. Output will be unaligned."
      )
    }

    sequences = stringr::str_remove_all(sequences, "-")
    sequences = align_end_to_codon(sequences)

    aa_sequences = Biostrings::translate(
      Biostrings::DNAStringSet(sequences),
      if.fuzzy.codon = "X"
    )

    aa_sequences = as.character(aa_sequences)
    names(aa_sequences) = names(sequences)

    if (!is.null(reference_aas)) {
      aa_sequences = mafft_align(aa_sequences, reference_aas)
    }

    aa_sequences
  }

  aa_sequences = setNames(
    rep(NA, length(sequences)),
    names(sequences)
  )

  aa_sequences[!is.na(sequences)] = do_translation(
    sequences[!is.na(sequences)],
    reference_aas
  )

  aa_sequences
}


count_flanking_char = function(string, char, leading = T) {
  string = strsplit(string, "")[[1]]
  if (!leading) {
    string = rev(string)
  }
  n = 0
  while (string[n + 1] == char) {
    n = n + 1
  }
  as.integer(n)
}

align_end_to_codon = function(sequences) {
  stringr::str_sub(
    sequences,
    1,
    floor(nchar(sequences) / 3) * 3
  )
}

align_dels_to_codons_for_translation = function(sequence) {
  if (is.na(sequence)) {
    return(sequence)
  }
  leading_dels = count_flanking_char(sequence, "-", leading = T)
  substr(sequence, 1, 3 * ceiling(leading_dels / 3)) = paste0(
    rep("-", 3 * ceiling(leading_dels / 3)),
    collapse = ""
  )

  sequence
}
