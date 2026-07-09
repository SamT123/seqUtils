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
#' If `aligned = TRUE` the input is assumed to be a fixed-column alignment whose
#' length is a multiple of 3. Gaps are then meaningful rather than noise, and no
#' realignment is performed: amino-acid coordinates are stable by construction.
#'
#' @param sequences Character vector of DNA/RNA sequences. May contain gap characters ("-").
#' @param reference_aas (optional) Reference amino acid sequence. If provided, translated
#'   sequences will be aligned to this reference using MAFFT. This is essential when
#'   input sequences contain deletions and you want to preserve alignment in the output.
#' @param aligned If TRUE, treat `sequences` as an in-frame alignment and
#'   translate it column-stably without MAFFT. A codon of three gaps becomes
#'   "-"; a codon with some but not all gaps, or with any ambiguity code,
#'   becomes "X"; a trailing stop codon is kept as "*". Incompatible with
#'   `reference_aas`.
#'
#' @return Character vector of amino acid sequences, with names preserved. If `reference_aas`
#'   is provided, sequences will be aligned. Otherwise, deletions are simply removed and
#'   sequences may have different lengths. If `aligned = TRUE`, every output has
#'   the same length, `nchar(sequences) / 3`.
#'
#' @examples
#' # Simple translation without deletions
#' seqs <- c("ATGGCC", "ATGTTT")
#' translate(seqs)
#'
#' # A fixed-column alignment: gaps are preserved, coordinates are stable
#' translate(c("ATGGCCTAA", "ATG---TAA"), aligned = TRUE)
#'
#' \dontrun{
#' # With deletions and reference alignment
#' seqs <- c("ATGGCC", "ATG---GCC")
#' ref_aa <- translate(seqs[1])
#' translate(seqs, reference_aas = ref_aa)
#' }
#'
#' @importFrom Biostrings translate DNAStringSet
#' @importFrom stringr str_detect str_remove_all
#'
#' @export
translate <- function(sequences, reference_aas = NULL, aligned = FALSE) {
  if (any(is.na(sequences))) {
    warning("Some sequences are NA!")
  }

  if (aligned) {
    if (!is.null(reference_aas)) {
      rlang::abort(
        paste(
          "`reference_aas` is meaningless when `aligned = TRUE`:",
          "an aligned translation preserves coordinates already."
        ),
        class = "seqUtils_error_aligned_with_reference"
      )
    }
    return(translateAligned(sequences))
  }

  do_translation <- function(sequences, reference_aas) {
    sequences[substr(sequences, 1, 1) == "-"] <- map_chr(
      sequences[substr(sequences, 1, 1) == "-"],
      align_dels_to_codons_for_translation
    )

    if (any(stringr::str_detect(sequences, "-")) & is.null(reference_aas)) {
      stop(
        "Deletions present but no reference for alignment provided. Output will be unaligned."
      )
    }

    sequences <- stringr::str_remove_all(sequences, "-")
    sequences <- align_end_to_codon(sequences)

    aa_sequences <- Biostrings::translate(
      Biostrings::DNAStringSet(sequences),
      if.fuzzy.codon = "X"
    )

    aa_sequences <- as.character(aa_sequences)
    names(aa_sequences) <- names(sequences)

    if (!is.null(reference_aas)) {
      aa_sequences <- mafft_align(aa_sequences, reference_aas)
    }

    aa_sequences
  }

  aa_sequences <- setNames(
    rep(NA, length(sequences)),
    names(sequences)
  )

  aa_sequences[!is.na(sequences)] <- do_translation(
    sequences[!is.na(sequences)],
    reference_aas
  )

  aa_sequences
}


translateAligned <- function(sequences) {
  aa_sequences <- setNames(
    rep(NA_character_, length(sequences)),
    names(sequences)
  )

  present <- !is.na(sequences)
  if (!any(present)) {
    return(aa_sequences)
  }

  nt <- toupper(sequences[present])
  width <- unique(nchar(nt))

  if (length(width) != 1L || width %% 3L != 0L) {
    rlang::abort(
      paste0(
        "`aligned = TRUE` requires a fixed-column alignment whose width is a ",
        "multiple of 3; got ",
        paste(sort(width), collapse = ", "),
        "."
      ),
      class = "seqUtils_error_not_an_alignment"
    )
  }

  # Gaps translate as X via N, then the all-gap codons are restored to "-".
  aa <- as.character(Biostrings::translate(
    Biostrings::DNAStringSet(chartr("-", "N", nt)),
    if.fuzzy.codon = "X"
  ))

  codon_starts <- seq(1L, width, by = 3L)

  for (i in which(stringr::str_detect(nt, "-"))) {
    all_gap <- which(
      substring(nt[[i]], codon_starts, codon_starts + 2L) == "---"
    )
    if (length(all_gap) > 0L) {
      residues <- strsplit(aa[[i]], "")[[1]]
      residues[all_gap] <- "-"
      aa[[i]] <- paste0(residues, collapse = "")
    }
  }

  aa_sequences[present] <- aa
  aa_sequences
}


count_flanking_char <- function(string, char, leading = TRUE) {
  string <- strsplit(string, "")[[1]]
  if (!leading) {
    string <- rev(string)
  }
  n <- 0
  while (string[n + 1] == char) {
    n <- n + 1
  }
  as.integer(n)
}

align_end_to_codon <- function(sequences) {
  stringr::str_sub(
    sequences,
    1,
    floor(nchar(sequences) / 3) * 3
  )
}

align_dels_to_codons_for_translation <- function(sequence) {
  if (is.na(sequence)) {
    return(sequence)
  }
  leading_dels <- count_flanking_char(sequence, "-", leading = TRUE)
  substr(sequence, 1, 3 * ceiling(leading_dels / 3)) <- paste0(
    rep("-", 3 * ceiling(leading_dels / 3)),
    collapse = ""
  )

  sequence
}
