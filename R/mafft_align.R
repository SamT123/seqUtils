
#' Align sequences to a reference using MAFFT
#'
#' Aligns sequences to a reference using MAFFT's --addfragments and --keeplength
#' flags. Output sequences preserve the reference length. Duplicate sequences are
#' aligned once for efficiency. Results are uppercase with names preserved.
#'
#' @param unaligned_sequences Character vector of unaligned sequences (can include NAs)
#' @param reference_sequence Single reference sequence to align to
#'
#' @return Character vector of aligned sequences (same length as reference), with
#'   names preserved. NA sequences remain NA.
#'
#' @note Requires MAFFT to be installed and available in PATH
#' @export
mafft_align = function(unaligned_sequences, reference_sequence){

  NA_sequence_locations = is.na(unaligned_sequences)

  temp_mafft_folder = tempdir()
  tempfile()

  reference_seq_file = tempfile(
    pattern = "reference",
    fileext = ".fasta",
    tmpdir = temp_mafft_folder
  )
  unaligned_seqs_file = tempfile(
    pattern = "sequences",
    fileext = ".fasta",
    tmpdir = temp_mafft_folder
  )
  aligned_seqs_file = tempfile(
    pattern = "aligned",
    fileext = ".fasta",
    tmpdir = temp_mafft_folder
  )

  unique_unaligned_sequences = unique(unaligned_sequences)
  all_to_unique_sequence_map = match(unaligned_sequences, unique_unaligned_sequences)
  message('Number to align = ', length(unique_unaligned_sequences))

  write_fast_fasta(
    unique_unaligned_sequences,
    names = seq_along(unique_unaligned_sequences),
    path  = unaligned_seqs_file
  )

  write_fast_fasta(reference_sequence, names = 'reference', path = reference_seq_file)

  system(
    paste0(
      'mafft --thread 8 --quiet --6merpair --keeplength --addfragments ',
      unaligned_seqs_file,
      ' ',
      reference_seq_file, ' > ',
      aligned_seqs_file
    )
  )

  unique_aligned_sequences = toupper(fast_fasta(aligned_seqs_file)[-1])
  all_aligned_sequences = unique_aligned_sequences[all_to_unique_sequence_map]

  all_aligned_sequences[NA_sequence_locations] = NA


  names(all_aligned_sequences) = names(unaligned_sequences)

  all_aligned_sequences
}
