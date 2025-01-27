
#' Align sequences using mafft
#'
#' @param unaligned_sequences Character vector of unaligned sequences
#' @param reference_sequence Single sequence to align `unaligned_sequences` to
#'
#' @return Character vector of aligned sequences, with names preserved
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
