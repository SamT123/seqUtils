#' Read FASTA files efficiently
#'
#' Reads FASTA format files, handling both single-line and multi-line sequences.
#' Sequence names have the leading ">" removed. Prints a message with the number
#' of sequences loaded.
#'
#' @param path Path to FASTA file
#'
#' @return Named character vector where names are sequence identifiers and values
#'   are the sequences (multi-line sequences are concatenated into single strings)
#'
#' @importFrom readr read_lines
#' @export
fast_fasta = function(path) {
  lines = readr::read_lines(path, skip_empty_rows = F)
  name_lines = which(substr(lines, 1, 1) == '>')
  message("Loaded ", length(name_lines), " sequences")

  if (
    length(name_lines) == length(lines) / 2 &&
      all(name_lines == seq(1, length(lines), 2))
  ) {
    # case: single-line sequences
    sequences = lines[-name_lines]
    names(sequences) = lines[name_lines]
  } else {
    # case: multi-line sequences
    sequences = rep(NA, length(name_lines))
    names(sequences) = rep(NA, length(name_lines))

    name_lines = c(name_lines, length(lines) + 1)

    for (i in seq_along(name_lines[-1])) {
      names(sequences)[[i]] = lines[[name_lines[[i]]]]
      sequences[[i]] = paste(
        lines[(name_lines[[i]] + 1):(name_lines[[i + 1]] - 1)],
        collapse = ""
      )
    }
  }

  names(sequences) = substring(names(sequences), 2)
  sequences
}


#' Write sequences to FASTA file
#'
#' @param seqs Character vector of sequences (named or unnamed)
#' @param names Sequence names (required if `seqs` is unnamed)
#' @param path Output file path
#' @param line_width Maximum characters per line. If NULL (default), sequences are
#'   written on single lines. Use 60 or 80 for standard wrapped format.
#'
#' @return NULL (invisibly)
#'
#' @importFrom readr write_lines
#' @export
write_fast_fasta = function(seqs, names = NULL, path, line_width = NULL) {
  if (is.null(names)) {
    names = names(seqs)
  }

  stopifnot(length(seqs) == length(names))

  if (!is.null(line_width)) {
    # Split sequences into chunks of line_width characters
    seqs = vapply(
      seqs,
      function(seq) {
        if (nchar(seq) <= line_width) {
          return(seq)
        }
        seq_chunks = character(ceiling(nchar(seq) / line_width))
        for (i in seq_along(seq_chunks)) {
          start_pos = (i - 1) * line_width + 1
          end_pos = min(i * line_width, nchar(seq))
          seq_chunks[i] = substr(seq, start_pos, end_pos)
        }
        paste(seq_chunks, collapse = "\n")
      },
      character(1),
      USE.NAMES = FALSE
    )
  }

  x <- vector(class(seqs), 2 * length(seqs))
  x[c(TRUE, FALSE)] <- paste0('>', names)
  x[c(FALSE, TRUE)] <- seqs

  readr::write_lines(x, file = path)

  invisible(NULL)
}
