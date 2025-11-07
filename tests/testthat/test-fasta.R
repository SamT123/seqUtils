test_that("fast_fasta reads single-line FASTA format", {
  # Create a temporary FASTA file with single-line sequences
  tmp_file = tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1",
      "ACGTACGT",
      ">seq2",
      "TGCATGCA"
    ),
    tmp_file
  )

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(length(result), 2)
  expect_equal(names(result), c("seq1", "seq2"))
  expect_equal(as.character(result[1]), "ACGTACGT")
  expect_equal(as.character(result[2]), "TGCATGCA")

  unlink(tmp_file)
})

test_that("fast_fasta reads multi-line FASTA format", {
  # Create a temporary FASTA file with multi-line sequences
  tmp_file = tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1",
      "ACGTACGT",
      "GGGGCCCC",
      ">seq2",
      "TGCATGCA",
      "AAAATTTT",
      "CCCCGGGG"
    ),
    tmp_file
  )

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(length(result), 2)
  expect_equal(names(result), c("seq1", "seq2"))
  expect_equal(as.character(result[1]), "ACGTACGTGGGGCCCC")
  expect_equal(as.character(result[2]), "TGCATGCAAAAATTTTCCCCGGGG")

  unlink(tmp_file)
})

test_that("fast_fasta handles mixed single-line and multi-line sequences", {
  tmp_file = tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1",
      "ACGTACGT",
      ">seq2",
      "TGCATGCA",
      "AAAATTTT"
    ),
    tmp_file
  )

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(length(result), 2)
  expect_equal(as.character(result[1]), "ACGTACGT")
  expect_equal(as.character(result[2]), "TGCATGCAAAAATTTT")

  unlink(tmp_file)
})

test_that("fast_fasta handles sequence names with spaces and special characters", {
  tmp_file = tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1 description with spaces",
      "ACGT",
      ">seq2|accession|description",
      "TGCA"
    ),
    tmp_file
  )

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(
    names(result),
    c("seq1 description with spaces", "seq2|accession|description")
  )

  unlink(tmp_file)
})

test_that("fast_fasta handles empty lines in sequences", {
  tmp_file = tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1",
      "ACGT",
      "",
      "TGCA"
    ),
    tmp_file
  )

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(length(result), 1)
  expect_equal(as.character(result[1]), "ACGTTGCA")

  unlink(tmp_file)
})

test_that("fast_fasta prints correct message about number of sequences", {
  tmp_file = tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1",
      "ACGT",
      ">seq2",
      "TGCA",
      ">seq3",
      "GGGG"
    ),
    tmp_file
  )

  expect_message(fast_fasta(tmp_file), "Loaded 3 sequences")

  unlink(tmp_file)
})

test_that("fast_fasta handles single sequence", {
  tmp_file = tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">single_seq",
      "ACGTACGTACGTACGT"
    ),
    tmp_file
  )

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(length(result), 1)
  expect_equal(names(result), "single_seq")
  expect_equal(as.character(result[1]), "ACGTACGTACGTACGT")

  unlink(tmp_file)
})

test_that("fast_fasta handles very long multi-line sequences", {
  tmp_file = tempfile(fileext = ".fasta")
  # Create a sequence split across many lines
  seq_lines = c(">long_seq")
  for (i in 1:10) {
    seq_lines = c(seq_lines, strrep("ACGT", 20))
  }
  writeLines(seq_lines, tmp_file)

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(length(result), 1)
  expect_equal(unname(nchar(result[1])), 80 * 10) # 80 chars * 10 lines

  unlink(tmp_file)
})

test_that("write_fast_fasta creates valid single-line FASTA", {
  seqs = c("ACGTACGT", "TGCATGCA")
  names(seqs) = c("seq1", "seq2")

  tmp_file = tempfile(fileext = ".fasta")
  write_fast_fasta(seqs, path = tmp_file)

  # Read it back
  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(result, seqs)

  unlink(tmp_file)
})

test_that("write_fast_fasta with line_width wraps sequences", {
  seqs = c("ACGTACGTACGTACGTACGTACGT")
  names(seqs) = c("seq1")

  tmp_file = tempfile(fileext = ".fasta")
  write_fast_fasta(seqs, path = tmp_file, line_width = 10)

  lines = readLines(tmp_file)

  expect_equal(lines[1], ">seq1")
  expect_equal(lines[2], "ACGTACGTAC")
  expect_equal(lines[3], "GTACGTACGT")
  expect_equal(lines[4], "ACGT")

  # Read it back and verify sequence is correct
  result = suppressMessages(fast_fasta(tmp_file))
  expect_equal(as.character(result[1]), seqs[[1]])

  unlink(tmp_file)
})

test_that("write_fast_fasta with line_width handles short sequences", {
  seqs = c("ACGT")
  names(seqs) = c("seq1")

  tmp_file = tempfile(fileext = ".fasta")
  write_fast_fasta(seqs, path = tmp_file, line_width = 10)

  lines = readLines(tmp_file)

  expect_equal(lines[1], ">seq1")
  expect_equal(lines[2], "ACGT")
  expect_equal(length(lines), 2)

  unlink(tmp_file)
})

test_that("write_fast_fasta with unnamed sequences and names argument", {
  seqs = c("ACGT", "TGCA")
  seq_names = c("seq1", "seq2")

  tmp_file = tempfile(fileext = ".fasta")
  write_fast_fasta(seqs, names = seq_names, path = tmp_file)

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(names(result), seq_names)
  expect_equal(as.character(result), seqs)

  unlink(tmp_file)
})

test_that("write_fast_fasta round-trip preserves data", {
  seqs = c("ACGTACGT", "TGCATGCATGCA", "GGGGCCCCAAAATTTT")
  names(seqs) = c("sequence_1", "sequence_2", "sequence_3")

  tmp_file = tempfile(fileext = ".fasta")
  write_fast_fasta(seqs, path = tmp_file)

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(result, seqs)

  unlink(tmp_file)
})

test_that("write_fast_fasta with line_width round-trip preserves data", {
  seqs = c("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA")
  names(seqs) = c("seq1", "seq2")

  tmp_file = tempfile(fileext = ".fasta")
  write_fast_fasta(seqs, path = tmp_file, line_width = 8)

  result = suppressMessages(fast_fasta(tmp_file))

  expect_equal(result, seqs)

  unlink(tmp_file)
})
