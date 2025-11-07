test_that("clean_sequences removes non-standard nucleotides", {
  seq = c("ATCGXYZ")
  result = clean_sequences(seq, type = "nt")

  # X, Y, Z should be replaced with N
  expect_equal(as.character(result[1]), "ATCGNNN")
})

test_that("clean_sequences removes non-standard amino acids", {
  seq = c("ACDEFG123")
  result = clean_sequences(seq, type = "aa")

  # 1, 2, 3 should be replaced with X
  expect_equal(as.character(result[1]), "ACDEFGXXX")
})

test_that("clean_sequences preserves sequence names", {
  seqs = c("ATCGXYZ", "GGGNNN")
  names(seqs) = c("seq1", "seq2")

  result = clean_sequences(seqs, type = "nt")

  expect_equal(names(result), c("seq1", "seq2"))
})

test_that("clean_sequences converts to uppercase", {
  seq = c("atcgxyz")
  result = clean_sequences(seq, type = "nt")

  # Should be uppercase and X, Y, Z replaced with N
  expect_equal(as.character(result[1]), "ATCGNNN")
})

test_that("clean_sequences handles gaps in nucleotides", {
  seq = c("ATC-G")
  result = clean_sequences(seq, type = "nt")

  # Gap should be preserved
  expect_equal(as.character(result[1]), "ATC-G")
})

test_that("clean_sequences handles gaps in amino acids", {
  seq = c("ACG-EF")
  result = clean_sequences(seq, type = "aa")

  # Gap should be preserved
  expect_equal(as.character(result[1]), "ACG-EF")
})

test_that("clean_sequences with custom replacement character for nucleotides", {
  seq = c("ATCGXYZ")
  result = clean_sequences(seq, type = "nt", replacement_character = "-")

  # X, Y, Z should be replaced with -
  expect_equal(as.character(result[1]), "ATCG---")
})

test_that("clean_sequences with custom replacement character for amino acids", {
  seq = c("ACDEF123")
  result = clean_sequences(seq, type = "aa", replacement_character = "-")

  # 1, 2, 3 should be replaced with -
  expect_equal(as.character(result[1]), "ACDEF---")
})

test_that("clean_sequences with custom alphabet", {
  seq = c("ATCGU")
  # Custom alphabet that includes U (for RNA)
  result = clean_sequences(
    seq,
    type = "nt",
    alphabet = c("A", "U", "C", "G", "-")
  )

  # T should be replaced with N (not in custom alphabet)
  expect_equal(as.character(result[1]), "ANCGU")
})

test_that("clean_sequences handles all standard amino acids", {
  # All 20 standard amino acids
  seq = c("ARNDCQEGHILKMFPSTWYV")
  result = clean_sequences(seq, type = "aa")

  # Should remain unchanged
  expect_equal(as.character(result[1]), "ARNDCQEGHILKMFPSTWYV")
})

test_that("clean_sequences handles all standard nucleotides", {
  seq = c("ATCG")
  result = clean_sequences(seq, type = "nt")

  # Should remain unchanged
  expect_equal(as.character(result[1]), "ATCG")
})

test_that("clean_sequences handles multiple sequences", {
  seqs = c("ATCGXYZ", "GGGNNN", "AAA123")
  names(seqs) = c("seq1", "seq2", "seq3")

  result = clean_sequences(seqs, type = "nt")

  expect_equal(as.character(result["seq1"]), "ATCGNNN")
  expect_equal(as.character(result["seq2"]), "GGGNNN")
  expect_equal(as.character(result["seq3"]), "AAANNN")
  expect_equal(names(result), c("seq1", "seq2", "seq3"))
})

test_that("clean_sequences handles empty sequence", {
  seq = c("")
  result = clean_sequences(seq, type = "nt")

  expect_equal(as.character(result[1]), "")
})

test_that("clean_sequences handles sequence with only invalid characters", {
  seq = c("XYZ123")
  result = clean_sequences(seq, type = "nt")

  # All should be replaced with N
  expect_equal(as.character(result[1]), "NNNNNN")
})

test_that("clean_sequences handles special characters", {
  seq = c("ATC*G@T#")
  result = clean_sequences(seq, type = "nt")

  # *, @, # should be replaced with N
  expect_equal(as.character(result[1]), "ATCNGNTN")
})

test_that("clean_sequences handles ambiguous nucleotides", {
  # R, Y, W, S are ambiguous nucleotide codes (not in default alphabet)
  seq = c("ATCGRYWSN")
  result = clean_sequences(seq, type = "nt")

  # R, Y, W, S, N should be replaced with N (only ATCG- are standard)
  expect_equal(as.character(result[1]), "ATCGNNNNN")
})

test_that("clean_sequences with custom alphabet and replacement", {
  seq = c("ATCGUUGG")
  result = clean_sequences(
    seq,
    type = "nt",
    alphabet = c("A", "U", "C", "G"),
    replacement_character = "X"
  )

  # T should be replaced with X
  expect_equal(as.character(result[1]), "AXCGUUGG")
})

test_that("clean_sequences handles lowercase mixed case", {
  seq = c("AtCgXyZ")
  result = clean_sequences(seq, type = "nt")

  expect_equal(as.character(result[1]), "ATCGNNN")
})

test_that("clean_sequences handles spaces and tabs", {
  seq = c("ATC G\tT")
  result = clean_sequences(seq, type = "nt")

  # Spaces and tabs should be replaced with N
  expect_equal(as.character(result[1]), "ATCNGNT")
})

test_that("clean_sequences with unnamed sequences", {
  seqs = c("ATCGXYZ", "GGGNNN")

  result = clean_sequences(seqs, type = "nt")

  expect_equal(length(result), 2)
  expect_equal(as.character(result[1]), "ATCGNNN")
  expect_equal(as.character(result[2]), "GGGNNN")
})
