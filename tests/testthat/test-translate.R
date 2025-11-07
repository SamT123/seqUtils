# Helper function to check if mafft is available
mafft_available = function() {
  result = suppressWarnings(system(
    "which mafft",
    ignore.stdout = TRUE,
    ignore.stderr = TRUE
  ))
  return(result == 0)
}

test_that("translate works for basic DNA sequence without deletions", {
  seq = c("ATGGCC")
  result = translate(seq)

  expect_equal(as.character(result[1]), "MA")
  expect_equal(length(result), 1)
})

test_that("translate preserves sequence names", {
  seqs = c("ATGGCC", "ATGTTT")
  names(seqs) = c("seq1", "seq2")

  result = translate(seqs)

  expect_equal(names(result), c("seq1", "seq2"))
  expect_equal(as.character(result["seq1"]), "MA")
  expect_equal(as.character(result["seq2"]), "MF")
})

test_that("translate errors when deletions present without reference", {
  seq = c("ATG---GCC")

  expect_error(
    translate(seq),
    "Deletions present but no reference for alignment provided"
  )
})

test_that("translate with reference_aas handles deletions correctly", {
  skip_if_not(mafft_available(), "mafft not installed")

  seqs = c("ATGGCCAAA", "ATG---AAA")
  names(seqs) = c("ref_seq", "del_seq")

  ref_aa = translate(seqs[1])

  result = suppressMessages(translate(seqs, reference_aas = ref_aa))

  expect_equal(names(result), c("ref_seq", "del_seq"))
  # Both should align to the reference
  expect_equal(as.character(result["ref_seq"]), "MAK")
})

test_that("translate handles leading deletions with reference", {
  skip_if_not(mafft_available(), "mafft not installed")

  seq = c("-ATGGCCAA")
  ref = c("GQ")

  result = suppressMessages(translate(seq, reference_aas = ref))

  expect_equal(as.character(result[1]), "GQ")
})

test_that("translate handles NA sequences with warning", {
  seqs = c("ATGGCC", NA, "ATGTTT")
  names(seqs) = c("seq1", "seq2", "seq3")

  expect_warning(
    result <- translate(seqs),
    "Some sequences are NA"
  )

  expect_true(is.na(result["seq2"]))
  expect_equal(as.character(result["seq1"]), "MA")
  expect_equal(as.character(result["seq3"]), "MF")
})

test_that("translate handles sequences not divisible by 3", {
  seq = c("ATGGCCT")

  result = translate(seq)

  # Should only translate first 6 bases: ATGGCC -> MA
  expect_equal(as.character(result[1]), "MA")
})

test_that("translate handles fuzzy codons", {
  seq = c("ATGNNN")

  result = translate(seq)

  # Should translate ATG to M, NNN to X (unknown)
  expect_equal(as.character(result[1]), "MX")
})

test_that("translate handles stop codons", {
  seq = c("ATGTAA")

  result = translate(seq)

  # Should translate ATG to M, TAA to * (stop)
  expect_equal(as.character(result[1]), "M*")
})

test_that("translate handles empty sequence vector", {
  seqs = character(0)

  result = translate(seqs)

  expect_equal(length(result), 0)
})
