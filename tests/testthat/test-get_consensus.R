test_that("get_consensus returns simple consensus from identical sequences", {
  sequences = c("ACGT", "ACGT", "ACGT")
  result = get_consensus(sequences)
  expect_equal(result, "ACGT")
})

test_that("get_consensus returns consensus from mixed sequences", {
  sequences = c("ACGT", "ACGT", "ACAT")
  result = get_consensus(sequences)
  # Position 3: G appears 2/3 times (67%), T appears 1/3 (33%)
  expect_equal(result, "ACGT")
})

test_that("get_consensus uses min_freq threshold correctly", {
  sequences = c("ACGT", "ACAT", "ATTT")
  # Position 2: C appears 2/3 (67%), position 3: G,A,T each appear 1/3 (33%)
  result = get_consensus(sequences, min_freq = 0.5)
  expect_equal(result, "AC?T")

  result = get_consensus(sequences, min_freq = 0.7)
  expect_equal(result, "A??T")
})

test_that("get_consensus excludes specified characters", {
  sequences = c("AC-T", "ACGT", "AC-T")
  # Without excluding "-": position 3 has "-" at 2/3 (67%)
  result = get_consensus(sequences)
  expect_equal(result, "AC-T")

  # With excluding "-": position 3 has only "G" at 1/1 (100%)
  result = get_consensus(sequences, excluded_characters = c("-"))
  expect_equal(result, "ACGT")
})

test_that("get_consensus handles multiple excluded characters", {
  sequences = c("ACX-", "ACGT", "ACNX")
  result = get_consensus(sequences, excluded_characters = c("-", "X", "N"))
  # After excluding X, -, N: pos 3 only has G, pos 4 only has T
  expect_equal(result, "ACGT")
})

test_that("get_consensus with no valid characters", {
  sequences = c("AC-T", "AC-T")
  result = get_consensus(
    sequences,
    excluded_characters = c("-"),
    min_freq = 0
  )
  expect_equal(result, "AC?T")
})


test_that("get_consensus with mixed case returns uppercase", {
  sequences = c("acgt", "ACGT", "AcGt")
  result = get_consensus(sequences)
  # Biostrings should handle this and return uppercase
  expect_equal(result, "ACGT")
})

test_that("get_consensus handles X in amino acids correctly", {
  sequences = c("QXLPFST", "QKLPFST", "QXLPFST")
  # Without excluding X, it appears 2/3 times
  result = get_consensus(sequences)
  expect_equal(result, "QXLPFST")

  # With excluding X, K is the only option
  result = get_consensus(sequences, excluded_characters = c("X"))
  expect_equal(result, "QKLPFST")
})
