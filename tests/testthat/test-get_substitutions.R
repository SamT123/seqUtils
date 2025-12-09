test_that("get_substitutions works with single and multiple differences", {
  # Single difference
  result1 = get_substitutions("ACDEF", "ACDEG")
  expect_equal(result1, "F5G")
  expect_type(result1, "character")

  # Multiple differences
  result2 = get_substitutions("ACDEF", "TCDXF")
  expect_equal(length(result2), 2)
  expect_true("A1T" %in% result2)
  expect_true("E4X" %in% result2)
})

test_that("get_substitutions returns list for multiple query sequences", {
  seq1 = "ACDEF"
  seq2 = c("ACDEG", "TCDEG", "ACDEF")

  result = get_substitutions(seq1, seq2)

  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_equal(result[[1]], "F5G")
  expect_equal(length(result[[2]]), 2)
  # Third is identical - should return character(0)
  expect_equal(length(result[[3]]), 0)
  expect_type(result[[3]], "character")
})

test_that("get_substitutions works with custom position_map", {
  seq1 = "ACE"
  seq2 = "TCE"
  position_map = c(145, 155, 156)

  result = get_substitutions(seq1, seq2, position_map = position_map)

  expect_equal(result, "A145T")
})

test_that("get_substitutions excludes specified characters", {
  # Exclude both positions with X
  result1 = get_substitutions("ACDEF", "XCDXF", exclude = "X")
  expect_equal(length(result1), 0)
  expect_type(result1, "character")

  # Exclude only affects positions with X
  result2 = get_substitutions("ACDXF", "ACDEG", exclude = "X")
  expect_equal(result2, "F5G")
})

test_that("get_substitutions handles edge cases", {
  # Identical sequences
  expect_equal(length(get_substitutions("ACDEF", "ACDEF")), 0)

  # Empty sequences
  result = get_substitutions("", "")
  expect_equal(length(result), 0)
  expect_type(result, "character")
})

test_that("get_substitutions trims to shortest length", {
  expect_message(
    result <- get_substitutions("ACDEF", "ACDEGXYZ"),
    "Trimming to shortest length: 5"
  )
  expect_equal(result, "F5G")
})

test_that("get_substitutions simplify parameter controls return type", {
  # simplify = TRUE (default) returns character vector
  result1 = get_substitutions("ACDEF", "ACDEG", simplify = TRUE)
  expect_false(is.list(result1))
  expect_type(result1, "character")

  # simplify = FALSE returns list
  result2 = get_substitutions("ACDEF", "ACDEG", simplify = FALSE)
  expect_type(result2, "list")
  expect_equal(result2[[1]], "F5G")
})

test_that("get_substitutions preserves names from sequence_2 (multiple sequences)", {
  seq1 = "ACDEF"
  seq2 = c(sample1 = "ACDEG", sample2 = "TCDEG", sample3 = "ACDEF")

  result = get_substitutions(seq1, seq2)

  expect_type(result, "list")
  expect_equal(names(result), c("sample1", "sample2", "sample3"))
  expect_equal(result[["sample1"]], "F5G")
  expect_equal(length(result[["sample3"]]), 0)
})

test_that("get_substitutions returns unnamed output when sequence_2 is unnamed", {
  seq1 = "ACDEF"
  seq2 = c("ACDEG", "TCDEG", "ACDEF")

  result = get_substitutions(seq1, seq2)

  expect_type(result, "list")
  expect_null(names(result))
  expect_equal(length(result), 3)
})

test_that("get_substitutions handles empty sequence_2 with allow_empty_sequence_2", {
  seq1 = "ACDEF"
  seq2 = character(0)

  # Should error by default
  expect_error(
    get_substitutions(seq1, seq2),
    "sequence_2 must contain at least one sequence"
  )

  # Should return empty list when allowed
  result = get_substitutions(seq1, seq2, allow_empty_sequence_2 = TRUE)
  expect_type(result, "list")
  expect_equal(length(result), 0)
})
