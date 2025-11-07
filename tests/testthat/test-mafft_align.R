# Helper function to check if mafft is available
mafft_available = function() {
  result = suppressWarnings(system(
    "which mafft",
    ignore.stdout = TRUE,
    ignore.stderr = TRUE
  ))
  return(result == 0)
}

test_that("mafft_align aligns sequences to reference", {
  skip_if_not(mafft_available(), "mafft not installed")

  reference = "ATCGATCGATCG"
  unaligned = c("ATCGATCG", "ATCGAT", "TCGATCGATCG")

  result = mafft_align(unaligned, reference)

  expect_length(result, 3)
  expect_true(all(nchar(result) == nchar(reference)))
})

test_that("mafft_align preserves sequence names", {
  skip_if_not(mafft_available(), "mafft not installed")

  reference = "ATCGATCGATCG"
  unaligned = c(seq1 = "ATCGATCG", seq2 = "ATCGAT", seq3 = "TCGATCGATCG")

  result = mafft_align(unaligned, reference)

  expect_equal(names(result), c("seq1", "seq2", "seq3"))
})

test_that("mafft_align handles NA sequences", {
  skip_if_not(mafft_available(), "mafft not installed")

  reference = "ATCGATCGATCG"
  unaligned = c("ATCGATCG", NA, "TCGATCGATCG", NA)

  result = mafft_align(unaligned, reference)

  expect_length(result, 4)
  expect_true(is.na(result[2]))
  expect_true(is.na(result[4]))
  expect_false(is.na(result[1]))
  expect_false(is.na(result[3]))
})

test_that("mafft_align handles duplicate sequences accurately", {
  skip_if_not(mafft_available(), "mafft not installed")

  reference = "ATCGATCGATCG"
  unaligned = c("ATCGATCG", "ATCGATCG", "ATCGAT", "ATCGATCG")

  result = mafft_align(unaligned, reference)

  expect_length(result, 4)
  # First, second, and fourth should be identical (same input)
  expect_equal(result[1], result[2])
  expect_equal(result[1], result[4])
})

test_that("mafft_align preserves names with NA sequences", {
  skip_if_not(mafft_available(), "mafft not installed")

  reference = "ATCGATCGATCG"
  unaligned = c(seq1 = "ATCGATCG", seq2 = NA, seq3 = "TCGATCGATCG")

  result = mafft_align(unaligned, reference)

  expect_equal(names(result), c("seq1", "seq2", "seq3"))
  expect_true(is.na(result["seq2"]))
})

test_that("mafft_align returns uppercase sequences", {
  skip_if_not(mafft_available(), "mafft not installed")

  reference = "ATCGATCGATCG"
  unaligned = c("atcgatcg", "ATCGAT")

  result = mafft_align(unaligned, reference)

  expect_true(all(result == toupper(result), na.rm = TRUE))
})

test_that("mafft_align handles single sequence", {
  skip_if_not(mafft_available(), "mafft not installed")

  reference = "ATCGATCGATCG"
  unaligned = c("ATCGATCG")

  result = mafft_align(unaligned, reference)

  expect_length(result, 1)
  expect_true(nchar(result) == nchar(reference))
})

test_that("mafft_align handles all NA sequences", {
  skip_if_not(mafft_available(), "mafft not installed")

  reference = "ATCGATCGATCG"
  unaligned = c(NA, NA, NA)

  result = mafft_align(unaligned, reference)

  expect_length(result, 3)
  expect_true(all(is.na(result)))
})

test_that("mafft_align errors when MAFFT is not installed", {
  # Temporarily clear PATH to simulate MAFFT not being installed
  old_path = Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path))
  Sys.setenv(PATH = "")

  reference = "ATCGATCGATCG"
  unaligned = c("ATCGATCG", "ATCGAT")

  expect_error(
    mafft_align(unaligned, reference),
    "MAFFT is not installed or not available in PATH"
  )
})
