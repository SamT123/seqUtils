test_that("apply_substitutions works with single amino acid substitution", {
  result = apply_substitutions("QKLPFST", c("K2N"))
  expect_equal(result, "QNLPFST")

  result = apply_substitutions("ACDEFGH", c("D3X"))
  expect_equal(result, "ACXEFGH")
})

test_that("apply_substitutions works with multiple amino acid substitutions", {
  result = apply_substitutions("QKLPFST", c("K2N", "L3M", "T7A"))
  expect_equal(result, "QNMPFSA")

  result = apply_substitutions("AAAAA", c("A1G", "A3R", "A5Q"))
  expect_equal(result, "GARAQ")
})

test_that("apply_substitutions works with nucleotide substitutions", {
  result = apply_substitutions("ATCGATCG", c("T2C"), type = "nt")
  expect_equal(result, "ACCGATCG")

  result = apply_substitutions("AAGGCCTT", c("A1G", "G3A", "T8A"), type = "nt")
  expect_equal(result, "GAAGCCTA")
})

test_that("apply_substitutions handles gaps and special characters", {
  result = apply_substitutions("ACD-FGH", c("-4X"), type = "aa")
  expect_equal(result, "ACDXFGH")

  result = apply_substitutions("ACDEFGH", c("E4-"), type = "aa")
  expect_equal(result, "ACD-FGH")

  result = apply_substitutions("ATNNGCT", c("N3A", "N4T"), type = "nt")
  expect_equal(result, "ATATGCT")
})

test_that("apply_substitutions errors on mismatched from character", {
  expect_error(
    apply_substitutions("QKLPFST", c("A2N")),
    "from_char mismatch"
  )

  expect_error(
    apply_substitutions("ATCG", c("G2T"), type = "nt"),
    "from_char mismatch"
  )
})

test_that("apply_substitutions errors on invalid substitution format", {
  expect_error(
    apply_substitutions("QKLPFST", c("Z2N")),
    "Invalid substitution"
  )

  expect_error(
    apply_substitutions("QKLPFST", c("KXN")),
    "Invalid substitution"
  )

  expect_error(
    apply_substitutions("ATCG", c("K2N"), type = "nt"),
    "Invalid substitution"
  )
})

test_that("apply_substitutions errors on invalid type parameter", {
  expect_error(
    apply_substitutions("QKLPFST", c("K2N"), type = "protein"),
    "type must be 'aa' or 'nt'"
  )

  expect_error(
    apply_substitutions("ATCG", c("A1T"), type = "dna"),
    "type must be 'aa' or 'nt'"
  )
})
