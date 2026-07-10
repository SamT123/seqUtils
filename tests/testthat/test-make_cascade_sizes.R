# Build a minimal alignment with a given per-year count profile. Only
# Collection_date is read by make_cascade_sizes; sequences are placeholders.
make_year_alignment <- function(counts_by_year) {
  years <- rep(
    as.integer(names(counts_by_year)),
    times = unname(counts_by_year)
  )
  data.frame(
    Isolate_unique_identifier = paste0("s", seq_along(years)),
    dna_sequence = "ACGT",
    Collection_date = as.Date(paste0(years, "-06-15"))
  )
}


test_that("make_cascade_sizes doubles until the densest year is covered", {
  aln <- make_year_alignment(c("2000" = 8, "2001" = 3)) # max_bin = 8
  expect_equal(make_cascade_sizes(aln, 1), c(2L, 4L))
})


test_that("make_cascade_sizes is empty when 2*start covers the densest year", {
  aln <- make_year_alignment(c("2000" = 8, "2001" = 3)) # max_bin = 8
  expect_equal(make_cascade_sizes(aln, 4), integer(0))
  expect_equal(make_cascade_sizes(aln, 20), integer(0))
})


test_that("make_cascade_sizes sizes are strictly increasing and above start", {
  aln <- make_year_alignment(c("2010" = 300)) # single year, max_bin = 300
  sizes <- make_cascade_sizes(aln, 10)
  expect_equal(sizes, c(20L, 40L, 80L, 160L))
  expect_true(all(diff(sizes) > 0))
  expect_true(all(sizes > 10))
  expect_true(all(sizes < 300))
})
