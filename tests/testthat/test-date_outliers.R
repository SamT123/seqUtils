makeClocklikeTree <- function(n_tips = 60, rate = 0.01, seed = 1) {
  withr::with_seed(seed, {
    days <- sample.int(2000, n_tips, replace = TRUE)
    jitter <- stats::rnorm(n_tips, sd = rate * 5)
  })

  tree <- ape::stree(n_tips, type = "star")
  tree$tip.label <- paste0("t", seq_len(n_tips))
  tree$edge.length <- pmax(rate * days + jitter, 0)

  list(
    tree = tree,
    dates = as.Date("2015-01-01") + days
  )
}

# 60 sparse old tips above the line set by 3,000 dense recent ones; t10 and
# t30 are misdated.
makeUnevenTree <- function(rate = 0.01, seed = 1) {
  withr::with_seed(seed, {
    days <- c(sample(0:2999, 60), sample(3000:3650, 3000, replace = TRUE))
    jitter <- stats::rnorm(length(days), sd = rate * 20)
  })
  divergence <- rate * days + rate * 0.4 * pmax(3000 - days, 0) + jitter
  divergence[c(10, 30)] <- divergence[c(10, 30)] + rate * 2000

  tree <- ape::stree(length(days), type = "star")
  tree$tip.label <- paste0("t", seq_along(days))
  tree$edge.length <- divergence

  list(
    tree = tree,
    dates = as.Date("2015-01-01") + days,
    even = seq_along(days) <= 60 |
      seq_along(days) %in% (60 + seq(1, 3000, by = 50))
  )
}

test_that("a tip far off the clock line is flagged", {
  sim <- makeClocklikeTree()
  sim$tree$edge.length[[7]] <- sim$tree$edge.length[[7]] + 30

  result <- find_date_outliers(sim$tree, sim$dates)

  expect_equal(result$outliers$label, "t7")
})

test_that("a clocklike tree flags nothing", {
  sim <- makeClocklikeTree()

  result <- find_date_outliers(sim$tree, sim$dates)

  expect_equal(nrow(result$outliers), 0L)
})

test_that("the reported residual is signed and in days", {
  sim <- makeClocklikeTree(rate = 0.01)
  sim$tree$edge.length[[7]] <- sim$tree$edge.length[[7]] + 0.01 * 400

  result <- find_date_outliers(sim$tree, sim$dates)

  expect_equal(result$outliers$root_to_tip_residual_days, 400, tolerance = 0.1)
})

test_that("fewer than 20 dated tips flags nothing", {
  sim <- makeClocklikeTree(n_tips = 19)

  result <- find_date_outliers(sim$tree, sim$dates)

  expect_equal(nrow(result$outliers), 0L)
  expect_true(is.na(result$slope_per_day))
})

test_that("all-NA dates returns an empty table rather than erroring", {
  sim <- makeClocklikeTree()

  result <- find_date_outliers(sim$tree, as.Date(rep(NA, 60)))

  expect_equal(nrow(result$outliers), 0L)
})

test_that("the flagged set is invariant to the branch length scale", {
  sim <- makeClocklikeTree()
  sim$tree$edge.length[[7]] <- sim$tree$edge.length[[7]] + 30
  scaled <- sim$tree
  scaled$edge.length <- scaled$edge.length * 1650

  expect_equal(
    find_date_outliers(sim$tree, sim$dates)$outliers$label,
    find_date_outliers(scaled, sim$dates)$outliers$label
  )
})

test_that("a dates vector of the wrong length is an error", {
  sim <- makeClocklikeTree()

  expect_error(
    find_date_outliers(sim$tree, sim$dates[-1]),
    class = "seqUtils_error_date_length"
  )
})

test_that("fitting on every tip of an uneven tree cuts the sparse old epoch", {
  sim <- makeUnevenTree()

  flagged <- find_date_outliers(sim$tree, sim$dates)$outliers$label

  expect_gt(sum(flagged %in% paste0("t", 1:60)), 40)
})

test_that("fitting on an even subsample flags only the misdated tips", {
  sim <- makeUnevenTree()

  result <- find_date_outliers(sim$tree, sim$dates, fit = sim$even)

  expect_setequal(result$outliers$label, c("t10", "t30"))
})

test_that("a fit vector of the wrong length is an error", {
  sim <- makeClocklikeTree()

  expect_error(
    find_date_outliers(sim$tree, sim$dates, fit = rep(TRUE, 59)),
    class = "seqUtils_error_fit_length"
  )
})

test_that("theilSenSlope resists contamination that shifts least squares", {
  # Mis-dated tips are reported earlier than they are, so they sit at low x
  # with undiminished divergence and drag a least-squares slope down.
  withr::with_seed(1, {
    x <- stats::runif(500, 0, 1000)
    y <- 0.01 * x + stats::rnorm(500, sd = 0.05)
  })
  contaminated <- order(x)[seq_len(25)]
  y[contaminated] <- y[contaminated] + 6

  least_squares <- unname(coef(stats::lm(y ~ x))[[2]])
  theil_sen <- theilSenSlope(x, y)

  expect_equal(theil_sen, 0.01, tolerance = 0.05)
  expect_lt(abs(theil_sen - 0.01), abs(least_squares - 0.01))
})

test_that("theilSenSlope is deterministic across calls", {
  withr::with_seed(2, {
    x <- stats::runif(8000, 0, 1000)
    y <- 0.01 * x + stats::rnorm(8000, sd = 0.05)
  })

  expect_identical(theilSenSlope(x, y), theilSenSlope(x, y))
})
