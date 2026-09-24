# Theil-Sen slope ----------

theilSenSlope <- function(x, y, n_points = 5000, seed = 1) {
  idx <- withr::with_seed(
    seed,
    sample.int(length(x), min(length(x), n_points))
  )

  numerator <- outer(y[idx], y[idx], "-")
  denominator <- outer(x[idx], x[idx], "-")

  stats::median((numerator / denominator)[denominator != 0])
}

# Date outliers ----------

#' Flag tips whose collection dates conflict with the phylogeny
#'
#' Flags a tip whose root-to-tip divergence sits more than `iqd` interquartile
#' ranges from a Theil-Sen root-to-tip regression line. The slope is robust so
#' that the outliers being hunted do not distort the line they are tested
#' against, and is returned alongside the flags so a caller can reuse it as a
#' clock rate.
#'
#' The test is invariant to the scale of the branch lengths, so the tree may be
#' in substitutions per site or in substitutions. Fewer than 20 dated tips
#' flags nothing, matching `chronumental`'s own guard.
#'
#' @param tree A rooted `phylo` object.
#' @param dates A `Date` vector, one per tip, in `tree$tip.label` order.
#' @param iqd Interquartile ranges from the regression line beyond which a tip
#'   is flagged.
#' @param fit Optional logical vector, one per tip: the tips that set the
#'   line and IQR (e.g. an even subsample of an unevenly sampled tree). All
#'   dated tips are still tested.
#'
#' @return A list with `outliers`, a tibble of `label`, `collection_date` and
#'   `root_to_tip_residual_days` for the flagged tips, and `slope_per_day`, the
#'   Theil-Sen slope in tree units per day.
#'
#' @export
find_date_outliers <- function(tree, dates, iqd = 3, fit = NULL) {
  if (length(dates) != ape::Ntip(tree)) {
    rlang::abort(
      "`dates` must have one element per tip.",
      class = "seqUtils_error_date_length"
    )
  }
  if (!is.null(fit) && (!is.logical(fit) || length(fit) != ape::Ntip(tree))) {
    rlang::abort(
      "`fit` must be a logical vector with one element per tip.",
      class = "seqUtils_error_fit_length"
    )
  }

  days <- as.numeric(dates)
  root_to_tip <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
  dated <- !is.na(days)
  fitting <- dated & (if (is.null(fit)) TRUE else fit)

  if (sum(fitting) < 20) {
    return(list(
      outliers = emptyOutlierTable(),
      slope_per_day = NA_real_
    ))
  }

  slope_per_day <- theilSenSlope(days[fitting], root_to_tip[fitting])

  residuals <- root_to_tip - slope_per_day * days
  residuals <- residuals - stats::median(residuals[fitting])
  cutoff <- iqd * diff(stats::quantile(residuals[fitting], c(0.25, 0.75)))

  flagged <- which(dated & abs(residuals) > cutoff)

  list(
    outliers = tibble::tibble(
      label = tree$tip.label[flagged],
      collection_date = dates[flagged],
      root_to_tip_residual_days = residuals[flagged] / slope_per_day
    ),
    slope_per_day = slope_per_day
  )
}

emptyOutlierTable <- function() {
  tibble::tibble(
    label = character(),
    collection_date = as.Date(character()),
    root_to_tip_residual_days = numeric()
  )
}
