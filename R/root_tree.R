#' Root a tree using an outgroup sequence
#'
#' Adds an outgroup sequence to the tree using CMAPLE, roots via that outgroup,
#' and optionally removes the outgroup afterward.
#'
#' @param tree Phylogenetic tree to root
#' @param sequences Named vector of sequences for tree tips
#' @param outsequence Outgroup sequence (will be truncated to match sequence length)
#' @param cmaple_path Path to CMAPLE executable
#' @param remove_outsequence Remove outgroup from final tree (default TRUE)
#'
#' @return Rooted phylogenetic tree
#' @export
root_tree_using_outsequence = function(
  tree,
  sequences,
  outsequence,
  cmaple_path = NULL,
  remove_outsequence = TRUE
) {
  sequences = c(
    sequences,
    c(
      outsequence = substr(
        outsequence,
        1,
        nchar(sequences[1])
      )
    )
  )

  temp_tree_without_outsequence = fs::file_temp()
  temp_tree_with_outsequence = fs::file_temp()

  castor::write_tree(
    tree = tree,
    file = temp_tree_without_outsequence
  )

  make_cmaple_tree(
    sequences,
    tree_path = temp_tree_with_outsequence,
    starting_tree_path = temp_tree_without_outsequence,
    freeze_starting_tree = TRUE,
    cmaple_path = cmaple_path,
    multifurcating = !castor::is_bifurcating(tree)
  )

  tree_with_outsequence = castor::read_tree(
    file = temp_tree_with_outsequence
  )

  rooted_tree_with_outsequence = castor::root_via_outgroup(
    tree_with_outsequence,
    outgroup = "outsequence"
  )

  if (!remove_outsequence) {
    return(rooted_tree_with_outsequence)
  }

  ape::drop.tip(
    rooted_tree_with_outsequence,
    tip = "outsequence"
  )
}
