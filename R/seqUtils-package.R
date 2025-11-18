#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom purrr map_chr
#' @importFrom stats setNames
## usethis namespace: end
NULL

.onAttach <- function(libname, pkgname) {
  if (system("mafft --version", ignore.stdout = T, ignore.stderr = T) != 0) {
    packageStartupMessage(
      "Note: MAFFT not found in PATH. Install MAFFT for alignment functions."
    )
  }
}
