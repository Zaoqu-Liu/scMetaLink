#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("scMetaLink v", utils::packageVersion("scMetaLink"), " loaded")
}

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.scMetaLink <- list(
    scMetaLink.verbose = TRUE,
    scMetaLink.n_cores = 1
  )
  toset <- !(names(op.scMetaLink) %in% names(op))
  if (any(toset)) options(op.scMetaLink[toset])
  invisible()
}

#' Show Package Citation
#' @export
citationScMetaLink <- function() {
  cat("
scMetaLink: Inferring metabolite-mediated cell-cell communication
from single-cell transcriptomic data.

Author: Zaoqu Liu
")
}
