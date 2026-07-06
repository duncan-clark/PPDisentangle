#' Find the PPDisentangle repository root
#'
#' Walks up from \code{start_dir} until a directory contains \code{DESCRIPTION}
#' and \code{NAMESPACE}, or a \code{.git} directory.
#'
#' @param start_dir Directory to start from (default: current working directory).
#' @return Normalized path to the repository root.
#' @export
find_repo_root <- function(start_dir = getwd()) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(cur, ".git"))) return(cur)
    if (file.exists(file.path(cur, "DESCRIPTION")) &&
        file.exists(file.path(cur, "NAMESPACE"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = TRUE)
}

#' Root directory for generated analysis outputs
#'
#' Outputs live **outside** the git repository by default, as a sibling folder
#' \code{PPDisentangle-output/} next to the package checkout. Override with the
#' environment variable \code{PPDISENTANGLE_OUTPUT_ROOT}.
#'
#' @param repo_root Package repository root; defaults to \code{find_repo_root()}.
#' @return Normalized path to the external output root.
#' @export
pp_output_root <- function(repo_root = NULL) {
  if (is.null(repo_root)) repo_root <- find_repo_root()
  env <- trimws(Sys.getenv("PPDISENTANGLE_OUTPUT_ROOT", ""))
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  normalizePath(
    file.path(dirname(repo_root), "PPDisentangle-output"),
    winslash = "/",
    mustWork = FALSE
  )
}

#' Path under the external output root
#'
#' @param ... Path components appended to \code{pp_output_root()}.
#' @param repo_root Package repository root; defaults to \code{find_repo_root()}.
#' @return Normalized path.
#' @export
pp_output_path <- function(..., repo_root = NULL) {
  file.path(pp_output_root(repo_root), ...)
}
