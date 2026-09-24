require_packages <- function(pkgs, context = "this script") {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      sprintf(
        paste0(
          "Missing package(s) for %s: %s\n",
          "Restore the pinned environment instead of installing ad hoc:\n",
          "  R_LIBS_USER=.Rlib Rscript -e 'renv::restore()'\n",
          "(renv.lock is committed at the repo root.)"
        ),
        context, paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
