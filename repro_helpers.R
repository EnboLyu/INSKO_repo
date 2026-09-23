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

# Write a per-run provenance record (sessionInfo + loaded/installed versions).
dump_session_info <- function(outdir, tag = NULL) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  fn <- file.path(
    outdir,
    if (is.null(tag)) "sessionInfo.txt" else sprintf("sessionInfo_%s.txt", tag)
  )
  con <- file(fn, "w")
  on.exit(close(con))
  writeLines(c(
    "# INSKO run provenance",
    paste("# generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("# working_dir:", getwd()),
    paste("# libPaths:", paste(.libPaths(), collapse = " ; ")),
    ""
  ), con)
  writeLines(capture.output(sessionInfo()), con)
  invisible(fn)
}
