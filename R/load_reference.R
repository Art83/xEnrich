# Package-level constant. change once for production
.ZENODO_RECORD_ID <- 19106109L


#' Load reference datasets from Zenodo cache
#'
#' Downloads one or more reference datasets from the package Zenodo record
#' and caches them locally. Subsequent calls with the same \code{datasets}
#' and \code{dest_dir} return the cached files without re-downloading unless
#' \code{overwrite = TRUE}.
#'
#' @param datasets Character vector of dataset names to load. If \code{NULL}
#'   or empty, all available datasets for \code{source} are downloaded.
#'   Use the manifest accessor functions to list available names.
#' @param source Character. \code{"tabula"} or \code{"hpa"}.
#' @param dest_dir Character. Local cache directory. Defaults to the
#'   platform-appropriate user cache via [tools::R_user_dir()].
#' @param overwrite Logical. Re-download even if a local copy exists.
#'   Default \code{FALSE}.
#' @param verify Logical. Verify SHA256 checksum after download.
#'   Default \code{TRUE}. Disable only for debugging.Corrupt downloads
#'   will silently propagate if verification is off.
#'
#' @return A named list of loaded R objects, one per requested dataset.
#'
#' @export
load_reference <- function(datasets  = NULL,
                           source    = c("tabula", "hpa"),
                           dest_dir  = xEnrich_cache_dir(),
                           overwrite = FALSE,
                           verify    = TRUE) {
  source   <- match.arg(source)
  manifest <- if (source == "tabula") tabula_manifest else hpa_manifest

  # Validate manifest structure upfront, catches broken sysdata early
  required_cols <- c("dataset", "filename", "sha256")
  missing_cols  <- setdiff(required_cols, names(manifest))
  if (length(missing_cols) > 0)
    stop("Manifest is missing required columns: ",
         paste(missing_cols, collapse = ", "), ".")

  available <- unique(manifest$dataset)

  if (is.null(datasets) || length(datasets) == 0L) {
    message(
      "No datasets specified. Downloading all available for '", source,
      "' (", length(available), " total).\n",
      "Note: loading all datasets simultaneously may require substantial memory."
    )
    to_download <- available
  } else {
    invalid <- setdiff(datasets, available)
    if (length(invalid) == length(datasets))
      stop("None of the requested datasets were found in the '", source,
           "' manifest.")
    if (length(invalid) > 0L)
      message("Datasets not found in manifest (skipped): ",
              paste(invalid, collapse = ", "), ".")
    to_download <- intersect(datasets, available)
  }

  n <- length(to_download)
  results <- lapply(seq_along(to_download), function(i) {
    ds <- to_download[i]
    message(sprintf("[%d/%d] %s", i, n, ds))
    .load_zenodo_rds(
      dataset   = ds,
      dest_dir  = dest_dir,
      source    = source,
      overwrite = overwrite,
      verify    = verify
    )
  })

  names(results) <- to_download
  results
}


#' Platform-appropriate cache directory for xEnrich
#'
#' @return Character scalar path.
#' @keywords internal
#' @noRd
xEnrich_cache_dir <- function() {
  tryCatch(
    tools::R_user_dir("xEnrich", which = "cache"),
    error = function(e) file.path(path.expand("~"), ".cache", "xEnrich")
  )
}


#' Download and cache a single Zenodo RDS file
#'
#' @param dataset Character. Dataset name (must exist in manifest).
#' @param dest_dir Character. Cache directory.
#' @param source Character. "tabula" or "hpa".
#' @param overwrite Logical.
#' @param verify Logical.
#'
#' @return The R object loaded from the cached RDS file.
#' @keywords internal
#' @noRd
.load_zenodo_rds <- function(dataset,
                             dest_dir  = xEnrich_cache_dir(),
                             source    = c("tabula", "hpa"),
                             overwrite = FALSE,
                             verify    = TRUE) {
  source   <- match.arg(source)
  manifest <- if (source == "tabula") tabula_manifest else hpa_manifest

  rows <- manifest[manifest$dataset == dataset, ]
  if (nrow(rows) == 0L)
    stop("Dataset '", dataset, "' not found in '", source, "' manifest.")
  if (nrow(rows) > 1L)
    stop("Dataset '", dataset, "' is not unique in '", source, "' manifest.")

  filename <- rows$filename
  expected <- rows$sha256

  if (verify && (is.na(expected) || !nzchar(expected)))
    stop("Missing or empty sha256 in manifest for: ", filename,
         ". Set verify = FALSE to skip (not recommended).")

  url <- paste0(
    "https://zenodo.org/api/records/",
    .ZENODO_RECORD_ID,
    "/files/", filename, "/content"
  )

  if (!dir.exists(dest_dir))
    dir.create(dest_dir, recursive = TRUE)

  dest_file <- file.path(dest_dir, filename)

  if (!file.exists(dest_file) || isTRUE(overwrite)) {
    message("  Downloading: ", filename)

    # Temporary timeout increase for large files
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = 120L)

    code <- tryCatch(
      utils::download.file(url, dest_file, mode = "wb", quiet = TRUE),
      error = function(e) NA_integer_
    )

    if (is.na(code) || code != 0L || !file.exists(dest_file))
      stop("Download failed for '", filename,
           "'. Check record ID, filename, and network connectivity.\n",
           "URL attempted: ", url)

    if (verify) {
      got <- digest::digest(dest_file, algo = "sha256", file = TRUE)
      if (!identical(tolower(got), tolower(expected)))
        stop("SHA256 mismatch for '", filename,
             "' - download may be corrupted.\n",
             "  Expected: ", expected, "\n",
             "  Got:      ", got)
    }

  } else {
    message("  Using cache: ", dest_file)
  }

  readRDS(dest_file)
}
