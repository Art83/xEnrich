#' Load one or more reference datasets shipped via Zenodo
#' @param datasets Character vector of dataset IDs to download/load. Use NULL (default)
#'   to load all datasets available for the given source.
#' @param source Which manifest to use: hpa or tabula.
#' @param dest_dir Where to cache downloaded files. Defaults to a persistent cache dir.
#' @param overwrite Re-download even if file exists locally.
#' @return A named list of loaded R objects (one per dataset).
#' @export
load_reference <- function(datasets = NULL,
                           source = c("tabula", "hpa"),
                           dest_dir = xEnrich_cache_dir(),
                           overwrite = FALSE,
                           verify = TRUE) {

  source <- match.arg(source)

  # Use internal manifests (stored in sysdata.rda)
  manifest <- if (source == "tabula") tabula_manifest else hpa_manifest

  if (!all(c("dataset", "filename") %in% names(manifest))) {
    stop("Manifest must contain columns: 'dataset' and 'filename'.")
  }

  available <- unique(manifest$dataset)

  # NULL or empty => load all
  if (is.null(datasets) || length(datasets) == 0) {
    message(
      "No datasets specified. Downloading all available for '", source,
      "' (", length(available), " total)."
    )
    to_download <- available
  } else {
    invalid <- setdiff(datasets, available)
    if (length(invalid) == length(datasets)) {
      stop("No dataset found in manifest for source='", source, "'.")
    }
    if (length(invalid) > 0) {
      message("Datasets not found: ", paste(invalid, collapse = ", "))
    }
    to_download <- intersect(datasets, available)
  }

  results <- lapply(to_download, function(ds) {
    load_zenodo_rds(
      dataset   = ds,
      dest_dir  = dest_dir,
      source    = source,
      overwrite = overwrite,
      verify = verify
    )
  })

  names(results) <- to_download
  results
}

#' @noRd
#' @keywords internal
xEnrich_cache_dir <- function() {
  tryCatch(
    tools::R_user_dir("xEnrich", which = "cache"),
    error = function(e) file.path(path.expand("~"), ".cache", "xEnrich")
  )
}


#' @noRd
#' @keywords internal
load_zenodo_rds <- function(dataset,
                            dest_dir = xEnrich_cache_dir(),
                            source = c("tabula", "hpa"),
                            overwrite = FALSE,
                            verify = TRUE) {
  record_id = 441892

  source <- match.arg(source)
  manifest <- if (source == "tabula") tabula_manifest else hpa_manifest


  filename <- manifest$filename[manifest$dataset == dataset]

  if (length(filename) == 0) stop("Dataset '", dataset, "' not found in manifest for source='", source, "'.")
  if (length(filename) > 1) stop("Dataset key '", dataset, "' is not unique in manifest for source='", source, "'.")

  expected <- if(verify) manifest$sha256[manifest$filename == filename] else NA_character_
  if (verify && (length(expected) != 1 || is.na(expected) || !nzchar(expected))) {
    stop("Missing/invalid sha256 in manifest for file: ", filename)
  }

  url <- paste0(
    "https://sandbox.zenodo.org/api/records/",
    record_id,
    "/files/",
    filename,
    "/content"
  )

  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest_file <- file.path(dest_dir, filename)

  if (!file.exists(dest_file) || isTRUE(overwrite)) {
    message("Fetching reference: ", filename, " ...")
    code <- tryCatch({
      utils::download.file(url, dest_file, mode = "wb", quiet = TRUE)
    }, error = function(e) NA_integer_)

    if (!identical(code, 0L) || !file.exists(dest_file)) {
      stop("Download failed for ", filename, ". Check record_id, filename, and connectivity.")
    }

    if(verify){
      got <- digest::digest(dest_file, algo = "sha256", file = TRUE)
      if(!identical(tolower(got), tolower(expected))) stop("SHA256 mismatch for ", url, " (download corrupted or wrong file).")
    }
  } else {
    message("Using cache: ", dest_file)
  }


  readRDS(dest_file)
}
