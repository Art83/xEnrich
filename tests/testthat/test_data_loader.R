testthat::test_that("load_reference loads all when datasets is NULL", {
  testthat::local_mocked_bindings(
    tabula_manifest = data.frame(
      dataset = c("blood", "eye"),
      filename = c("blood.rds", "eye.rds"),
      stringsAsFactors = FALSE
    ),
    hpa_manifest = data.frame(
      dataset = c("rna_consensus"),
      filename = c("rna_consensus.rds"),
      stringsAsFactors = FALSE
    ),
    load_zenodo_rds = function(dataset, dest_dir, source, overwrite, verify) {
      paste0("loaded_", dataset)
    },
    .package = "xEnrich"
  )

  out <- xEnrich::load_reference(datasets = NULL, source = "tabula")
  testthat::expect_type(out, "list")
  testthat::expect_named(out, c("blood", "eye"))
  testthat::expect_equal(out$blood, "loaded_blood")
})

testthat::test_that("load_reference errors if all datasets invalid", {
  testthat::local_mocked_bindings(
    tabula_manifest = data.frame(
      dataset = c("blood"),
      filename = c("blood.rds"),
      stringsAsFactors = FALSE
    ),
    .package = "xEnrich"
  )

  testthat::expect_error(
    xEnrich::load_reference(datasets = c("nope"), source = "tabula"),
    "No dataset found"
  )
})

testthat::test_that("load_zenodo_rds errors on non-unique dataset key", {
  testthat::local_mocked_bindings(
    tabula_manifest = data.frame(
      dataset = c("blood", "blood"),
      filename = c("a.rds", "b.rds"),
      stringsAsFactors = FALSE
    ),
    .package = "xEnrich"
  )

  testthat::expect_error(
    xEnrich:::load_zenodo_rds("blood", source = "tabula", dest_dir = tempdir()),
    "not unique"
  )
})

testthat::test_that("load_zenodo_rds uses cache when file exists and overwrite=FALSE", {
  cache <- tempfile("xEnrich-cache-")
  dir.create(cache)

  # mock manifest inside xEnrich namespace
  testthat::local_mocked_bindings(
    tabula_manifest = data.frame(
      dataset = c("blood"),
      filename = c("blood.rds"),
      stringsAsFactors = FALSE
    ),
    .package = "xEnrich"
  )

  # mock download.file inside utils namespace (because code calls utils::download.file)
  testthat::local_mocked_bindings(
    download.file = function(...) stop("download should not be called"),
    .package = "utils"
  )

  saveRDS(list(ok = TRUE), file.path(cache, "blood.rds"))

  obj <- xEnrich:::load_zenodo_rds("blood", source = "tabula", dest_dir = cache, overwrite = FALSE)
  testthat::expect_true(is.list(obj))
  testthat::expect_true(isTRUE(obj$ok))
})

testthat::test_that("load_zenodo_rds calls download.file when file missing", {
  cache <- tempfile("xEnrich-cache-")
  dir.create(cache)

  called <- 0L

  testthat::local_mocked_bindings(
    tabula_manifest = data.frame(
      dataset = c("blood"),
      filename = c("blood.rds"),
      stringsAsFactors = FALSE
    ),
    .package = "xEnrich"
  )

  testthat::local_mocked_bindings(
    download.file = function(url, destfile, mode, quiet) {
      called <<- called + 1L
      saveRDS(list(downloaded = TRUE, url = url), destfile)
      0
    },
    .package = "utils"
  )

  obj <- xEnrich:::load_zenodo_rds("blood", source = "tabula", dest_dir = cache, overwrite = FALSE)
  testthat::expect_equal(called, 1L)
  testthat::expect_true(isTRUE(obj$downloaded))
})


testthat::test_that("load_zenodo_rds errors when verify=TRUE and sha256 missing", {
  testthat::local_mocked_bindings(
    tabula_manifest = data.frame(
      dataset = "blood",
      filename = "blood.rds",
      sha256 = NA_character_,
      stringsAsFactors = FALSE
    ),
    .package = "xEnrich"
  )

  testthat::expect_error(
    xEnrich:::load_zenodo_rds("blood", source = "tabula", dest_dir = tempdir(), verify = TRUE),
    "Missing/invalid sha256"
  )
})

testthat::test_that("load_zenodo_rds skips sha256 validation when verify=FALSE", {
  cache <- tempfile("xEnrich-cache-")
  dir.create(cache)

  testthat::local_mocked_bindings(
    tabula_manifest = data.frame(
      dataset = "blood",
      filename = "blood.rds",
      sha256 = NA_character_,
      stringsAsFactors = FALSE
    ),
    .package = "xEnrich"
  )

  testthat::local_mocked_bindings(
    download.file = function(url, destfile, mode, quiet) {
      saveRDS(list(downloaded = TRUE), destfile)
      0
    },
    .package = "utils"
  )

  obj <- xEnrich:::load_zenodo_rds("blood", source = "tabula", dest_dir = cache, verify = FALSE)
  testthat::expect_true(isTRUE(obj$downloaded))
})
