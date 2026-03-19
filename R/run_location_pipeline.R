# =============================================================================
# Location-aware enrichment pipeline
# =============================================================================


#' Extract top organ hits from an HPA enrichment result
#'
#' Selects rows from an [enrich_by_hpa_ihc()] or [enrich_by_hpa()] result
#' that pass a significance threshold, returning their tissue/group labels.
#'
#' @param hpa_result Data frame returned by [enrich_by_hpa_ihc()] or
#'   [enrich_by_hpa()].
#' @param group_col Character. Column containing organ/tissue labels.
#' @param top_n Integer. Maximum organs to return.
#' @param max_padj Numeric. Adjusted p-value ceiling. Default 0.05.
#'
#' @return Character vector of organ labels, in ascending padj order.
#' @keywords internal
#' @noRd
.top_organs <- function(hpa_result, group_col, top_n, max_padj) {
  if (!"padj" %in% colnames(hpa_result))
    stop("HPA result has no `padj` column. Was it produced by enrich_by_hpa_ihc() ",
         "or enrich_by_hpa()?")
  if (!group_col %in% colnames(hpa_result))
    stop("Column '", group_col, "' not found in HPA result.")

  hits <- hpa_result[hpa_result$padj <= max_padj, ]
  hits <- hits[order(hits$padj), ]
  head(unique(hits[[group_col]]), top_n)
}


#' Map HPA tissue names to Tabula dataset names via the bridge table
#'
#' Returns a data frame of matched pairs, flagging approximate mappings and
#' noting organs with no Tabula equivalent.
#'
#' @param organs Character vector of HPA tissue names.
#'
#' @return A data frame with columns \code{hpa_tissue}, \code{tabula},
#'   \code{approximate}.
#' @keywords internal
#' @noRd
.map_organs_to_tabula <- function(organs) {
  # tabula_hpa_bridge and tabula_no_hpa_bridge live in sysdata.rda
  organs <- tolower(organs)
  matched   <- tabula_hpa_bridge[tolower(tabula_hpa_bridge$hpa_tissue) %in% organs, ]
  no_bridge <- organs[
    !organs %in% tolower(tabula_hpa_bridge$hpa_tissue) &
      !organs %in% tolower(tabula_no_hpa_bridge)
  ]
  known_absent <- organs[organs %in% tabula_no_hpa_bridge]

  if (length(known_absent) > 0L)
    message("No Tabula equivalent for: ",
            paste(known_absent, collapse = ", "), " (known gap.skipped).")
  if (length(no_bridge) > 0L)
    message("No bridge entry found for: ",
            paste(no_bridge, collapse = ", "),
            " - skipped. Consider updating tabula_hpa_bridge.")

  # Deduplicate tabula datasets (l.intestine can appear via Colon + Rectum)
  matched[!duplicated(matched$tabula), ]
}


#' Location-aware enrichment pipeline
#'
#' Implements a two-step enrichment strategy that moves from organ-level
#' signal to cell-type resolution. Step 1 uses Human Protein Atlas (HPA)
#' data (IHC or RNA) to identify enriched organs. Step 2 drills into the
#' top organs using Tabula Sapiens scRNA data to identify which cell types
#' drive the signal.
#'
#' This pipeline operationalises the central dogma-aware enrichment concept:
#' IHC data provides protein-level organ evidence; scRNA data resolves
#' cellular origin. Together they answer not just "what pathway" but "where
#' in the body and in which cell type."
#'
#' @section Pipeline steps:
#' \enumerate{
#'   \item Run [enrich_by_hpa_ihc()] (if \code{hpa_type = "ihc"}) or
#'     [enrich_by_hpa()] (if \code{hpa_type = "rna"}) on \code{hpa_data}.
#'   \item Select up to \code{top_n_organs} organs with
#'     \code{padj <= max_organ_padj}, ordered by adjusted p-value.
#'   \item Map organ names to Tabula dataset names via the internal
#'     \code{tabula_hpa_bridge} lookup. Approximate mappings (e.g.
#'     \code{trachea} → \code{"Bronchus"}) are flagged in the output.
#'   \item For each matched organ, run [enrich_by_tabula()] on the
#'     corresponding dataset (loaded from \code{tabula_data}).
#'   \item Return a structured list with organ-level and cell-type-level
#'     results.
#' }
#'
#' @section Providing Tabula data:
#' Pass pre-loaded datasets as a named list via \code{tabula_data}:
#' \code{list(kidney = load_reference("kidney", source = "tabula")[[1]], ...)}.
#' Dataset names must match the \code{tabula} column of the bridge table
#' (e.g. \code{"kidney"}, \code{"liver"}). If a required organ is absent
#' from \code{tabula_data}, it is skipped with a message.
#'
#' @section Argument passing:
#' Arguments for the HPA step and the Tabula step are passed as named lists
#' via \code{hpa_args} and \code{tabula_args} respectively. This avoids
#' ambiguity when both steps share parameter names (e.g. \code{method},
#' \code{universe_type}). Example:
#' \preformatted{
#' run_location_pipeline(
#'   gene_list   = my_genes,
#'   hpa_data    = ihc_data,
#'   tabula_data = tabula_list,
#'   hpa_args    = list(method = "hypergeometric", n_perm = 1000),
#'   tabula_args = list(method = "gsea", min_pct_expr = 0.15)
#' )
#' }
#'
#' @param gene_list Character vector of gene symbols (the query set).
#' @param hpa_data Data frame of HPA data passed to the HPA enrichment step.
#' @param tabula_data Named list of Tabula Sapiens data frames, one per organ.
#'   Names must match the \code{tabula} column of the bridge table. See
#'   [load_reference()] to load datasets.
#' @param hpa_type Character. \code{"ihc"} (default, uses
#'   [enrich_by_hpa_ihc()]) or \code{"rna"} (uses [enrich_by_hpa()]).
#' @param top_n_organs Integer. Maximum number of significant organs to drill
#'   into. Default \code{3}.
#' @param max_organ_padj Numeric. Adjusted p-value ceiling for organ
#'   selection. Only organs with \code{padj <= max_organ_padj} proceed to
#'   the Tabula step. Default \code{0.05}.
#' @param hpa_args Named list of arguments passed to [enrich_by_hpa_ihc()]
#'   or [enrich_by_hpa()]. Do not include \code{gene_list} or
#'   \code{reference_df} these are supplied internally.
#' @param tabula_args Named list of arguments passed to [enrich_by_tabula()].
#'   Do not include \code{gene_list} or \code{reference_df}.
#' @param seed Optional integer. Applied once before the pipeline starts,
#'   ensuring reproducibility across both steps.
#'
#' @return An object of class \code{xEnrich_pipeline}, a list containing:
#'   \describe{
#'     \item{\code{organ_results}}{Data frame of HPA enrichment results across
#'       all groups (full output of the HPA step).}
#'     \item{\code{celltype_results}}{Named list of data frames, one per
#'       organ that was successfully mapped and tested. Each data frame is
#'       the output of [enrich_by_tabula()] for that organ.}
#'     \item{\code{selected_organs}}{Character vector of organs that passed
#'       \code{max_organ_padj} and were mapped to Tabula.}
#'     \item{\code{bridge_used}}{Data frame of the organ→Tabula mappings
#'       actually used, including the \code{approximate} flag.}
#'     \item{\code{hpa_type}}{Character. HPA step type used.}
#'     \item{\code{call_args}}{List recording \code{top_n_organs},
#'       \code{max_organ_padj}, \code{hpa_args}, \code{tabula_args} for
#'       reproducibility.}
#'   }
#'
#' @seealso [enrich_by_hpa_ihc()], [enrich_by_hpa()], [enrich_by_tabula()],
#'   [load_reference()], [print.xEnrich_pipeline()],
#'   [summary.xEnrich_pipeline()]
#'
#' @examples
#' \dontrun{
#' # Load data
#' hpa_ihc   <- load_reference("ihc",    source = "hpa")[[1]]
#' tabula_kd <- load_reference("kidney", source = "tabula")[[1]]
#' tabula_lv <- load_reference("liver",  source = "tabula")[[1]]
#'
#' gene_list <- c("UMOD", "SLC12A1", "NPHS1", "PODXL", "CDH16")
#'
#' result <- run_location_pipeline(
#'   gene_list   = gene_list,
#'   hpa_data    = hpa_ihc,
#'   tabula_data = list(kidney = tabula_kd, liver = tabula_lv),
#'   hpa_type    = "ihc",
#'   top_n_organs = 3,
#'   hpa_args    = list(method = "hypergeometric"),
#'   tabula_args = list(method = "hypergeometric", min_pct_expr = 0.1),
#'   seed        = 42
#' )
#'
#' print(result)
#' result$organ_results
#' result$celltype_results$kidney
#' }
#'
#' @export
run_location_pipeline <- function(
    gene_list,
    hpa_data,
    tabula_data,
    hpa_type       = c("ihc", "rna"),
    top_n_organs   = 3L,
    max_organ_padj = 0.05,
    hpa_args       = list(),
    tabula_args    = list(),
    seed           = NULL
) {
  # --- Validation -------------------------------------------------------------
  if (!is.character(gene_list) || length(gene_list) == 0L)
    stop("`gene_list` must be a non-empty character vector.")
  if (!is.data.frame(hpa_data))
    stop("`hpa_data` must be a data frame.")
  if (!is.list(tabula_data) || is.null(names(tabula_data)))
    stop("`tabula_data` must be a named list of data frames.")
  if (!is.numeric(top_n_organs) || top_n_organs < 1L)
    stop("`top_n_organs` must be a positive integer.")
  if (!is.numeric(max_organ_padj) || max_organ_padj <= 0 ||
      max_organ_padj > 1)
    stop("`max_organ_padj` must be in (0, 1].")
  if (!is.list(hpa_args))
    stop("`hpa_args` must be a named list.")
  if (!is.list(tabula_args))
    stop("`tabula_args` must be a named list.")

  # Disallow passing gene_list / reference_df through the arg lists
  reserved_hpa    <- c("gene_list", "reference_df")
  reserved_tabula <- c("gene_list", "reference_df")
  if (any(reserved_hpa %in% names(hpa_args)))
    stop("Do not pass `gene_list` or `reference_df` via `hpa_args` - ",
         "these are supplied by the pipeline.")
  if (any(reserved_tabula %in% names(tabula_args)))
    stop("Do not pass `gene_list` or `reference_df` via `tabula_args` - ",
         "these are supplied by the pipeline.")

  hpa_type <- match.arg(hpa_type)
  if (!is.null(seed)) set.seed(seed)

  # --- Step 1: HPA enrichment -------------------------------------------------
  message("=== Step 1: HPA ", toupper(hpa_type), " enrichment ===")

  hpa_fn <- if (hpa_type == "ihc") enrich_by_hpa_ihc else enrich_by_hpa

  organ_results <- tryCatch(
    do.call(hpa_fn, c(
      list(gene_list    = gene_list,
           reference_df = hpa_data),
      hpa_args
    )),
    error = function(e)
      stop("HPA enrichment step failed: ", conditionMessage(e))
  )

  if (nrow(organ_results) == 0L) {
    message("No significant organ results from HPA step. Pipeline stopping.")
    return(structure(
      list(
        organ_results    = organ_results,
        celltype_results = list(),
        selected_organs  = character(0),
        bridge_used      = data.frame(),
        hpa_type         = hpa_type,
        call_args        = list(
          top_n_organs   = top_n_organs,
          max_organ_padj = max_organ_padj,
          hpa_args       = hpa_args,
          tabula_args    = tabula_args
        )
      ),
      class = "xEnrich_pipeline"
    ))
  }

  # Detect group column from HPA result
  candidate_cols <- c("Tissue", "IHC.tissue.name", "Cell.type",
                      "Brain.region", "Subregion", "Immune.cell")
  group_col <- candidate_cols[candidate_cols %in% colnames(organ_results)][1L]
  if (is.na(group_col))
    stop("Cannot detect group column in HPA result. Expected one of: ",
         paste(candidate_cols, collapse = ", "), ".")

  # Select top organs
  selected_organs <- .top_organs(
    hpa_result = organ_results,
    group_col  = group_col,
    top_n      = as.integer(top_n_organs),
    max_padj   = max_organ_padj
  )

  if (length(selected_organs) == 0L) {
    message("No organs passed max_organ_padj = ", max_organ_padj,
            ". Pipeline stopping after HPA step.")
    return(structure(
      list(
        organ_results    = organ_results,
        celltype_results = list(),
        selected_organs  = character(0),
        bridge_used      = data.frame(),
        hpa_type         = hpa_type,
        call_args        = list(
          top_n_organs   = top_n_organs,
          max_organ_padj = max_organ_padj,
          hpa_args       = hpa_args,
          tabula_args    = tabula_args
        )
      ),
      class = "xEnrich_pipeline"
    ))
  }

  message("Selected organs (padj <= ", max_organ_padj, "): ",
          paste(selected_organs, collapse = ", "))

  # --- Step 2: Map organs to Tabula datasets ----------------------------------
  message("\n=== Step 2: Mapping organs to Tabula datasets ===")
  bridge_used <- .map_organs_to_tabula(selected_organs)

  if (nrow(bridge_used) == 0L) {
    message("No selected organs could be mapped to Tabula datasets.")
    return(structure(
      list(
        organ_results    = organ_results,
        celltype_results = list(),
        selected_organs  = selected_organs,
        bridge_used      = bridge_used,
        hpa_type         = hpa_type,
        call_args        = list(
          top_n_organs   = top_n_organs,
          max_organ_padj = max_organ_padj,
          hpa_args       = hpa_args,
          tabula_args    = tabula_args
        )
      ),
      class = "xEnrich_pipeline"
    ))
  }

  # Warn about approximate mappings
  approx <- bridge_used$tabula[bridge_used$approximate]
  if (length(approx) > 0L)
    warning("Approximate organ mappings used for: ",
            paste(approx, collapse = ", "),
            ". Check `result$bridge_used` for details.")

  # --- Step 3: Tabula enrichment per organ ------------------------------------
  message("\n=== Step 3: Tabula cell-type enrichment ===")

  celltype_results <- list()

  for (i in seq_len(nrow(bridge_used))) {
    tab_name  <- bridge_used$tabula[i]
    hpa_name  <- bridge_used$hpa_tissue[i]
    is_approx <- bridge_used$approximate[i]

    # Check dataset is available
    if (!tab_name %in% names(tabula_data)) {
      message("  [", tab_name, "] not found in `tabula_data` - skipped. ",
              "Load with load_reference('", tab_name, "', source = 'tabula').")
      next
    }

    approx_note <- if (is_approx) {
      paste0(" [approximate mapping from '", hpa_name, "']")
    } else ""

    message("  [", i, "/", nrow(bridge_used), "] ",
            tab_name, approx_note)

    res <- tryCatch(
      do.call(enrich_by_tabula, c(
        list(gene_list    = gene_list,
             reference_df = tabula_data[[tab_name]]),
        tabula_args
      )),
      error = function(e) {
        message("  ERROR for ", tab_name, ": ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(res) && nrow(res) > 0L) {
      res$organ          <- tab_name
      res$hpa_tissue     <- hpa_name
      res$approx_mapping <- is_approx
      celltype_results[[tab_name]] <- res
    } else {
      message("  No cell types passed filters for ", tab_name, ".")
    }
  }

  # --- Assemble output --------------------------------------------------------
  n_ct_hits <- sum(vapply(celltype_results,
                          function(x) sum(x$padj <= 0.05, na.rm = TRUE),
                          integer(1L)))

  message("\n=== Pipeline complete ===")
  message("  Organs tested:       ", nrow(bridge_used))
  message("  Organs with results: ", length(celltype_results))
  message("  Cell type hits (padj <= 0.05): ", n_ct_hits)

  structure(
    list(
      organ_results    = organ_results,
      celltype_results = celltype_results,
      selected_organs  = selected_organs,
      bridge_used      = bridge_used,
      hpa_type         = hpa_type,
      call_args        = list(
        top_n_organs   = top_n_organs,
        max_organ_padj = max_organ_padj,
        hpa_args       = hpa_args,
        tabula_args    = tabula_args
      )
    ),
    class = "xEnrich_pipeline"
  )
}


# =============================================================================
# S3 methods for xEnrich_pipeline
# =============================================================================

#' Print method for xEnrich_pipeline
#'
#' @param x An \code{xEnrich_pipeline} object.
#' @param ... Ignored.
#' @export
print.xEnrich_pipeline <- function(x, ...) {
  cat("xEnrich location pipeline result\n")
  cat(rep("-", 40), "\n", sep = "")
  cat("HPA step:       ", toupper(x$hpa_type), "\n")
  cat("Organs tested:  ", nrow(x$organ_results), "\n")
  cat("Organs selected:", length(x$selected_organs), "\n")

  if (length(x$selected_organs) > 0L) {
    cat("  ", paste(x$selected_organs, collapse = ", "), "\n")
  }

  if (nrow(x$bridge_used) > 0L) {
    approx <- x$bridge_used$tabula[x$bridge_used$approximate]
    if (length(approx) > 0L)
      cat("  Approximate mappings: ",
          paste(approx, collapse = ", "), "\n")
  }

  cat("Tabula datasets with results:", length(x$celltype_results), "\n")

  if (length(x$celltype_results) > 0L) {
    for (nm in names(x$celltype_results)) {
      df      <- x$celltype_results[[nm]]
      n_sig   <- sum(df$padj <= 0.05, na.rm = TRUE)
      top_ct  <- if (nrow(df) > 0L) df$cell_type[1L] else "-"
      cat(sprintf("  %-20s %d sig cell types  (top: %s)\n",
                  nm, n_sig, top_ct))
    }
  }
  invisible(x)
}


#' Summary method for xEnrich_pipeline
#'
#' Returns a flat data frame of the top cell type hit per organ, suitable
#' for reporting or downstream visualisation.
#'
#' @param object An \code{xEnrich_pipeline} object.
#' @param max_padj Numeric. Adjusted p-value ceiling for inclusion.
#'   Default \code{0.05}.
#' @param ... Ignored.
#'
#' @return A data frame with columns \code{organ}, \code{hpa_tissue},
#'   \code{approx_mapping}, \code{cell_type}, \code{p_value}, \code{padj},
#'   and method-specific columns.
#' @export
summary.xEnrich_pipeline <- function(object, max_padj = 0.05, ...) {
  if (length(object$celltype_results) == 0L) {
    message("No cell-type results available.")
    return(invisible(data.frame()))
  }

  rows <- lapply(object$celltype_results, function(df) {
    df[df$padj <= max_padj, ]
  })
  rows <- Filter(function(df) nrow(df) > 0L, rows)

  if (length(rows) == 0L) {
    message("No cell types pass max_padj = ", max_padj, ".")
    return(invisible(data.frame()))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$organ, out$padj), ]
}
