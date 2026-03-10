#' Summarise expression distribution of a reference data frame
#'
#' Inspects a reference data frame (HPA RNA, Tabula Sapiens, or similar)
#' to help choose \code{cutoff_type} and \code{cutoff_value} before calling
#' \code{\link{enrich_by_hpa}} or \code{\link{enrich_by_tabula}}. Returns
#' per-group expression summaries and a preview of how many genes each
#' cutoff combination would select.
#'
#' @param data Data frame in long format. Must contain \code{group_col},
#'   \code{gene_col}, and \code{expr_col}.
#' @param group_col Character. Column name for the grouping variable
#'   (e.g. \code{"Tissue"}, \code{"Brain.region"}, \code{"cell_type"}).
#' @param gene_col Character. Column name for gene identifiers.
#'   Default \code{"Gene.name"}.
#' @param expr_col Character. Column name for expression values.
#'   Default \code{"pTPM"}.
#' @param cutoff_preview Named list of cutoff specifications. Each element
#'   is a list with \code{type} (one of \code{"quantile"}, \code{"threshold"},
#'   \code{"top_k"}) and \code{value}. Defaults cover the three most common
#'   choices.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{distribution}}{Data frame with one row per group:
#'       \code{group}, \code{n_genes}, \code{median_expr},
#'       \code{n_above_p75}, \code{n_above_p90}, \code{n_above_p95}.}
#'     \item{\code{cutoff_preview}}{Data frame with one row per
#'       group-cutoff combination: \code{group}, \code{cutoff_type},
#'       \code{cutoff_value}, \code{n_selected},
#'       \code{pct_of_universe}.}
#'   }
#'
#' @examples
#' # Simulate a minimal HPA-style long data frame
#' set.seed(1)
#' ref <- data.frame(
#'   Gene.name = rep(paste0("G", 1:100), times = 3),
#'   Tissue    = rep(c("Liver", "Kidney", "Brain"), each = 100),
#'   pTPM      = c(rexp(100, 0.05), rexp(100, 0.02), rexp(100, 0.1))
#' )
#'
#' sm <- summarize_expression(ref, group_col = "Tissue")
#' sm$distribution
#' sm$cutoff_preview
#'
#' @export
summarize_expression <- function(
    data,
    group_col,
    gene_col       = "Gene.name",
    expr_col       = "pTPM",
    cutoff_preview = list(
      list(type = "quantile",  value = 0.9),
      list(type = "threshold", value = 1),
      list(type = "top_k",     value = 200)
    )
) {
  ## --- Validation ------------------------------------------------------------
  if (!is.data.frame(data))
    stop("`data` must be a data frame.")
  missing_cols <- setdiff(c(group_col, gene_col, expr_col), colnames(data))
  if (length(missing_cols) > 0L)
    stop("Missing columns in `data`: ", paste(missing_cols, collapse = ", "))
  if (!is.numeric(data[[expr_col]]))
    stop("`expr_col` ('", expr_col, "') must be numeric.")

  # Validate cutoff_preview entries
  valid_types <- c("quantile", "threshold", "top_k")
  for (i in seq_along(cutoff_preview)) {
    cp <- cutoff_preview[[i]]
    if (!all(c("type", "value") %in% names(cp)))
      stop("Each `cutoff_preview` entry must have 'type' and 'value'.")
    if (!cp$type %in% valid_types)
      stop("cutoff_preview[[", i, "]]$type must be one of: ",
           paste(valid_types, collapse = ", "))
    if (cp$type == "quantile" && (cp$value <= 0 || cp$value >= 1))
      stop("cutoff_preview[[", i, "]]: quantile value must be in (0, 1).")
    if (cp$type == "top_k" && (cp$value %% 1 != 0 || cp$value < 1))
      stop("cutoff_preview[[", i, "]]: top_k value must be a positive integer.")
  }

  ## --- Setup -----------------------------------------------------------------
  universe      <- unique(data[[gene_col]])
  n_universe    <- length(universe)
  data_split    <- split(data, data[[group_col]])

  ## --- Per-group computation -------------------------------------------------
  dist_rows    <- vector("list", length(data_split))
  preview_rows <- vector("list", length(data_split))

  for (i in seq_along(data_split)) {
    group_name <- names(data_split)[i]
    df         <- data_split[[i]]
    genes      <- unique(df[[gene_col]])
    expr       <- df[[expr_col]]

    # Distribution summary
    qs <- stats::quantile(expr, probs = c(0.5, 0.75, 0.90, 0.95), na.rm = TRUE)

    n_above <- function(q) {
      thresh <- stats::quantile(expr, q, na.rm = TRUE)
      length(unique(df[[gene_col]][df[[expr_col]] > thresh]))
    }

    dist_rows[[i]] <- data.frame(
      group        = group_name,
      n_genes      = length(genes),
      median_expr  = unname(round(qs["50%"], 3)),
      n_above_p75  = n_above(0.75),
      n_above_p90  = n_above(0.90),
      n_above_p95  = n_above(0.95),
      stringsAsFactors = FALSE
    )

    # Cutoff preview
    cp_rows <- vector("list", length(cutoff_preview))

    for (j in seq_along(cutoff_preview)) {
      cp   <- cutoff_preview[[j]]
      type <- cp$type
      val  <- cp$value

      gene_set <- switch(type,
                         quantile = {
                           thresh <- stats::quantile(expr, val, na.rm = TRUE)
                           unique(df[[gene_col]][df[[expr_col]] > thresh])
                         },
                         threshold = {
                           unique(df[[gene_col]][df[[expr_col]] > val])
                         },
                         top_k = {
                           # deduplicate by taking max expression per gene first
                           gene_max <- tapply(df[[expr_col]], df[[gene_col]], max, na.rm = TRUE)
                           names(sort(gene_max, decreasing = TRUE))[
                             seq_len(min(val, length(gene_max)))
                           ]
                         }
      )

      n_selected   <- length(gene_set)
      pct_universe <- round(
        100 * length(intersect(gene_set, universe)) / n_universe, 1
      )

      cp_rows[[j]] <- data.frame(
        group         = group_name,
        cutoff_type   = type,
        cutoff_value  = as.character(val),
        n_selected    = n_selected,
        pct_of_universe = pct_universe,
        stringsAsFactors = FALSE
      )
    }

    preview_rows[[i]] <- do.call(rbind, cp_rows)
  }

  list(
    distribution   = do.call(rbind, dist_rows),
    cutoff_preview = do.call(rbind, preview_rows)
  )
}
