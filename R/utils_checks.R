# =============================================================================
# Internal defensive-check helpers shared across the IT pipeline
# =============================================================================
# None of these functions are exported. They emit warnings or messages so
# that the caller always continues; stopping on these conditions would be
# too aggressive for interactive use.


#' Warn if a character vector contains duplicate values
#'
#' Duplicates in gene_list inflate observed overlap counts; duplicates in
#' universe distort background frequencies. This helper warns once with a
#' count and the first few offenders.
#'
#' @param x Character vector to check.
#' @param arg_name String used in the warning message (e.g. "gene_list").
#'
#' @return Invisibly returns \code{x} (unchanged).
#' @keywords internal
#' @noRd
.check_duplicates <- function(x, arg_name) {
  dups <- x[duplicated(x)]
  if (length(dups) > 0L) {
    examples <- paste(head(unique(dups), 5L), collapse = ", ")
    warning(
      sprintf(
        "`%s` contains %d duplicate value(s) (e.g. %s). ",
        arg_name, length(dups), examples
      ),
      "Duplicates will be retained but may inflate overlap counts. ",
      "Consider calling unique() before passing to this function.",
      call. = FALSE
    )
  }
  invisible(x)
}


#' Warn if gene_list and universe have suspiciously low overlap
#'
#' Less than 1% overlap almost always means a case mismatch (TP53 vs tp53)
#' or identifier type mismatch (symbol vs Ensembl ID). This is silent in
#' the main functions and causes confusingly empty results.
#'
#' @param gene_list Character vector of query genes.
#' @param universe Character vector of background genes.
#' @param threshold Fraction below which a warning is emitted. Default 0.01.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
#' @noRd
.check_case_mismatch <- function(gene_list, universe, threshold = 0.01) {
  n_overlap <- length(intersect(gene_list, universe))
  frac      <- if (length(gene_list) > 0L) n_overlap / length(gene_list) else 0

  if (frac < threshold && length(gene_list) >= 10L) {
    warning(
      sprintf(
        "Only %.1f%% of `gene_list` (%d / %d genes) found in `universe`. ",
        frac * 100, n_overlap, length(gene_list)
      ),
      "This may indicate a gene identifier or case mismatch ",
      "(e.g. 'TP53' vs 'tp53', symbols vs Ensembl IDs). ",
      "Check that both use the same identifier type and case.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Warn if pathway names contain characters that commonly break subsetting
#'
#' Characters such as brackets, slashes, and colons can cause silent failures
#' when pathway names are used as list or data frame indices.
#'
#' @param nms Character vector of pathway names (e.g. names(gene_sets)).
#' @param arg_name String used in the warning message.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
#' @noRd
.check_pathway_names <- function(nms, arg_name = "gene_sets") {
  bad <- grep("[\\[\\]\\(\\)/\\\\:@$]", nms, value = TRUE)
  if (length(bad) > 0L) {
    examples <- paste(head(bad, 3L), collapse = '", "')
    warning(
      sprintf(
        '`%s` has %d name(s) containing special characters (e.g. "%s"). ',
        arg_name, length(bad), examples
      ),
      "This may cause issues when subsetting results by pathway name. ",
      "Consider sanitising names with make.names() or gsub().",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Warn and remove zero-variance columns from an expression matrix
#'
#' scale() produces NaN for zero-variance columns, which propagates silently
#' into activity scores. This helper removes offending columns and warns.
#'
#' @param expr Numeric matrix (samples x genes).
#'
#' @return The matrix with zero-variance columns removed.
#' @keywords internal
#' @noRd
.remove_zero_variance <- function(expr) {
  col_var  <- apply(expr, 2L, var, na.rm = TRUE)
  zero_var <- which(col_var == 0 | is.na(col_var))

  if (length(zero_var) > 0L) {
    examples <- paste(head(colnames(expr)[zero_var], 5L), collapse = ", ")
    warning(
      sprintf(
        "%d gene(s) in `expr` have zero variance and will be removed ",
        length(zero_var)
      ),
      sprintf("(e.g. %s). ", examples),
      "These cannot contribute to pathway activity scores.",
      call. = FALSE
    )
    expr <- expr[, -zero_var, drop = FALSE]
  }
  expr
}


#' Warn if expression matrix appears to be transposed
#'
#' run_info_assoc expects samples x genes. If nrow >> ncol the matrix is
#' likely genes x samples (a common mistake). Warns when nrow > 5 * ncol
#' and ncol < 500 (a reasonable upper bound on sample count).
#'
#' @param expr Numeric matrix.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
#' @noRd
.check_matrix_orientation <- function(expr) {
  nr <- nrow(expr)
  nc <- ncol(expr)
  if (nr > 5L * nc && nc < 500L) {
    warning(
      sprintf(
        "`expr` has %d rows and %d columns. ", nr, nc
      ),
      "run_info_assoc expects a samples x genes matrix (rows = samples). ",
      "If your matrix is genes x samples, transpose it with t(expr).",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Warn about NA values in an expression matrix
#'
#' NAs in expr are not an error -- scale() and rowMeans handle them -- but
#' they can silently reduce effective sample size per pathway. This helper
#' warns once with a count so the user is aware.
#'
#' @param expr Numeric matrix.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
#' @noRd
.check_na_expr <- function(expr) {
  n_na <- sum(is.na(expr))
  if (n_na > 0L) {
    pct <- round(100 * n_na / length(expr), 1)
    warning(
      sprintf(
        "`expr` contains %d NA value(s) (%.1f%% of entries). ", n_na, pct
      ),
      "NAs will be excluded from activity score computation per pathway. ",
      "High NA rates may reduce reliability of MI estimates.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Warn if sample size or group size is below recommended minimums
#'
#' MI estimation from equal-frequency bins is unreliable with very few
#' samples. Minimum of 50 total samples and 30 per group are practical
#' lower bounds for stable estimates.
#'
#' @param expr Numeric matrix (samples x genes).
#' @param phenotype Vector of phenotype values/labels.
#'
#' @return Invisibly returns NULL.
#' @keywords internal
#' @noRd
.check_sample_size <- function(expr, phenotype) {
  n <- nrow(expr)

  if (n < 50L) {
    warning(
      sprintf("Only %d samples in `expr`. ", n),
      "MI estimation from equal-frequency bins is unreliable with fewer ",
      "than ~50 samples. Results should be interpreted with caution.",
      call. = FALSE
    )
  }

  if (!is.numeric(phenotype)) {
    grp_sizes <- table(phenotype)
    small_grp <- grp_sizes[grp_sizes < 30L]
    if (length(small_grp) > 0L) {
      examples <- paste(
        sprintf("%s (n=%d)", names(small_grp), as.integer(small_grp)),
        collapse = ", "
      )
      warning(
        sprintf(
          "Some phenotype groups have fewer than 30 samples: %s. ", examples
        ),
        "Small group sizes reduce MI estimation stability and may inflate ",
        "false positive rates.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}
