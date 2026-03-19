folder_in  <- "D:/tabula/processed_long/"
folder_out <- "D:/tabula/rds/"
dir.create(folder_out, showWarnings = FALSE, recursive = TRUE)

files <- list.files(folder_in, pattern="\\.csv$", full.names=TRUE)

for (f in files) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)

  # basic type hygiene
  df$GeneID <- as.character(df$GeneID)
  if ("GeneSymbol" %in% names(df)) df$GeneSymbol <- as.character(df$GeneSymbol)
  df$CellType <- as.character(df$CellType)

  num_cols <- intersect(c("MeanLogNorm","PctExpr","Expression"), names(df))
  for (cc in num_cols) df[[cc]] <- as.numeric(df[[cc]])

  # Strip Ensembl version suffix (e.g. ENSG00000000003.15 -> ENSG00000000003)
  df$GeneID <- sub("\\.\\d+$", "", df$GeneID)

  out <- file.path(folder_out, paste0(tools::file_path_sans_ext(basename(f)), ".rds"))
  saveRDS(df, out, compress = "xz")
}


library(digest)

folder <- "D:/tabula/rds"
files  <- list.files(folder, pattern="\\.rds$", full.names=TRUE)

sha <- data.frame(
  file   = basename(files),
  sha256 = vapply(files, digest, algo="sha256", file=TRUE, FUN.VALUE=character(1)),
  stringsAsFactors = FALSE
)

write.csv(sha, file.path(folder, "checksums_sha256.csv"), row.names=FALSE)

tab_man        <- read.csv("D:/tabula/manifest_tabula.csv", stringsAsFactors = FALSE)
tab_man$sha256 <- NULL  # drop stale hashes before merge
tab_man2       <- merge(tab_man, sha, by.x = "filename", by.y = "file", all.x = TRUE)

write.csv(tab_man2, "D:/tabula/manifest_tabula.csv", row.names = FALSE)

