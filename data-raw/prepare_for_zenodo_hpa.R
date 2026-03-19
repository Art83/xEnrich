folder_in  <- "D:/hpa/raw/"
folder_out <- "D:/hpa/rds"
dir.create(folder_out, showWarnings = FALSE, recursive = TRUE)

files <- list.files(folder_in, pattern="\\.tsv$", full.names=TRUE)

for (f in files) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, sep='\t')
  out <- file.path(folder_out, paste0(tools::file_path_sans_ext(basename(f)), ".rds"))
  saveRDS(df, out, compress = "xz")
}

# =============================================================================
# IHC-specific preprocessing
# =============================================================================
# HPA IHC data contains sub-region entries (e.g. "Endometrium 1",
# "Endometrium 2") under the same Tissue label, and multiple Cell.type
# rows per gene per region. Before Zenodo upload we:
#   1. Aggregate to max IHC level per gene per Tissue
#      (a gene expressed High in any cell type/sub-region is High in that tissue)
#   2. Drop IHC.tissue.name — redundant after aggregation
# This produces one row per Gene x Tissue x Cell.type combination at the
# highest observed expression level.

ihc_path <- file.path(folder_out, "normal_ihc_data.rds")

ihc <- readRDS(ihc_path)

levels_ordered <- c("Not detected", "Low", "Medium", "High")

  # Recode non-ordinal levels as NA and drop
  # Ascending, Descending, Not representative, N/A do not fit the
  # ordinal scale and represent <0.2% of entries combined
non_ordinal <- c("Ascending", "Descending", "Not representative", "N/A")
n_dropped_levels <- sum(ihc$Level %in% non_ordinal)
ihc <- ihc[!ihc$Level %in% non_ordinal, ]
message(sprintf("Dropped %d rows with non-ordinal levels (%s).",
                  n_dropped_levels, paste(non_ordinal, collapse = ", ")))

ihc$Level <- factor(ihc$Level, levels = levels_ordered, ordered = TRUE)

# Max Level per Gene x Tissue x Cell.type
# (collapses sub-regions like Endometrium 1 / Endometrium 2)
library(data.table)
dt <- as.data.table(ihc)
dt[, Level := factor(Level, 
                     levels = c("Not detected", "Low", "Medium", "High"), 
                     ordered = TRUE)]
df_clean <- dt[dt[, .I[which.max(Level)], by = .(Gene, Tissue)]$V1]
df_clean[, IHC.tissue.name := NULL]
df_clean <- as.data.frame(df_clean)


saveRDS(df_clean, ihc_path, compress = "xz")

message(sprintf(
  "IHC preprocessing complete: %d rows -> %d rows, IHC.tissue.name dropped.",
    nrow(ihc), nrow(df_clean)
))


# =============================================================================
# Regenerate SHA256 checksums for all rds files
# =============================================================================
library(digest)

folder <- "D:/hpa/rds"
files  <- list.files(folder, pattern="\\.rds$", full.names=TRUE)

sha <- data.frame(
  file   = basename(files),
  sha256 = vapply(files, digest, algo = "sha256", file = TRUE,
                  FUN.VALUE = character(1)),
  stringsAsFactors = FALSE
)

write.csv(sha, file.path(folder, "checksums_sha256.csv"), row.names = FALSE)

hpa_man  <- read.csv("D:/hpa/manifest_hpa.csv", stringsAsFactors = FALSE)
hpa_man$sha256 <- NULL  # drop old hashes before merge
hpa_man2 <- merge(hpa_man, sha, by.x = "filename", by.y = "file", all.x = TRUE)
write.csv(hpa_man2, "D:/hpa/manifest_hpa.csv", row.names = FALSE)

message("All SHA256 checksums updated in manifest_hpa.csv")
