# Input/output paths are configured via environment variables: PD_TRPV1_GSE8397_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""))) stop("Set PD_TRPV1_GSE8397_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 18_GSE8397_GO_Description.R
# Requires original project datasets; full R execution not validated here.

# ============================================================
# GSE8397 GO term description annotation
# ============================================================
# Input:
#   pain_inflam_GO_GSEA_by_region.csv
#
# Output:
#   pain_inflam_GO_GSEA_by_region_with_desc.csv
#
# Purpose:
#   Add GO term descriptions without changing GSEA statistics.
# ============================================================

# 1. Required packages
library(readr)
library(dplyr)
library(GO.db)

# 2. File paths
input_csv <- file.path(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""), 'pain_inflam_GO_GSEA_by_region.csv')

output_csv <- file.path(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""), 'pain_inflam_GO_GSEA_by_region_with_desc.csv')

# 3. Read GSEA results
df <- readr::read_csv(
  input_csv,
  show_col_types = FALSE
)

# 4. Check required columns
required <- c(
  "region",
  "pathway",
  "NES",
  "pval",
  "padj",
  "size"
)

missing_cols <- setdiff(required, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# 5. GO ID to GO description
get_go_description <- function(id) {
  
  if (is.na(id) ||
      !grepl("^GO:[0-9]{7}$", id)) {
    return(NA_character_)
  }
  
  term <- GO.db::GOTERM[[id]]
  
  if (is.null(term)) {
    return(NA_character_)
  }
  
  return(
    as.character(GO.db::Term(term))
  )
}

# 6. Add GO descriptions
df <- df %>%
  dplyr::mutate(
    GO_ID = as.character(pathway),
    
    GO_Description = vapply(
      GO_ID,
      get_go_description,
      character(1)
    )
  )

# 7. Check unmapped GO terms
missing_go <- df %>%
  dplyr::filter(is.na(GO_Description)) %>%
  dplyr::distinct(GO_ID)

cat(
  "Number of GO IDs without descriptions:",
  nrow(missing_go),
  "\n"
)

if (nrow(missing_go) > 0) {
  print(missing_go)
}

# 8. Save annotated GSEA results
readr::write_csv(
  df,
  output_csv
)

cat(
  "Saved:",
  output_csv,
  "\n"
)

# 9. Preview
print(
  df %>%
    dplyr::select(
      region,
      GO_ID,
      GO_Description,
      NES,
      pval,
      padj,
      size
    ) %>%
    head(10)
)
