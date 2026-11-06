# Input/output paths are configured via environment variables: PD_TRPV1_GSE8397_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""))) stop("Set PD_TRPV1_GSE8397_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 19_GSE8397_GSEA_Dotplot.R
# Requires original project datasets; full R execution not validated here.
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# =========================
# 1. 데이터 불러오기
# =========================
gsea_df <- read_csv(file.path(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""), 'pain_inflam_GO_GSEA_by_region_with_desc.csv'))

# =========================
# 2. 저장 폴더
# =========================
out_dir <- file.path(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""), 'GSEA_region_dotplots')

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# =========================
# 3. region 목록
# =========================
regions <- unique(gsea_df$region)

# =========================
# 4. region별 plot 저장
# =========================
for (reg in regions) {
  
  plot_df <- gsea_df %>%
    filter(
      region == reg,
      !is.na(NES),
      !is.na(pval),
      !is.na(padj),
      !is.na(size),
      !is.na(GO_Description)
    ) %>%
    arrange(padj) %>%
    slice_head(n = 15) %>%
    mutate(
      GO_Description_clean = paste0(str_wrap(GO_Description, width = 45), ifelse(padj < 0.05, " *", "")),
      GO_Description_clean = factor(
        GO_Description_clean,
        levels = GO_Description_clean[order(NES)]
      )
    )
  
  p <- ggplot(
    plot_df,
    aes(
      x = GO_Description_clean,
      y = NES,
      size = size,
      color = NES
    )
  ) +
    geom_point(alpha = 0.9) +
    coord_flip() +
    scale_color_gradient2(
      low = "#648FFF",
      mid = "white",
      high = "#FA8072",
      midpoint = 0,
      name = "NES"
    ) +
    scale_size_continuous(
      range = c(3, 9),
      name = "Size"
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = "NES",
      caption = "* Benjamini-Hochberg FDR < 0.05"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      # 글자
      axis.title.x = element_text(face = "bold", size = 15, color = "black"),
      axis.title.y = element_text(face = "bold", size = 15, color = "black"),
      axis.text.x  = element_text(face = "bold", size = 12, color = "black"),
      axis.text.y  = element_text(face = "bold", size = 11, color = "black"),
      legend.title = element_text(face = "bold", size = 12, color = "black"),
      legend.text  = element_text(face = "bold", size = 10, color = "black"),
      
      # guide line
      panel.grid.major.x = element_line(color = "grey75", linewidth = 0.5),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5),
      panel.grid.minor.x = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor.y = element_blank(),
      
      # 바깥 테두리
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      
      # 배경
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(p)
  
  # 파일명 정리
  reg_clean <- gsub("[^A-Za-z0-9_\\-]", "_", reg)
  
  ggsave(
    filename = file.path(out_dir, paste0("GO_GSEA_dotplot_top15_", reg_clean, ".pdf")),
    plot = p,
    width = 9,
    height = 6.5,
    device = cairo_pdf
  )
}
