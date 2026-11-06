# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 08_PPMI_Gene_Expression_Longitudinal.R
# Requires original project datasets; full R execution not validated here.
# =========================
# 0. Package
# =========================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)


# =========================
# 1. 경로 설정
# =========================

base_dir <- file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), "PPMI")

out_dir <- file.path(base_dir, "BL_V02_paired_TRPV1_analysis")

plot_dir <- file.path(out_dir, "TRPV1_related_gene_BL_V02_lineplots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)


# =========================
# 2. 파일 불러오기
# =========================

sample_info <- readr::read_csv(
  file.path(out_dir, "sample_info_BL_V02_paired.csv"),
  show_col_types = FALSE
)

tpm_df <- readr::read_csv(
  file.path(out_dir, "TPM_matrix_BL_V02_paired.csv"),
  show_col_types = FALSE
)

gene_annot <- readr::read_csv(
  file.path(out_dir, "gene_annotation_ensembl_to_symbol.csv"),
  show_col_types = FALSE
)


# =========================
# 3. target genes
# =========================

target_genes <- c(
  "NGF",
  "PLCB4",
  "MAPK10",
  "MAPK12",
  "TRPV1",
  "MAPK13",
  "IL1R1",
  "MAPK11",
  "IL1RAP",
  "PLCB3",
  "PLCG2",
  "PLCB1",
  "IL1B",
  "MAPK14",
  "MAPK9",
  "PLCB2",
  "MAPK8",
  "PLCG1"
)


# =========================
# 4. sample_info 정리
# =========================

sample_info <- sample_info %>%
  dplyr::mutate(
    PATNO = as.character(PATNO),
    sample_id = as.character(sample_id),
    Visit = factor(Visit, levels = c("BL", "V02")),
    Visit_num = dplyr::case_when(
      Visit == "BL" ~ 0,
      Visit == "V02" ~ 1
    ),
    Final_Group = factor(
      Final_Group,
      levels = c("Control", "No Pain PD", "Low Pain PD", "High Pain PD")
    )
  )


# =========================
# 5. TPM long format 만들기
# =========================

target_tpm_long <- tpm_df %>%
  dplyr::rename(Ensembl_ID = 1) %>%
  dplyr::left_join(gene_annot, by = "Ensembl_ID") %>%
  dplyr::filter(GeneSymbol %in% target_genes) %>%
  tidyr::pivot_longer(
    cols = -c(Ensembl_ID, GeneSymbol, GeneName),
    names_to = "sample_id",
    values_to = "TPM"
  ) %>%
  dplyr::left_join(sample_info, by = "sample_id") %>%
  dplyr::mutate(
    TPM = as.numeric(TPM),
    log2TPM = log2(TPM + 1)
  ) %>%
  dplyr::filter(
    !is.na(Final_Group),
    !is.na(Visit),
    !is.na(log2TPM)
  )


# Verify one value per target gene and biological sample before summary.
# Never treat multiple Ensembl rows as independent patients.
replicate_qc <- target_tpm_long %>%
  dplyr::count(GeneSymbol, sample_id, name = "n_rows") %>%
  dplyr::filter(n_rows > 1L)
if (nrow(replicate_qc) > 0L) {
  readr::write_csv(replicate_qc, file.path(plot_dir, "QC_duplicate_gene_symbol_per_sample.csv"))
  stop("Multiple Ensembl rows for one GeneSymbol/sample; inspect QC_duplicate_gene_symbol_per_sample.csv. No patient pseudo-replication allowed.")
}

# =========================
# 6. Group-visit means and 95% confidence intervals (NOT mean +/- SE)
# Pointwise two-sided 95% Student-t CI; n = non-missing observations.
# One patient should contribute one row per gene and visit (check mapping QC).
# If n < 2, show no CI ribbon for that group/visit.
# =========================

line_summary <- target_tpm_long %>%
  dplyr::group_by(GeneSymbol, Final_Group, Visit, Visit_num) %>%
  dplyr::summarise(
    n = sum(!is.na(log2TPM)),
    mean_log2TPM = mean(log2TPM, na.rm = TRUE),
    sd_log2TPM = sd(log2TPM, na.rm = TRUE),
    se_log2TPM = ifelse(n > 1L, sd_log2TPM / sqrt(n), NA_real_),
    ci95_halfwidth = ifelse(n > 1L, stats::qt(0.975, df = pmax(n - 1L, 1L)) * se_log2TPM, NA_real_),
    ymin = mean_log2TPM - ci95_halfwidth,
    ymax = mean_log2TPM + ci95_halfwidth,
    .groups = "drop"
  )

readr::write_csv(
  line_summary,
  file.path(plot_dir, "TRPV1_related_genes_BL_V02_line_summary.csv")
)


# =========================
# 7. 색상 설정
# =========================

group_colors <- c(
  "Control" = "#648FFF",
  "No Pain PD" = "#f5ff96",
  "Low Pain PD" = "#FA8072",
  "High Pain PD" = "#ff5b47"
)

group_fill <- c(
  "Control" = "#648FFF33",
  "No Pain PD" = "#f5ff9633",
  "Low Pain PD" = "#FA807233",
  "High Pain PD" = "#ff5b4733"
)

# 그래프 표시용 label만 변경 (원본 Final_Group level은 그대로 유지)
group_labels <- c(
  "Control" = " HC",
  "No Pain PD" = "No-pain PD",
  "Low Pain PD" = "Low-pain PD",
  "High Pain PD" = "High-pain PD"
)


# =========================
# 8. gene별 line plot 함수
# =========================

plot_gene_line <- function(gene_name) {
  
  df_gene <- line_summary %>%
    dplyr::filter(GeneSymbol == gene_name)
  
  p <- ggplot(
    df_gene,
    aes(
      x = Visit_num,
      y = mean_log2TPM,
      color = Final_Group,
      fill = Final_Group,
      group = Final_Group
    )
  ) +
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax),
      alpha = 0.20,
      color = NA
    ) +
    geom_line(
      linewidth = 1.5
    ) +
    geom_point(
      size = 3.2
    ) +
    scale_color_manual(values = group_colors, labels = group_labels) +
    scale_fill_manual(values = group_fill, labels = group_labels) +
    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("Baseline", "6 months"),
      limits = c(-0.08, 1.08)
    ) +
    labs(
      title = paste0(gene_name),
      x = "Time",
      y = "log2(TPM + 1)",
      color = "",
      fill = ""
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.9, color = "black"),
      axis.ticks = element_line(linewidth = 0.8, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
      panel.grid.minor = element_blank()
    )
  
  safe_gene <- stringr::str_replace_all(gene_name, "[^A-Za-z0-9_]", "_")
  
  ggsave(
    filename = file.path(plot_dir, paste0("Lineplot_", safe_gene, "_BL_to_V02.png")),
    plot = p,
    width = 8.2,
    height = 4.6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0("Lineplot_", safe_gene, "_BL_to_V02.pdf")),
    plot = p,
    width = 8.2,
    height = 4.6
  )
  
  print(p)
  
  return(p)
}


# =========================
# 9. 모든 target gene 따로 저장
# =========================

gene_list <- sort(unique(target_tpm_long$GeneSymbol))

line_plots <- lapply(gene_list, plot_gene_line)
names(line_plots) <- gene_list

cat("완료: target gene별 BL-V02 line plot 저장 완료\n")
