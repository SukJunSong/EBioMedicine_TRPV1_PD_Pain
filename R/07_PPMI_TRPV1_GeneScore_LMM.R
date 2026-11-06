# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 07_PPMI_TRPV1_GeneScore_LMM.R
# Requires original project datasets; full R execution not validated here.
# =========================
# 0. Package
# =========================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)

if (!requireNamespace("lme4", quietly = TRUE)) {
  install.packages("lme4")
}
if (!requireNamespace("lmerTest", quietly = TRUE)) {
  install.packages("lmerTest")
}
if (!requireNamespace("emmeans", quietly = TRUE)) {
  install.packages("emmeans")
}

library(lme4)
library(lmerTest)
library(emmeans)


# =========================
# 1. 경로 설정
# =========================

base_dir <- file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), "PPMI")
out_dir <- file.path(base_dir, "BL_V02_paired_TRPV1_analysis")

score_dir <- file.path(out_dir, "TRPV1_related_gene_score_BL_V02_lineplot")
dir.create(score_dir, showWarnings = FALSE, recursive = TRUE)


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
# Visit 값은 BL/V02 유지
# =========================

sample_info <- sample_info %>%
  dplyr::mutate(
    PATNO = as.character(PATNO),
    sample_id = as.character(sample_id),
    Visit = as.character(Visit),
    Visit = factor(Visit, levels = c("BL", "V02")),
    Visit_num = dplyr::case_when(
      Visit == "BL" ~ 0,
      Visit == "V02" ~ 1,
      TRUE ~ NA_real_
    ),
    Final_Group = factor(
      Final_Group,
      levels = c("Control", "No Pain PD", "Low Pain PD", "High Pain PD")
    ),
    GENDER = factor(GENDER),
    ENROLL_AGE = as.numeric(ENROLL_AGE)
  )

print(table(sample_info$Final_Group, sample_info$Visit, useNA = "ifany"))


# =========================
# 5. target gene expression long format
# =========================

tpm_annotated <- tpm_df %>%
  dplyr::rename(Ensembl_ID = 1) %>%
  dplyr::left_join(gene_annot, by = "Ensembl_ID")

id_cols <- c("Ensembl_ID", "GeneSymbol")
if ("GeneName" %in% colnames(tpm_annotated)) {
  id_cols <- c("Ensembl_ID", "GeneSymbol", "GeneName")
}

target_expr_long <- tpm_annotated %>%
  dplyr::filter(GeneSymbol %in% target_genes) %>%
  tidyr::pivot_longer(
    cols = -dplyr::all_of(id_cols),
    names_to = "sample_id",
    values_to = "TPM"
  ) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    TPM = as.numeric(TPM),
    log2TPM = log2(TPM + 1)
  ) %>%
  dplyr::left_join(sample_info, by = "sample_id") %>%
  dplyr::filter(
    !is.na(Final_Group),
    !is.na(Visit),
    !is.na(log2TPM)
  )

cat("Target expression check:\n")
print(dim(target_expr_long))
print(table(target_expr_long$Final_Group, target_expr_long$Visit))


# Prevent duplicated gene symbols from giving one target gene extra weight.
replicate_qc <- target_expr_long %>%
  dplyr::count(GeneSymbol, sample_id, name = "n_rows") %>%
  dplyr::filter(n_rows > 1L)
if (nrow(replicate_qc) > 0L) {
  readr::write_csv(replicate_qc, file.path(score_dir, "QC_duplicate_gene_symbol_per_sample.csv"))
  stop("Multiple Ensembl rows for one GeneSymbol/sample; inspect QC_duplicate_gene_symbol_per_sample.csv. Do not silently change gene score weighting.")
}

# =========================
# 6. gene별 z-score 변환 후 sample별 score 계산
# =========================

target_expr_z <- target_expr_long %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::mutate(
    gene_z = as.numeric(scale(log2TPM))
  ) %>%
  dplyr::ungroup()

score_long <- target_expr_z %>%
  dplyr::group_by(
    PATNO,
    sample_id,
    Visit,
    Visit_num,
    Final_Group,
    ENROLL_AGE,
    GENDER
  ) %>%
  dplyr::summarise(
    TRPV1_related_score = mean(gene_z, na.rm = TRUE),
    n_genes_used = sum(!is.na(gene_z)),
    .groups = "drop"
  )

readr::write_csv(
  score_long,
  file.path(score_dir, "TRPV1_related_gene_score_BL_V02_long.csv")
)

print(table(score_long$Final_Group, score_long$Visit))


# =========================
# 7. Group-visit means and 95% confidence intervals (NOT mean +/- SE)
# Pointwise two-sided 95% Student-t CI; n is non-missing sample observations.
# If n < 2, CI cannot be estimated and will be NA (no misleading ribbon).
# =========================

score_summary <- score_long %>%
  dplyr::group_by(Final_Group, Visit, Visit_num) %>%
  dplyr::summarise(
    n = sum(!is.na(TRPV1_related_score)),
    mean_score = mean(TRPV1_related_score, na.rm = TRUE),
    sd_score = sd(TRPV1_related_score, na.rm = TRUE),
    se_score = ifelse(n > 1L, sd_score / sqrt(n), NA_real_),
    ci95_halfwidth = ifelse(n > 1L, stats::qt(0.975, df = pmax(n - 1L, 1L)) * se_score, NA_real_),
    ymin = mean_score - ci95_halfwidth,
    ymax = mean_score + ci95_halfwidth,
    .groups = "drop"
  )

readr::write_csv(
  score_summary,
  file.path(score_dir, "TRPV1_related_gene_score_BL_V02_summary.csv")
)

print(score_summary)


# =========================
# 8. Linear mixed model
# 전체적으로 BL -> V02 변화 확인
# =========================

fit_score_lmm <- lmerTest::lmer(
  TRPV1_related_score ~ Visit * Final_Group + ENROLL_AGE + GENDER + (1 | PATNO),
  data = score_long
)

anova_score <- anova(fit_score_lmm)

readr::write_csv(
  as.data.frame(anova_score) %>%
    tibble::rownames_to_column("term"),
  file.path(score_dir, "TRPV1_related_gene_score_LMM_anova.csv")
)

print(anova_score)


# =========================
# 9. 그룹별 BL vs V02 paired comparison
# =========================

emm_visit_by_group <- emmeans::emmeans(
  fit_score_lmm,
  ~ Visit | Final_Group
)

visit_contrast_by_group <- summary(
  emmeans::contrast(
    emm_visit_by_group,
    method = "revpairwise"
  ),
  infer = TRUE,
  adjust = "none"
) %>%
  as.data.frame() %>%
  dplyr::mutate(
    p_adj_BH = p.adjust(p.value, method = "BH"),
    significance_raw = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    significance_BH = dplyr::case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01 ~ "**",
      p_adj_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

readr::write_csv(
  visit_contrast_by_group,
  file.path(score_dir, "TRPV1_related_gene_score_BL_vs_V02_by_group_stats.csv")
)

print(visit_contrast_by_group)


# =========================
# 10. 색상 설정
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

# 그래프 표시용 label만 변경
# 실제 Final_Group 값은 그대로 유지됨
group_labels <- c(
  "Control" = "HC",
  "No Pain PD" = "No-pain PD",
  "Low Pain PD" = "Low-pain PD",
  "High Pain PD" = "High-pain PD"
)


# =========================
# 11. BL -> V02 전체 score line graph
# =========================

p_score_line <- ggplot(
  score_summary,
  aes(
    x = Visit_num,
    y = mean_score,
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
    linewidth = 1.6
  ) +
  geom_point(
    size = 3.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey45"
  ) +
  scale_color_manual(values = group_colors, labels = group_labels) +
  scale_fill_manual(values = group_fill, labels = group_labels) +
  scale_x_continuous(
    breaks = c(0, 1),
    labels = c("Baseline", "6 months"),
    limits = c(-0.08, 1.08)
  ) +
  labs(
    title = "TRPV1-associated gene score",
x = "Time",
    y = " mean gene-wise z-scored log2(TPM + 1)",
    color = "",
    fill = ""
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.9, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

print(p_score_line)

ggsave(
  filename = file.path(score_dir, "TRPV1_related_gene_score_BL_to_V02_lineplot.png"),
  plot = p_score_line,
  width = 8.2,
  height = 4.8,
  dpi = 300
)

ggsave(
  filename = file.path(score_dir, "TRPV1_related_gene_score_BL_to_V02_lineplot.pdf"),
  plot = p_score_line,
  width = 8.2,
  height = 4.8
)

cat("완료: TRPV1-related gene score BL-to-V02 line graph 저장 완료\n")

