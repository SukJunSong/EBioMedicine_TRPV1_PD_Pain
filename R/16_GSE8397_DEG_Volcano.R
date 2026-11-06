# Input/output paths are configured via environment variables: PD_TRPV1_ROOT, PD_TRPV1_GSE8397_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)
if (!nzchar(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""))) stop("Set PD_TRPV1_GSE8397_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 16_GSE8397_DEG_Volcano.R
# Requires original project datasets; full R execution not validated here.
############################################################
# GSE8397 Volcano plot
# Final clean version
#
# Rules:
# 1. Expression file path fixed
# 2. Contrast fixed as PD - Control
# 3. If multiple probes map to the same gene,
#    keep only the probe with the lowest P.Value
# 4. Label only KEGG hsa04750 genes with P.Value < 0.05
# 5. Avoid duplicate gene labels
############################################################

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(GEOquery)
library(limma)
library(grid)
library(scales)

############################################################
# 0. File paths
############################################################

expr_candidates <- c(
  file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'GSE8397-GPL96_series_matrix.txt','GSE8397-GPL96_series_matrix.txt'),
  file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'GSE8397-GPL96_series_matrix.txt'),
  file.path(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""), 'GSE8397-GPL96_series_matrix.txt'),
  file.path(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""), 'GSE8397-GPL96_series_matrix.txt','GSE8397-GPL96_series_matrix.txt')
)
expr_existing <- expr_candidates[file.exists(expr_candidates)]
if (!length(expr_existing)) stop("GSE8397 series matrix missing. Checked: ", paste(expr_candidates, collapse=" | "))
expr_file <- expr_existing[1]

out_dir <- file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'GSE8397_volcano_final')

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

############################################################
# 1. Load phenotype from GSE8397
############################################################

gse <- getGEO("GSE8397", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

cat("\nPhenotype columns:\n")
print(colnames(pheno))

############################################################
# 2. Load GPL96 annotation from GEO
############################################################

gpl <- getGEO("GPL96")
gpl_table <- Table(gpl)

cat("\nGPL96 annotation columns:\n")
print(colnames(gpl_table))

possible_gene_cols <- c(
  "Gene Symbol",
  "Gene.symbol",
  "Gene.Symbol",
  "GENE_SYMBOL",
  "Symbol",
  "gene_symbol"
)

gene_col <- possible_gene_cols[possible_gene_cols %in% colnames(gpl_table)][1]

if (is.na(gene_col)) {
  cat("\nPossible gene-related columns:\n")
  print(grep("symbol|gene", colnames(gpl_table), ignore.case = TRUE, value = TRUE))
  stop("Gene symbol column was not found in GPL96 annotation.")
}

cat("\nDetected gene symbol column:", gene_col, "\n")

probe2gene <- gpl_table %>%
  dplyr::select(
    probe_id = ID,
    gene_symbol = all_of(gene_col)
  ) %>%
  dplyr::filter(!is.na(probe_id), probe_id != "") %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
  dplyr::distinct(probe_id, gene_symbol, .keep_all = TRUE)

cat("\nprobe2gene head:\n")
print(head(probe2gene))

############################################################
# 3. Read expression matrix
############################################################

if (!file.exists(expr_file)) {
  stop(paste0("Expression file not found:\n", expr_file))
}

con <- file(expr_file, open = "rt")
lines <- readLines(con, n = 5000)
close(con)

header_line <- grep("!series_matrix_table_begin", lines)

if (length(header_line) == 0) {
  stop("!series_matrix_table_begin was not found. Check expression file format.")
}

cat("\nseries_matrix_table_begin line:", header_line, "\n")

expr_con <- file(expr_file, open = "rt")

expr <- read.table(
  expr_con,
  header = TRUE,
  sep = "\t",
  skip = header_line,
  comment.char = "!",
  quote = "\"",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

close(expr_con)

rownames(expr) <- expr$ID_REF
expr <- expr[, -1, drop = FALSE]

colnames(expr) <- gsub('"', "", colnames(expr))
colnames(expr) <- trimws(colnames(expr))

expr <- as.data.frame(expr)
expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))
expr <- as.matrix(expr)

cat("\nExpression matrix dimension:\n")
print(dim(expr))

cat("\nExpression sample IDs head:\n")
print(head(colnames(expr), 20))

############################################################
# 4. Match phenotype and expression samples
############################################################

cat("\nPhenotype ID-related columns:\n")
print(grep("geo|accession|sample|gsm", colnames(pheno), ignore.case = TRUE, value = TRUE))

if ("geo_accession" %in% colnames(pheno)) {
  pheno$sample_match_id <- as.character(pheno$geo_accession)
} else if ("sample_id" %in% colnames(pheno)) {
  pheno$sample_match_id <- as.character(pheno$sample_id)
} else {
  pheno$sample_match_id <- as.character(rownames(pheno))
}

pheno$sample_match_id <- gsub('"', "", pheno$sample_match_id)
pheno$sample_match_id <- trimws(pheno$sample_match_id)

cat("\nPhenotype sample_match_id head:\n")
print(head(pheno$sample_match_id, 20))

valid_samples <- intersect(colnames(expr), pheno$sample_match_id)

cat("\nMatched sample number:", length(valid_samples), "\n")
print(valid_samples)

if (length(valid_samples) == 0) {
  stop("No matched sample IDs between phenotype and expression matrix.")
}

# expression column order 기준으로 phenotype 정렬
expr2 <- expr[, valid_samples, drop = FALSE]
pheno2 <- pheno[match(valid_samples, pheno$sample_match_id), , drop = FALSE]

if (!all(colnames(expr2) == pheno2$sample_match_id)) {
  stop("Sample order mismatch after matching.")
}

cat("\nSample matching confirmed:\n")
print(all(colnames(expr2) == pheno2$sample_match_id))

############################################################
# 5. Extract region and disease
############################################################

pheno2$region <- ifelse(
  grepl("^Superior frontal gyrus", pheno2$title),
  "Superior frontal gyrus",
  ifelse(
    grepl("^Lateral substantia nigra", pheno2$title),
    "Lateral substantia nigra",
    ifelse(
      grepl("^Medial substantia nigra", pheno2$title),
      "Medial substantia nigra",
      NA
    )
  )
)

pheno2$disease <- ifelse(
  grepl("Parkinson's disease", pheno2$title),
  "Parkinson's disease",
  ifelse(
    grepl("control", pheno2$title, ignore.case = TRUE),
    "control",
    NA
  )
)

cat("\nRegion table:\n")
print(table(pheno2$region, useNA = "ifany"))

cat("\nDisease table:\n")
print(table(pheno2$disease, useNA = "ifany"))

cat("\nRegion x Disease table:\n")
print(table(pheno2$region, pheno2$disease, useNA = "ifany"))

# Export exact inputs for the separate GSEA script without changing DEG results.
if (!all(c("sample_match_id", "region", "disease") %in% names(pheno2))) {
  stop("GSE8397 phenotype parsing incomplete; cannot export GSEA inputs.")
}
saveRDS(list(expr=expr2, pheno=pheno2, probe2gene=probe2gene),
        file.path(out_dir, "GSE8397_GSEA_input_objects.rds"))
cat("Prepared GSEA inputs saved in: ", out_dir, "\n")

############################################################
# 6. Define KEGG hsa04750 TRPV1-related genes
############################################################

kegg_04750_genes <- c(
  "MAPK9",
  "PLCB4",
  "MAPK8",
  "MAPK10",
  "MAPK14",
  "IL1R1",
  "PLCG1",
  "IL1B",
  "NGF",
  "PLCB2",
  "PLCB3",
  "MAPK13",
  "PLCB1",
  "MAPK11",
  "MAPK12",
  "IL1RAP",
  "TRPV1",
  "PLCG2"
)

############################################################
# 7. Analysis by region
############################################################

regions <- unique(na.omit(pheno2$region))

all_kegg_best_results <- list()
all_kegg_label_results <- list()
all_trpv1_check_results <- list()

for (reg in regions) {
  
  cat("\n==================================================\n")
  cat("Region:", reg, "\n")
  cat("==================================================\n")
  
  idx <- which(
    pheno2$region == reg &
      pheno2$disease %in% c("control", "Parkinson's disease")
  )
  
  pheno_reg <- pheno2[idx, , drop = FALSE]
  expr_reg <- expr2[, idx, drop = FALSE]
  
  if (!all(colnames(expr_reg) == pheno_reg$sample_match_id)) {
    stop(paste0("Sample order mismatch in region: ", reg))
  }
  
  group <- factor(
    pheno_reg$disease,
    levels = c("control", "Parkinson's disease")
  )
  
  cat("\nSample number:\n")
  print(table(group))
  
  cat("\nSamples used in this region:\n")
  print(data.frame(
    sample = pheno_reg$sample_match_id,
    title = pheno_reg$title,
    region = pheno_reg$region,
    disease = pheno_reg$disease
  ))
  
  if (length(unique(group)) < 2 || length(group) < 4) {
    cat("Skipped: not enough samples or only one group.\n")
    next
  }
  
  ##########################################################
  # 7-1. Limma DEG analysis
  # Explicit contrast: PD - Control
  ##########################################################
  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("Control", "PD")
  
  cat("\nDesign matrix columns:\n")
  print(colnames(design))
  
  contrast_matrix <- limma::makeContrasts(
    PD_vs_Control = PD - Control,
    levels = design
  )
  
  fit <- limma::lmFit(expr_reg, design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  
  tt <- limma::topTable(
    fit2,
    coef = "PD_vs_Control",
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )
  
  tt$probe <- rownames(tt)
  tt$negLogP <- -log10(tt$P.Value)
  
  tt$color <- ifelse(
    tt$P.Value < 0.05 & tt$logFC > 0,
    "#FA8072",
    ifelse(
      tt$P.Value < 0.05 & tt$logFC < 0,
      "#648FFF",
      "grey"
    )
  )
  
  ##########################################################
  # 7-2. Merge annotation
  ##########################################################
  
  tt_annot <- tt %>%
    dplyr::left_join(probe2gene, by = c("probe" = "probe_id"))
  
  ##########################################################
  # 7-3. KEGG gene filtering
  ##########################################################
  
  kegg_long <- tt_annot %>%
    dplyr::select(
      probe,
      gene_symbol,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B,
      negLogP,
      color
    ) %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::mutate(
      gene_symbol_sep = strsplit(as.character(gene_symbol), " /// ")
    ) %>%
    tidyr::unnest(gene_symbol_sep) %>%
    dplyr::mutate(
      gene_symbol_sep = trimws(as.character(gene_symbol_sep))
    ) %>%
    dplyr::filter(gene_symbol_sep %in% kegg_04750_genes)
  
  ##########################################################
  # 7-4. Select one best probe per gene
  # Rule:
  # If multiple probes map to the same gene,
  # keep the probe with the lowest P.Value.
  ##########################################################
  
  kegg_top <- kegg_long %>%
    dplyr::group_by(gene_symbol_sep) %>%
    dplyr::arrange(P.Value, dplyr::desc(abs(logFC)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(gene_symbol_sep, .keep_all = TRUE)
  
  ##########################################################
  # 7-5. Label only genes with P.Value < 0.05
  ##########################################################
  
  label_df <- kegg_top %>%
    dplyr::filter(P.Value < 0.05) %>%
    dplyr::mutate(
      label_gene = gene_symbol_sep,
      lab_col = dplyr::case_when(
        logFC > 0 ~ "#FA8072",
        logFC < 0 ~ "#648FFF",
        TRUE ~ "grey"
      )
    )
  
  cat("\nKEGG best probe per gene:\n")
  print(
    kegg_top %>%
      dplyr::select(gene_symbol_sep, probe, gene_symbol, logFC, P.Value, adj.P.Val)
  )
  
  cat("\nGenes shown on volcano plot, P.Value < 0.05 only:\n")
  print(
    label_df %>%
      dplyr::select(label_gene, probe, gene_symbol, logFC, P.Value, adj.P.Val)
  )
  
  ##########################################################
  # 7-6. TRPV1 direction check
  ##########################################################
  
  trpv1_all <- kegg_long %>%
    dplyr::filter(gene_symbol_sep == "TRPV1") %>%
    dplyr::arrange(P.Value)
  
  if (nrow(trpv1_all) > 0) {
    
    trpv1_check <- lapply(trpv1_all$probe, function(pb) {
      
      expr_probe <- as.numeric(expr_reg[pb, ])
      
      control_values <- expr_probe[group == "control"]
      pd_values <- expr_probe[group == "Parkinson's disease"]
      
      control_mean <- mean(control_values, na.rm = TRUE)
      pd_mean <- mean(pd_values, na.rm = TRUE)
      pd_minus_control <- pd_mean - control_mean
      
      limma_logfc <- trpv1_all$logFC[trpv1_all$probe == pb][1]
      pval <- trpv1_all$P.Value[trpv1_all$probe == pb][1]
      adjp <- trpv1_all$adj.P.Val[trpv1_all$probe == pb][1]
      raw_symbol <- trpv1_all$gene_symbol[trpv1_all$probe == pb][1]
      
      data.frame(
        region = reg,
        probe = pb,
        gene_symbol = raw_symbol,
        control_n = length(control_values),
        pd_n = length(pd_values),
        control_mean = control_mean,
        pd_mean = pd_mean,
        pd_minus_control = pd_minus_control,
        limma_logFC_PD_vs_Control = limma_logfc,
        P.Value = pval,
        adj.P.Val = adjp,
        control_values = paste(round(control_values, 4), collapse = ", "),
        pd_values = paste(round(pd_values, 4), collapse = ", "),
        direction_by_mean = ifelse(
          pd_minus_control > 0,
          "PD higher",
          ifelse(pd_minus_control < 0, "PD lower", "same")
        ),
        direction_by_limma = ifelse(
          limma_logfc > 0,
          "PD higher",
          ifelse(limma_logfc < 0, "PD lower", "same")
        ),
        direction_match = sign(pd_minus_control) == sign(limma_logfc)
      )
    }) %>%
      dplyr::bind_rows()
    
    cat("\n******** TRPV1 DIRECTION CHECK ********\n")
    print(trpv1_check)
    cat("****************************************\n")
    
    all_trpv1_check_results[[reg]] <- trpv1_check
  }
  
  ##########################################################
  # 7-7. Store results
  ##########################################################
  
  all_kegg_best_results[[reg]] <- kegg_top %>%
    dplyr::mutate(region = reg)
  
  all_kegg_label_results[[reg]] <- label_df %>%
    dplyr::mutate(region = reg)
  
  ##########################################################
  # 7-8. Volcano plot
  ##########################################################
  
  p <- ggplot(tt_annot, aes(x = logFC, y = negLogP)) +
    
    geom_point(
      aes(color = color),
      size = 2,
      alpha = 0.75
    ) +
    
    geom_point(
      data = label_df,
      aes(x = logFC, y = negLogP, fill = lab_col),
      color = "black",
      shape = 21,
      size = 2.6,
      stroke = 0.6,
      alpha = 1,
      show.legend = FALSE
    ) +
    
    scale_color_identity() +
    scale_fill_identity() +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      color = "black"
    ) +
    
    labs(
      title = paste0("Volcano plot: ", reg),
      x = "Log2 Fold Change (PD vs Control)",
      y = "-log10(P value)"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    
    ggrepel::geom_label_repel(
      data = label_df,
      aes(x = logFC, y = negLogP, label = label_gene),
      color = "black",
      size = 4,
      fontface = "bold",
      box.padding = 0.4,
      label.size = 0.15,
      segment.color = "grey40",
      min.segment.length = 0.1,
      label.r = grid::unit(0.12, "lines"),
      label.padding = grid::unit(0.25, "lines"),
      fill = scales::alpha("white", 0.7),
      max.overlaps = Inf,
      force = 6
    )
  
  print(p)
  
  ##########################################################
  # 7-9. Save outputs
  ##########################################################
  
  reg_clean <- gsub("[^A-Za-z0-9]+", "_", reg)
  
  pdf_file <- file.path(
    out_dir,
    paste0("volcano_kegg04750_", reg_clean, "_PD_vs_Control_p005_label.pdf")
  )
  
  png_file <- file.path(
    out_dir,
    paste0("volcano_kegg04750_", reg_clean, "_PD_vs_Control_p005_label.png")
  )
  
  all_deg_csv <- file.path(
    out_dir,
    paste0("all_DEG_", reg_clean, "_PD_vs_Control.csv")
  )
  
  best_probe_csv <- file.path(
    out_dir,
    paste0("kegg04750_best_probe_all_", reg_clean, ".csv")
  )
  
  labeled_csv <- file.path(
    out_dir,
    paste0("kegg04750_labeled_p005_", reg_clean, ".csv")
  )
  
  sample_csv <- file.path(
    out_dir,
    paste0("samples_used_", reg_clean, ".csv")
  )
  
  ggsave(
    filename = pdf_file,
    plot = p,
    width = 7.5,
    height = 5.5
  )
  
  ggsave(
    filename = png_file,
    plot = p,
    width = 7.5,
    height = 5.5,
    dpi = 300
  )
  
  write.csv(
    tt_annot,
    all_deg_csv,
    row.names = FALSE
  )
  
  write.csv(
    kegg_top %>%
      dplyr::select(
        gene_symbol_sep,
        probe,
        gene_symbol,
        logFC,
        AveExpr,
        t,
        P.Value,
        adj.P.Val,
        B
      ),
    best_probe_csv,
    row.names = FALSE
  )
  
  write.csv(
    label_df %>%
      dplyr::select(
        label_gene,
        probe,
        gene_symbol,
        logFC,
        AveExpr,
        t,
        P.Value,
        adj.P.Val,
        B
      ),
    labeled_csv,
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      sample = pheno_reg$sample_match_id,
      title = pheno_reg$title,
      region = pheno_reg$region,
      disease = pheno_reg$disease
    ),
    sample_csv,
    row.names = FALSE
  )
  
  cat("\nSaved files for region:", reg, "\n")
  cat(pdf_file, "\n")
  cat(png_file, "\n")
  cat(all_deg_csv, "\n")
  cat(best_probe_csv, "\n")
  cat(labeled_csv, "\n")
  cat(sample_csv, "\n")
}

############################################################
# 8. Save combined results
############################################################

if (length(all_kegg_best_results) > 0) {
  
  all_kegg_best_df <- dplyr::bind_rows(all_kegg_best_results)
  
  write.csv(
    all_kegg_best_df %>%
      dplyr::select(
        region,
        gene_symbol_sep,
        probe,
        gene_symbol,
        logFC,
        AveExpr,
        t,
        P.Value,
        adj.P.Val,
        B
      ),
    file.path(out_dir, "GSE8397_all_regions_kegg04750_best_probe_all.csv"),
    row.names = FALSE
  )
}

if (length(all_kegg_label_results) > 0) {
  
  all_kegg_label_df <- dplyr::bind_rows(all_kegg_label_results)
  
  write.csv(
    all_kegg_label_df %>%
      dplyr::select(
        region,
        label_gene,
        probe,
        gene_symbol,
        logFC,
        AveExpr,
        t,
        P.Value,
        adj.P.Val,
        B
      ),
    file.path(out_dir, "GSE8397_all_regions_kegg04750_labeled_p005.csv"),
    row.names = FALSE
  )
}

if (length(all_trpv1_check_results) > 0) {
  
  all_trpv1_check_df <- dplyr::bind_rows(all_trpv1_check_results)
  
  write.csv(
    all_trpv1_check_df,
    file.path(out_dir, "GSE8397_all_regions_TRPV1_probe_direction_check.csv"),
    row.names = FALSE
  )
}

cat("\nDone.\n")
cat("Output directory:\n")
cat(out_dir, "\n")
