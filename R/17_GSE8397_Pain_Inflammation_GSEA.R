# Input/output paths are configured via environment variables: PD_TRPV1_ROOT, PD_TRPV1_GSE8397_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)
if (!nzchar(Sys.getenv("PD_TRPV1_GSE8397_DIR", unset=""))) stop("Set PD_TRPV1_GSE8397_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 17_GSE8397_Pain_Inflammation_GSEA.R
# Requires original project datasets; full R execution not validated here.
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(dplyr)
library(tibble)
library(fgsea)
library(purrr)
library(tidyr)
library(limma)

# Read exact matched samples/region/disease/probe annotation from 16.
gse_input <- file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'GSE8397_volcano_final', 'GSE8397_GSEA_input_objects.rds')
if (!file.exists(gse_input)) stop("Run 16 first; missing: ", gse_input)
gse_objects <- readRDS(gse_input)
expr <- gse_objects$expr
pheno <- gse_objects$pheno
probe2gene <- gse_objects$probe2gene
if (!all(c("sample_match_id", "region", "disease") %in% names(pheno))) stop("GSE8397 phenotype fields missing")
pheno$geo_accession <- trimws(pheno$sample_match_id)
colnames(expr) <- trimws(colnames(expr))
if (!identical(colnames(expr), pheno$geo_accession)) stop("GSE8397 expression and phenotype sample order mismatch")

# 1. Pain/Inflammation 관련 GO term 추출
pain_keywords <- c(
  "pain", "nociception", "hyperalgesia", "allodynia",
  "inflammation", "inflammatory", "injury", "immune",
  "cytokine", "neuroinflammation", "edema", "swelling", "hypersensitivity"
)
all_go <- AnnotationDbi::toTable(GO.db::GOTERM)
if (!"go_id" %in% names(all_go) && "GOID" %in% names(all_go)) names(all_go)[names(all_go)=="GOID"] <- "go_id"
if (!all(c("go_id", "Term", "Ontology") %in% names(all_go))) stop("Unexpected GO.db annotation columns: ", paste(names(all_go),collapse=", "))
all_go <- all_go %>%
  distinct(go_id, .keep_all = TRUE) %>%
  filter(
    !is.na(go_id), !is.na(Term), !is.na(Ontology),
    grepl(paste(pain_keywords, collapse = "|"), Term, ignore.case = TRUE)
  )
pain_go_ids <- all_go$go_id

# 2. Pain 관련 GO term의 gene set 리스트 생성
go2gene <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = pain_go_ids,
  columns = c("SYMBOL"),
  keytype = "GOALL"
)
go2gene <- go2gene[!is.na(go2gene$GOALL) & !is.na(go2gene$SYMBOL), ]
pain_gmt <- split(go2gene$SYMBOL, go2gene$GOALL)

# 3. probe2gene에서 probe_id ~ gene_symbol 1:1 매핑 (/// 분리)
probe2gene_split <- probe2gene %>%
  mutate(gene_symbol = strsplit(as.character(gene_symbol), " /// ")) %>%
  tidyr::unnest(gene_symbol)

# 4. region 리스트 자동 추출(모든 brain region)
region_list <- sort(unique(pheno$region))
gsea_results <- list()

for (reg in region_list) {
  region_samples <- pheno$geo_accession[pheno$region == reg & pheno$disease %in% c("control", "Parkinson's disease")]
  region_samples <- intersect(region_samples, colnames(expr))
  region_pheno <- pheno[match(region_samples, pheno$geo_accession), ]
  group_counts <- table(region_pheno$disease)
  
  cat("\n======", reg, "======\n")
  print(group_counts)
  
  # 분석 조건: 두 그룹 모두 2개 이상일 때만 진행
  if (all(c("control", "Parkinson's disease") %in% names(group_counts)) && all(group_counts[c("control", "Parkinson's disease")] >= 2)) {
    cat(sprintf("region %s: 분석 진행\n", reg))
    region_expr <- expr[, region_samples, drop = FALSE]
    group <- factor(region_pheno$disease, levels = c("control", "Parkinson's disease"))
    design <- model.matrix(~ group)
    fit <- limma::lmFit(region_expr, design)
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    tt$probe <- rownames(tt)
    tt <- left_join(
      tt,
      probe2gene_split[, c("probe_id", "gene_symbol")],   # 여기서 [, c(...)]처럼 base R로!
      by = c("probe" = "probe_id")
    )
    tt <- tt[!is.na(tt$gene_symbol) & tt$gene_symbol != "" & !is.na(tt$logFC), ]
    tt <- tt[order(-abs(tt$logFC)), ]
    tt <- tt[!duplicated(tt$gene_symbol), ]
    ranks <- tt$logFC
    names(ranks) <- tt$gene_symbol
    ranks <- sort(ranks, decreasing = TRUE)
    gsea <- fgsea(
      pathways = pain_gmt,
      stats = ranks,
      minSize = 5, maxSize = 500
    )
    gsea$region <- reg
    if (nrow(gsea) > 0L) gsea_results[[reg]] <- gsea
  } else {
    cat(sprintf("region %s: 샘플 수 부족, 분석 생략\n", reg))
    next
  }
}

# 5. 결과 통합 및 주요 컬럼만 표로 정리
if (!length(gsea_results)) stop("No GSEA results; inspect region grouping, GO gene mapping and rank identifiers")
all_gsea <- bind_rows(gsea_results)
final_res <- dplyr::select(all_gsea, region, pathway, NES, pval, padj, size, leadingEdge)
cat("\n==== GSEA 결과 미리보기 ====\n")
print(head(final_res))
final_res$leadingEdge <- map_chr(final_res$leadingEdge, ~ paste(.x, collapse = ","))
gse_out <- Sys.getenv("PD_TRPV1_GSE8397_DIR", unset="")
dir.create(gse_out, recursive=TRUE, showWarnings=FALSE)
write.csv(final_res, file.path(gse_out, "pain_inflam_GO_GSEA_by_region.csv"), row.names = FALSE)
