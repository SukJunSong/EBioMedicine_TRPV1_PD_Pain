# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 15_PPMI_HighPain_vs_NoPain.R
# Requires original project datasets; full R execution not validated here.
# PPMI V02: High Pain PD vs No Pain PD (DIRECT contrast)
# Two distinct analyses:
# A. ORIGINAL saved DESeq2 V02 count-based model (no re-fitting)
# B. age/sex-adjusted V02 18-gene composite score (separate outcome)
# Input paths as used in previously uploaded analysis scripts.

root <- Sys.getenv("PD_TRPV1_ROOT", unset="")
original_rds <- file.path(root, "V02_DESeq2_Original_vs_Medication", "A_original_dds.rds")
score_candidates <- c(
  file.path(root, "TRPV1_related_gene_score_BL_V02_long.csv"),
  file.path(root, "V02_original_PRESERVED_medication_comparison", "TRPV1_related_gene_score_BL_V02_long.csv"),
  file.path(root, "PPMI", "TRPV1_related_gene_score_BL_V02_long.csv"),
  file.path(root, "PPMI", "BL_V02_paired_TRPV1_analysis", "TRPV1_related_gene_score_BL_V02_lineplot", "TRPV1_related_gene_score_BL_V02_long.csv")
)
score_file <- score_candidates[file.exists(score_candidates)][1]
out <- file.path(root,"V02_HighPain_vs_NoPain_DIRECT")
dir.create(out, recursive=TRUE,showWarnings=FALSE)

# A. Exactly the original fitted DESeq2 object and size factors/dispersions.
if (!requireNamespace("DESeq2",quietly=TRUE)) stop("Install Bioconductor DESeq2")
if (!requireNamespace("SummarizedExperiment",quietly=TRUE)) stop("Install SummarizedExperiment")
if (!file.exists(original_rds)) stop("Missing original fitted DESeq2 object: ", original_rds,
 "\nThis code intentionally refuses to re-create a different model.")
dds <- readRDS(original_rds)
if (!methods::is(dds,"DESeqDataSet")) stop("RDS is not a DESeqDataSet")
cd <- SummarizedExperiment::colData(dds)
if (!"group" %in% names(cd)) stop("No group in cached colData")
groups <- as.character(cd$group)
if (!all(c("High pain PD","No pain PD") %in% groups)) stop("Groups missing: ",paste(unique(groups),collapse=", "))
if (!("group" %in% all.vars(DESeq2::design(dds)))) stop("Original DESeq2 design doesn't include group")
message("Original DESeq2 design: ",paste(deparse(DESeq2::design(dds)),collapse=""))
message("Original cohort:");print(table(groups))
res <- DESeq2::results(dds,
 contrast=c("group","High pain PD","No pain PD"),alpha=.05)
tab <- as.data.frame(res)
tab$Ensembl_versioned <- rownames(tab)
tab$Ensembl <- sub("\\.[0-9]+$","",tab$Ensembl_versioned)
tab$BH_manual_all_genes <- p.adjust(tab$pvalue,method="BH")
tab <- tab[,c("Ensembl_versioned","Ensembl","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","BH_manual_all_genes")]
utils::write.csv(tab,file.path(out,"01_V02_DESeq2_HighPain_vs_NoPain_ALL_genes.csv"),row.names=FALSE)
trpv1 <- tab[tab$Ensembl=="ENSG00000196689",,drop=FALSE]
if(nrow(trpv1)!=1) stop("TRPV1 Ensembl row missing or duplicated (n=",nrow(trpv1),")")
utils::write.csv(trpv1,file.path(out,"02_V02_DESeq2_HighPain_vs_NoPain_TRPV1.csv"),row.names=FALSE)
print(trpv1)

# 18-gene panel IDs from original hsa04750 analysis
panel_ids <- c("ENSG00000050748","ENSG00000101333","ENSG00000107643",
"ENSG00000109339","ENSG00000112062","ENSG00000115594",
"ENSG00000124181","ENSG00000125538","ENSG00000134259",
"ENSG00000137841","ENSG00000149782","ENSG00000156711",
"ENSG00000182621","ENSG00000185386","ENSG00000188130",
"ENSG00000196083","ENSG00000196689","ENSG00000197943")
panel <- tab[tab$Ensembl %in% panel_ids,,drop=FALSE]
if(length(unique(panel$Ensembl))!=18) warning("Only ",length(unique(panel$Ensembl)),"/18 panel genes found")
# Gene names are optional, from the existing original panel mapping file.
mapping_file <- file.path(root,"V02_original_PRESERVED_medication_comparison",
                           "A_original_V02_hsa04750_18_genes.csv")
if(file.exists(mapping_file)) {
 m <- utils::read.csv(mapping_file,stringsAsFactors=FALSE)
 if(all(c("Ensembl","gene_symbol") %in% names(m))) {
  m <- unique(m[,c("Ensembl","gene_symbol")])
  m <- m[!duplicated(m$Ensembl),,drop=FALSE]
  panel$gene_symbol <- m$gene_symbol[match(panel$Ensembl,m$Ensembl)]
 }
}
utils::write.csv(panel,file.path(out,"03_V02_DESeq2_HighPain_vs_NoPain_18genes.csv"),row.names=FALSE)

# B. Separate V02 composite-score model: No Pain PD as reference
if(is.na(score_file)) {
 message("Composite score CSV not found in known paths. DESeq2 outputs saved successfully.")
} else {
 d <- utils::read.csv(score_file,stringsAsFactors=FALSE)
 req <- c("PATNO","Visit","Final_Group","TRPV1_related_score","ENROLL_AGE","GENDER")
 missing <- setdiff(req,names(d))
 if(length(missing)) stop("Score file missing columns: ",paste(missing,collapse=", "))
 d <- d[d$Visit=="V02" & d$Final_Group %in% c("No Pain PD","High Pain PD"),,drop=FALSE]
 if(anyDuplicated(d$PATNO)) stop("Duplicate V02 PATNO in score file")
 d$Final_Group <- factor(d$Final_Group,levels=c("No Pain PD","High Pain PD"))
 d$GENDER <- factor(d$GENDER)
 d$ENROLL_AGE <- as.numeric(d$ENROLL_AGE)
 d$TRPV1_related_score <- as.numeric(d$TRPV1_related_score)
 d <- d[stats::complete.cases(d[,c("TRPV1_related_score","Final_Group","ENROLL_AGE","GENDER")]),,drop=FALSE]
 if(nrow(d)<5) stop("Not enough complete-case score records")
 model <- stats::lm(TRPV1_related_score ~ Final_Group + ENROLL_AGE + GENDER,data=d)
 sm <- summary(model)$coefficients
 coef_name <- "Final_GroupHigh Pain PD"
 if(!coef_name %in% rownames(sm)) stop("Expected High vs No Pain coefficient not found")
 ci <- stats::confint(model,level=.95)[coef_name,]
 x <- data.frame(contrast="High Pain PD - No Pain PD",n_NoPain=sum(d$Final_Group=="No Pain PD"),
 n_HighPain=sum(d$Final_Group=="High Pain PD"),estimate=unname(sm[coef_name,"Estimate"]),
 SE=unname(sm[coef_name,"Std. Error"]),t=unname(sm[coef_name,"t value"]),
 pvalue=unname(sm[coef_name,"Pr(>|t|)"]),CI_lower=unname(ci[1]),CI_upper=unname(ci[2]),
 adjustment="age + sex",stringsAsFactors=FALSE)
 utils::write.csv(x,file.path(out,"04_V02_composite_HighPain_vs_NoPain_age_sex.csv"),row.names=FALSE)
 print(x)
}
message("DONE. All outputs: ",out)
