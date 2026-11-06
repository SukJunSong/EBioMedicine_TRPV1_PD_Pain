# Input/output paths are configured via environment variables: PD_TRPV1_ROOT, PD_TRPV1_CLINICAL_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)
if (!nzchar(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""))) stop("Set PD_TRPV1_CLINICAL_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 06_PPMI_V02_DESeq2_Age_Sex_Adjusted.R
# Requires original project datasets; full R execution not validated here.
# V02 DESeq2: SAME samples/counts as saved original, add age and sex to design
# Original: ~ group; Adjusted: ~ age + sex + group
# Run entire script after setting PD_TRPV1_ROOT and PD_TRPV1_CLINICAL_DIR.
# Dependencies: DESeq2, readxl, readr
# Input: original A_original_dds.rds from previous run, V02 metadata, age-at-visit CSV.
# Duplicate V02 AGE_AT_VISIT conflict resolution: median of finite values + CSV audit.
# Do not modify original RDS or previously generated figures.

run_adjusted_V02 <- function() {
  suppressPackageStartupMessages({library(DESeq2); library(readxl); library(readr)})
  root <- Sys.getenv("PD_TRPV1_ROOT", unset="")
  prior <- file.path(root,'V02_DESeq2_Original_vs_Medication')
  out <- file.path(root,'V02_DESeq2_AgeSex_Adjusted')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)

  choose <- function(paths,what) {
    hit <- paths[file.exists(paths)]
    if(!length(hit)) stop('Cannot find ',what,'; looked in:\n',paste(paths,collapse='\n'))
    hit[1]
  }
  dds_file <- choose(c(file.path(prior,'A_original_dds.rds'),file.path(root,'A_original_dds.rds')), 'original DESeq2 RDS')
  meta_file <- choose(c(file.path(root,'PPMI','metaData V02.xlsx'),file.path(root,'metaData V02.xlsx'),
                        file.path(root,'PPMI','metaData BL.xlsx'),file.path(root,'metaData BL.xlsx')),'sex metadata xlsx')
  age_file <- choose(c(file.path(root,'Age_at_visit_20Mar2025.csv'),
                       file.path(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""), 'Age_at_visit_20Mar2025.csv'),
                       file.path(root,'PPMI','Age_at_visit_20Mar2025.csv')),'Age_at_visit_20Mar2025.csv')

  original <- readRDS(dds_file)
  if(!inherits(original,'DESeqDataSet')) stop('RDS is not a DESeqDataSet: ',dds_file)
  if(!'group' %in% names(SummarizedExperiment::colData(original))) stop('Original RDS lacks group field')
  ids <- colnames(original)
  if(anyDuplicated(ids)) stop('Original RDS has duplicated PATNO column names')
  original_group <- factor(as.character(original$group),
    levels=c('Control','No pain PD','Low pain PD','High pain PD'))
  if(anyNA(original_group)) stop('Unexpected group labels in original RDS')

  key <- function(x) trimws(as.character(x))
  meta <- as.data.frame(readxl::read_excel(meta_file))
  if(!all(c('PATNO','GENDER') %in% names(meta))) stop('metadata must contain PATNO and GENDER; columns: ',paste(names(meta),collapse=', '))
  meta$PATNO <- key(meta$PATNO)
  sex <- meta[,c('PATNO','GENDER')]
  sex <- sex[!is.na(sex$PATNO) & nzchar(sex$PATNO),,drop=FALSE]
  # Duplicates are allowed only if they have identical known sex.
  badsex <- vapply(split(as.character(sex$GENDER),sex$PATNO),function(v)length(unique(v[!is.na(v)&nzchar(v)]))>1,logical(1))
  if(any(badsex)) stop('Conflicting sex for PATNO: ',paste(names(badsex)[badsex],collapse=', '))
  sex <- sex[!duplicated(sex$PATNO),,drop=FALSE]
  sex$GENDER <- tolower(trimws(as.character(sex$GENDER)))
  sex$GENDER[sex$GENDER %in% c('m','1')] <- 'male'
  sex$GENDER[sex$GENDER %in% c('f','2')] <- 'female'

  ages <- as.data.frame(readr::read_csv(age_file,show_col_types=FALSE,
       col_types=readr::cols(PATNO=readr::col_character(), EVENT_ID=readr::col_character(),
                              AGE_AT_VISIT=readr::col_double(),.default=readr::col_guess())))
  if(!all(c('PATNO','EVENT_ID','AGE_AT_VISIT') %in% names(ages))) stop('Age file needs PATNO, EVENT_ID, AGE_AT_VISIT')
  ages$PATNO <- key(ages$PATNO)
  ages$EVENT_ID <- trimws(as.character(ages$EVENT_ID))
  ages <- ages[ages$EVENT_ID=='V02',c('PATNO','AGE_AT_VISIT'),drop=FALSE]
  # PPMI age-at-visit source occasionally contains TWO different ages at the SAME V02 visit.
  # Keep an explicit audit trail. Resolve finite duplicate values using the within-V02 median;
  # do not arbitrarily take the first row, and do not drop these original participants.
  ages <- ages[!is.na(ages$PATNO) & nzchar(ages$PATNO),,drop=FALSE]
  age_splits <- split(ages$AGE_AT_VISIT,ages$PATNO)
  age_qc <- do.call(rbind,lapply(names(age_splits),function(id) {
    v <- age_splits[[id]]; f <- v[is.finite(v)]
    data.frame(PATNO=id, n_records=length(v), n_finite=length(f),
       distinct_age_values=paste(sort(unique(f)),collapse='; '),
       age_range=if(length(f)) max(f)-min(f) else NA_real_,
       age_used=if(length(f)) median(f) else NA_real_,
       conflict=length(unique(f))>1, stringsAsFactors=FALSE)
  }))
  readr::write_csv(age_qc[age_qc$n_records>1 | age_qc$conflict,,drop=FALSE],
                   file.path(out,'00b_V02_duplicate_age_audit.csv'))
  if(any(age_qc$conflict)) {
    message('V02 duplicate ages resolved by within-visit median for ',sum(age_qc$conflict),
            ' PATNOs. Review 00b_V02_duplicate_age_audit.csv.')
  }
  ages <- data.frame(PATNO=age_qc$PATNO, AGE_AT_VISIT=age_qc$age_used,
                     stringsAsFactors=FALSE)

  info <- data.frame(PATNO=ids, group=original_group,
       age=ages$AGE_AT_VISIT[match(ids,ages$PATNO)],
       sex=sex$GENDER[match(ids,sex$PATNO)],stringsAsFactors=FALSE)
  readr::write_csv(info,file.path(out,'00_original_cohort_age_sex_QC.csv'),na='')
  missing <- is.na(info$age) | !is.finite(info$age) | is.na(info$sex) | !(info$sex %in% c('female','male'))
  if(any(missing)) stop(sum(missing),' of ',length(ids),' samples have missing/invalid V02 age or sex. ',
      'Inspect ',file.path(out,'00_original_cohort_age_sex_QC.csv'),
      '. Stopping to PRESERVE original cohort; no silent sample exclusion.')
  if(any(info$age<18 | info$age>110)) stop('Unexpected ages; inspect QC CSV')

  info$sex <- factor(info$sex,levels=c('female','male'))
  info$group <- factor(info$group,levels=c('Control','No pain PD','Low pain PD','High pain PD'))
  rownames(info) <- info$PATNO
  stopifnot(identical(ids,rownames(info)))
  X <- model.matrix(~age+sex+group,info)
  if(qr(X)$rank != ncol(X)) stop('Rank-deficient design: age/sex/group cannot be jointly estimated')
  print(table(info$group))
  print(table(info$group,info$sex))
  message('Same ',length(ids),' subjects as original; no changes to counts or groups.')

  # Fresh DESeqDataSet from stored exact integer count matrix for reproducibility.
  counts <- DESeq2::counts(original,normalized=FALSE)
  if(any(counts<0) || any(abs(counts-round(counts))>1e-8)) stop('Noninteger original counts')
  adjusted <- DESeq2::DESeqDataSetFromMatrix(countData=counts,colData=info,
                                             design=~age+sex+group)
  adjusted <- adjusted[rowSums(DESeq2::counts(adjusted))>0,]
  adjusted <- DESeq2::DESeq(adjusted)
  saveRDS(adjusted,file.path(out,'B_age_sex_adjusted_dds.rds'))

  # Extract all HC-reference contrasts + optional original model re-fit on same samples.
  extract <- function(dds,comparison) {
    rr <- DESeq2::results(dds,contrast=c('group',comparison,'Control'),alpha=.05)
    z <- as.data.frame(rr)
    z$gene_id <- rownames(z)
    z$comparison <- paste(comparison,'vs Control')
    z[,c('gene_id','comparison','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
  }
  targets <- c('No pain PD','Low pain PD','High pain PD')
  adjusted_results <- do.call(rbind,lapply(targets,function(g)extract(adjusted,g)))
  readr::write_csv(adjusted_results,file.path(out,'01_DESeq2_AgeSex_adjusted_ALL_genes_vs_HC.csv'),na='')
  original_results <- do.call(rbind,lapply(targets,function(g)extract(original,g)))
  readr::write_csv(original_results,file.path(out,'02_DESeq2_ORIGINAL_ALL_genes_vs_HC.csv'),na='')
  merged <- merge(original_results,adjusted_results,by=c('gene_id','comparison'),
                  suffixes=c('_original','_ageSex'),sort=FALSE)
  readr::write_csv(merged,file.path(out,'03_Original_vs_AgeSex_ALL_genes.csv'),na='')
  group_n <- as.data.frame(table(info$group,info$sex))
  names(group_n) <- c('group','sex','n')
  readr::write_csv(group_n,file.path(out,'04_Group_by_Sex_N.csv'))

  # IDs are those from .sf; could be ENSEMBL, no unverified symbol mapping.
  interesting <- c('ENSG00000196689')  # TRPV1 canonical Ensembl ID, only if IDs match
  out_trpv1 <- merged[grepl('^ENSG00000196689(\\.|$)',merged$gene_id),,drop=FALSE]
  readr::write_csv(out_trpv1,file.path(out,'05_TRPV1_original_vs_ageSex.csv'),na='')
  if(!nrow(out_trpv1)) message('TRPV1 row not found by Ensembl ID; inspect gene_id naming in 03 CSV.')
  summary <- aggregate(cbind(n_original=as.integer(!is.na(merged$padj_original)&merged$padj_original<.05),
                            n_adjusted=as.integer(!is.na(merged$padj_ageSex)&merged$padj_ageSex<.05)),
                       by=list(comparison=merged$comparison),FUN=sum)
  readr::write_csv(summary,file.path(out,'06_FDR_significant_count_original_vs_ageSex.csv'))
  message('Done: ',out)
  invisible(list(dds=adjusted,comparison=merged,TRPV1=out_trpv1))
}
run_adjusted_V02()
