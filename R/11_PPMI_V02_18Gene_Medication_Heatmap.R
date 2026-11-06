# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# Supplementary Figure S2: original vs medication-exclusion sensitivity analysis.
# 4 columns: original, analgesic possible, analgesic strict, PD medication possible.
# Uses saved script-05 original DESeq2 (group-only) raw counts and its V02 medication flags.
# For sensitivity cohorts: retain HC as reference, retain PD ONLY when V02 exposure
# status is known and 0, and refit DESeq2 group-only. Never label unknown as unexposed.
# Numbers = high-pain PD vs HC log2FC; stars = DESeq2 padj (BH FDR).

run_supplementary_s2 <- function() {
  pkgs <- c('DESeq2','SummarizedExperiment','readr','dplyr','tidyr',
            'ggplot2','AnnotationDbi','org.Hs.eg.db')
  missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop('Missing packages: ',paste(missing,collapse=', '),call.=FALSE)
  root <- Sys.getenv("PD_TRPV1_ROOT", unset="")
  src <- file.path(root,'V02_DESeq2_Original_vs_Medication')
  ddsfile <- file.path(src,'A_original_dds.rds')
  flagsfile <- file.path(src,'05_V02_medication_patient_flags_ALL_groups.csv')
  out <- file.path(root,'Supplementary_Figure_S2_Medication_18Gene')
  for(f in c(ddsfile,flagsfile)) if(!file.exists(f)) stop('Run 05 first; missing input: ',f,call.=FALSE)

  ids <- c('ENSG00000050748','ENSG00000101333','ENSG00000107643',
    'ENSG00000109339','ENSG00000112062','ENSG00000115594',
    'ENSG00000124181','ENSG00000125538','ENSG00000134259',
    'ENSG00000137841','ENSG00000149782','ENSG00000156711',
    'ENSG00000182621','ENSG00000185386','ENSG00000188130',
    'ENSG00000196083','ENSG00000196689','ENSG00000197943')
  gene_order <- c('NGF','PLCB4','MAPK10','MAPK12','TRPV1','MAPK13','IL1R1',
    'MAPK11','IL1RAP','PLCB3','PLCG2','PLCB1','IL1B','MAPK14','MAPK9',
    'PLCB2','MAPK8','PLCG1')
  gm <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=ids,
    keytype='ENSEMBL',column='SYMBOL',multiVals='first')
  gm <- as.character(gm);names(gm)<-ids
  if(anyNA(gm) || anyDuplicated(unname(gm)) || !setequal(unname(gm),gene_order))
    stop('18-gene annotation does not match expected target set; inspect mappings.',call.=FALSE)
  if(!identical(unname(gm['ENSG00000196689']),'TRPV1')) stop('TRPV1 ID mismatch')
  id_by_gene <- setNames(names(gm),unname(gm))

  original <- readRDS(ddsfile)
  if(!methods::is(original,'DESeqDataSet')) stop('Invalid original DESeq2 RDS')
  if(!'group' %in% names(SummarizedExperiment::colData(original))) stop('Original DDS lacks group')
  if (!identical(all.vars(DESeq2::design(original)), 'group')) {
    stop('Original DDS design must be ~group: ',
         paste(deparse(DESeq2::design(original)),collapse=' '))
  }
  sample_ids <- colnames(original)
  if(anyDuplicated(sample_ids)) stop('Duplicate PATNO in original DDS')
  gr <- as.character(original$group)
  if(anyNA(gr)) stop('Missing group in original DDS')
  lev <- c('Control','No pain PD','Low pain PD','High pain PD')
  if(!all(lev %in% gr)) stop('Original DDS missing one of the four groups')

  flags <- readr::read_csv(flagsfile,show_col_types=FALSE,col_types=readr::cols(.default='c'))
  flagvars <- c('analgesic_possible','analgesic_strict','PD_drug_possible')
  if(!all(c('PATNO','V02_month',flagvars) %in% names(flags)))
    stop('05 medication flags missing columns: ',paste(setdiff(c('PATNO','V02_month',flagvars),names(flags)),collapse=', '))
  flags$PATNO <- trimws(as.character(flags$PATNO))
  if(anyDuplicated(flags$PATNO)) stop('Duplicate PATNO in medication flags CSV')
  flagged <- flags[match(sample_ids,flags$PATNO),,drop=FALSE]
  if(!identical(flagged$PATNO,sample_ids)) {
    stop('PATNO mismatch between original DDS and 05 medication flags; do not silently match or drop')
  }
  for(k in flagvars) {
    flagged[[k]] <- suppressWarnings(as.integer(flagged[[k]]))
    if(any(!is.na(flagged[[k]]) & !(flagged[[k]] %in% 0:1))) stop('Nonbinary exposure flag: ',k)
  }
  # Absent/unresolved V02 visit dates must NEVER be assumed unexposed in PD.
  visit_known <- !is.na(flagged$V02_month) & nzchar(trimws(flagged$V02_month))
  names_scenarios <- c('Original','Exclude analgesic (possible)',
                       'Exclude analgesic (strict)','Exclude PD drugs (possible)')
  masks <- list(rep(TRUE,length(sample_ids)))
  for(k in flagvars) {
    # Retain all HC. For PD, require known visit month + explicit zero flag.
    masks[[length(masks)+1L]] <- gr=='Control' |
      (gr!='Control' & visit_known & !is.na(flagged[[k]]) & flagged[[k]]==0L)
  }
  names(masks) <- names_scenarios
  dir.create(out,recursive=TRUE,showWarnings=FALSE)

  cnt <- DESeq2::counts(original,normalized=FALSE)
  if(anyNA(cnt) || any(cnt<0) || any(cnt!=round(cnt))) stop('Original counts not valid integers')
  raw_ids <- sub('\\.[0-9]+$','',rownames(original))
  if(anyDuplicated(raw_ids)) stop('Duplicated version-stripped Ensembl IDs in original DDS')
  results_list <- list(); qc_list <- list()
  for(i in seq_along(masks)) {
    nm <- names(masks)[i]
    keep <- masks[[i]]
    counts_group <- table(factor(gr[keep],levels=lev))
    qc_list[[i]] <- data.frame(Scenario=nm,Group=lev,n=as.integer(counts_group),
      pd_unknown_visit_n=sum(gr!='Control' & !visit_known),
      stringsAsFactors=FALSE)
    if(any(counts_group<2L)) stop('Insufficient group counts for ',nm,': ',paste(names(counts_group),counts_group,collapse=', '))
    message('SCENARIO ',nm,' group N: ',paste(names(counts_group),counts_group,collapse='; '))
    if(i==1L) {
      fit <- original   # already fitted in 05, avoid changing original results
    } else {
      cd <- data.frame(group=factor(gr[keep],levels=lev),row.names=sample_ids[keep])
      xx <- cnt[,keep,drop=FALSE]
      fit <- DESeq2::DESeqDataSetFromMatrix(countData=xx,colData=cd,design=~group)
      fit <- fit[rowSums(DESeq2::counts(fit))>0,,drop=FALSE]
      fit <- DESeq2::DESeq(fit,quiet=TRUE)
    }
    fit_ids <- sub('\\.[0-9]+$','',rownames(fit))
    rr <- as.data.frame(DESeq2::results(fit,contrast=c('group','High pain PD','Control'),alpha=.05))
    rr$Ensembl <- fit_ids
    if(anyDuplicated(rr$Ensembl)) stop('Duplicate Ensembl after stripping in ',nm)
    rr <- rr[rr$Ensembl %in% ids,c('Ensembl','log2FoldChange','pvalue','padj'),drop=FALSE]
    rr$Scenario <- nm
    results_list[[i]] <- rr
    readr::write_csv(rr,file.path(out,sprintf('%02d_%s_high_vs_HC_18genes.csv',i,gsub('[^a-z0-9]+','_',tolower(nm)))))
  }
  qc <- dplyr::bind_rows(qc_list)
  values <- dplyr::bind_rows(results_list)
  grid <- expand.grid(Gene=gene_order,Scenario=names_scenarios,stringsAsFactors=FALSE)
  grid$Ensembl <- unname(id_by_gene[grid$Gene])
  values <- dplyr::left_join(grid,values,by=c('Ensembl','Scenario'))
  values$Gene <- factor(values$Gene,levels=rev(gene_order))
  values$Scenario <- factor(values$Scenario,levels=names_scenarios)
  values$stars <- ifelse(is.na(values$padj),'',ifelse(values$padj<.001,'***',ifelse(values$padj<.01,'**',ifelse(values$padj<.05,'*',''))))
  values$cell <- ifelse(is.finite(values$log2FoldChange),paste0(sprintf('%.2f',values$log2FoldChange),values$stars),'NA')
  lim <- suppressWarnings(max(abs(values$log2FoldChange),na.rm=TRUE));if(!is.finite(lim)||lim<1)lim<-1
  p <- ggplot2::ggplot(values,ggplot2::aes(x=Scenario,y=Gene,fill=log2FoldChange))+
    ggplot2::geom_tile(color='grey90',linewidth=.6)+
    ggplot2::geom_text(ggplot2::aes(label=cell),size=2.8,color='black')+
    ggplot2::scale_fill_gradient2(low='#648FFF',mid='white',high='#FA8072',
      midpoint=0,limits=c(-lim,lim),na.value='grey90',name='log2FC')+
    ggplot2::labs(x=NULL,y=NULL,caption='* BH-FDR < 0.05; ** < 0.01; *** < 0.001')+
    ggplot2::theme_minimal(base_size=12)+ggplot2::theme(
      panel.grid=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=20,hjust=1),
      axis.text.y=ggplot2::element_text(face='bold'))
  readr::write_csv(qc,file.path(out,'S2_group_N_by_scenario.csv'))
  readr::write_csv(values,file.path(out,'S2_18genes_all_scenarios_values.csv'))
  ggplot2::ggsave(file.path(out,'Supplementary_Figure_S2_medication_heatmap.png'),p,width=11.5,height=9,dpi=320)
  ggplot2::ggsave(file.path(out,'Supplementary_Figure_S2_medication_heatmap.pdf'),p,width=11.5,height=9)
  message('SUPPLEMENTARY FIGURE S2 COMPLETE: ',out)
  invisible(values)
}
run_supplementary_s2()
