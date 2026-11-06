# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# Figure 2D: V02, age/sex-adjusted DESeq2, 18 hsa04750 genes
# Input is the ACTUAL output of script 06. No refitting and no averaging of DESeq2 results.
# Plot numbers = log2FC; stars = DESeq2 padj (BH FDR).

run_figure2d <- function() {
  pkgs <- c('readr','dplyr','tidyr','ggplot2','AnnotationDbi','org.Hs.eg.db')
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop('Missing packages: ',paste(missing,collapse=', '),call.=FALSE)

  root <- Sys.getenv("PD_TRPV1_ROOT", unset="")
  inp <- file.path(root,'V02_DESeq2_AgeSex_Adjusted','01_DESeq2_AgeSex_adjusted_ALL_genes_vs_HC.csv')
  out <- file.path(root,'Figure2D_AgeSex_18Gene_Heatmap')
  if (!file.exists(inp)) stop('Run 06 first. Missing file: ',inp,call.=FALSE)

  # Same 18 Ensembl IDs used in script 05 and in the medication sensitivity analysis.
  target_ids <- c('ENSG00000050748','ENSG00000101333','ENSG00000107643',
    'ENSG00000109339','ENSG00000112062','ENSG00000115594',
    'ENSG00000124181','ENSG00000125538','ENSG00000134259',
    'ENSG00000137841','ENSG00000149782','ENSG00000156711',
    'ENSG00000182621','ENSG00000185386','ENSG00000188130',
    'ENSG00000196083','ENSG00000196689','ENSG00000197943')
  wanted_genes <- c('NGF','PLCB4','MAPK10','MAPK12','TRPV1','MAPK13','IL1R1',
    'MAPK11','IL1RAP','PLCB3','PLCG2','PLCB1','IL1B','MAPK14','MAPK9',
    'PLCB2','MAPK8','PLCG1')
  gene_map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
    keys=target_ids,keytype='ENSEMBL',column='SYMBOL',multiVals='first')
  gene_map <- as.character(gene_map); names(gene_map) <- target_ids
  if (anyNA(gene_map) || anyDuplicated(unname(gene_map))) {
    stop('18-gene annotation missing or duplicated: inspect org.Hs.eg.db mapping.',call.=FALSE)
  }
  if (!identical(unname(gene_map['ENSG00000196689']),'TRPV1')) stop('TRPV1 ID mismatch')
  if (!setequal(unname(gene_map),wanted_genes)) {
    stop('Ensembl IDs do not map to expected 18 target symbols. Do not silently relabel genes. ',
      paste(paste(names(gene_map),gene_map,sep='='),collapse='; '),call.=FALSE)
  }
  gene_id_by_name <- setNames(names(gene_map),unname(gene_map))

  raw <- readr::read_csv(inp,show_col_types=FALSE)
  needed <- c('gene_id','comparison','log2FoldChange','pvalue','padj')
  if (!all(needed %in% names(raw))) stop('06 output missing columns: ',paste(setdiff(needed,names(raw)),collapse=', '))
  comparisons <- c('No pain PD vs Control','Low pain PD vs Control','High pain PD vs Control')
  comp_labels <- c('No-pain PD','Low-pain PD','High-pain PD')
  raw$Ensembl <- sub('\\.[0-9]+$','',as.character(raw$gene_id))
  raw$comparison <- as.character(raw$comparison)
  if (!all(comparisons %in% raw$comparison)) stop('Required contrasts absent: ',paste(setdiff(comparisons,raw$comparison),collapse=', '))
  tab <- raw[raw$Ensembl %in% target_ids & raw$comparison %in% comparisons,,drop=FALSE]
  if (anyDuplicated(tab[,c('Ensembl','comparison')])) stop('Duplicate Ensembl + comparison after version stripping; do not average DESeq2 estimates.')
  grid <- expand.grid(Gene=wanted_genes,comparison=comparisons,stringsAsFactors=FALSE)
  grid$Ensembl <- unname(gene_id_by_name[grid$Gene])
  tab <- dplyr::left_join(grid,tab,by=c('Ensembl','comparison'))
  missing_results <- tab[is.na(tab$gene_id),c('Gene','comparison')]
  if (nrow(missing_results)) warning(nrow(missing_results),' target/contrast entries missing in DESeq2 file; plotted as NA')
  tab$Gene <- factor(tab$Gene,levels=rev(wanted_genes))
  tab$comparison <- factor(tab$comparison,levels=comparisons,labels=comp_labels)
  tab$FDR_stars <- ifelse(is.na(tab$padj),'',ifelse(tab$padj<.001,'***',ifelse(tab$padj<.01,'**',ifelse(tab$padj<.05,'*',''))))
  tab$Cell <- ifelse(is.finite(tab$log2FoldChange),paste0(sprintf('%.2f',tab$log2FoldChange),tab$FDR_stars),'NA')
  maxfc <- suppressWarnings(max(abs(tab$log2FoldChange),na.rm=TRUE))
  if(!is.finite(maxfc) || maxfc<1) maxfc<-1
  p <- ggplot2::ggplot(tab,ggplot2::aes(x=comparison,y=Gene,fill=log2FoldChange))+
    ggplot2::geom_tile(color='grey90',linewidth=.6)+
    ggplot2::geom_text(ggplot2::aes(label=Cell),size=3.1,color='black')+
    ggplot2::scale_fill_gradient2(low='#648FFF',mid='white',high='#FA8072',
      limits=c(-maxfc,maxfc),midpoint=0,na.value='grey90',name='log2FC')+
    ggplot2::labs(x=NULL,y=NULL,caption='* BH-FDR < 0.05; ** < 0.01; *** < 0.001')+
    ggplot2::theme_minimal(base_size=12)+ggplot2::theme(
      panel.grid=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=20,hjust=1),
      axis.text.y=ggplot2::element_text(face='bold'),plot.caption=ggplot2::element_text(hjust=0))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  readr::write_csv(tab,file.path(out,'Figure2D_18genes_age_sex_adjusted_values.csv'))
  ggplot2::ggsave(file.path(out,'Figure2D_18genes_age_sex_adjusted.png'),p,width=8,height=9,dpi=320)
  ggplot2::ggsave(file.path(out,'Figure2D_18genes_age_sex_adjusted.pdf'),p,width=8,height=9)
  message('FIGURE 2D COMPLETE: ',out)
  invisible(tab)
}
run_figure2d()
