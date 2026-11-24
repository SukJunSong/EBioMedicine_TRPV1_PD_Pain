# GSE8397 GPL96 downloaded GEO Series -> PD vs Control -> custom 18-gene fgsea
# FIXED: Control sample PDC1 was misclassified as PD by previous broad PD regex.
# Outputs intentionally cannot be compared to old incorrect results; rerun from scratch.
# Use this after saving gse to GSE8397/GSE8397_downloaded.rds
pkgs <- c('Biobase','GEOquery','limma','fgsea','dplyr','ggplot2','tidyr')
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install missing packages first: ', paste(missing,collapse=', '))

input <- 'GSE8397/GSE8397_downloaded.rds'
if(!file.exists(input)) stop('Cannot find ', input, '. Run getwd() and check where the RDS was saved.')
outdir <- 'GSE8397/hsa04750_custom_results_CORRECTED'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

G <- readRDS(input)
if(inherits(G,'ExpressionSet')) G <- list(G)
if(!is.list(G) || !length(G)) stop('Unexpected GEO RDS format')
platforms <- vapply(G,function(z) Biobase::annotation(z),character(1))
print(platforms)
idx <- which(platforms=='GPL96')
if(length(idx)!=1L) stop('Expected one GPL96 ExpressionSet; inspect platform IDs')
E <- G[[idx]]
expr <- Biobase::exprs(E)
ph <- Biobase::pData(E)
if(!'title' %in% colnames(ph)) stop('Missing title in phenotype metadata')
ph$region <- dplyr::case_when(
 grepl('lateral substantia nigra',ph$title,ignore.case=TRUE) ~ 'Lateral substantia nigra',
 grepl('medial substantia nigra',ph$title,ignore.case=TRUE) ~ 'Medial substantia nigra',
 grepl('superior frontal gyrus',ph$title,ignore.case=TRUE) ~ 'Superior frontal gyrus',
 TRUE ~ NA_character_)
# IMPORTANT: GEO controls include the label "control case PDC1".
# Searching for "PD" before "control" MISCLASSIFIES that control as PD.
# Use unambiguous title patterns, inspect any unmatched or ambiguous rows.
control_hit <- grepl("\\bcontrol\\b", ph$title, ignore.case = TRUE, perl = TRUE)
pd_hit <- grepl("Parkinson.s disease|parkinsonian", ph$title, ignore.case = TRUE, perl = TRUE)
if (any(control_hit & pd_hit)) {
  print(data.frame(accession = rownames(ph)[control_hit & pd_hit],
                   title = ph$title[control_hit & pd_hit]))
  stop("Ambiguous disease labels: both control and PD detected. Review titles.")
}
ph$disease <- ifelse(control_hit, "Control", ifelse(pd_hit, "PD", NA_character_))

# Use original GEO series design sample counts as a critical data-integrity check.
expected_regions <- c("Lateral substantia nigra",
                      "Medial substantia nigra",
                      "Superior frontal gyrus")
count_check <- table(factor(ph$region, levels=expected_regions),
                     factor(ph$disease, levels=c("Control","PD")))
print(count_check)
expected_counts <- matrix(c(7,8,3,9,15,5), nrow=3,
                          dimnames=list(expected_regions,c("Control","PD")))
if (!all(as.matrix(count_check) == expected_counts)) {
  print(data.frame(accession=rownames(ph), title=ph$title,
                   region=ph$region, disease=ph$disease),row.names=FALSE)
  stop("GEO design sample-count mismatch. Expected L-SN 7/9, M-SN 8/15, SFG 3/5 (Control/PD).")
}
if (any(is.na(ph$region) | is.na(ph$disease))) {
  print(ph[is.na(ph$region) | is.na(ph$disease), "title", drop=FALSE])
  stop("Unclassified samples. Check GEO title metadata.")
}
utils::write.csv(data.frame(accession=rownames(ph), title=ph$title,
                            region=ph$region, disease=ph$disease),
                 file.path(outdir,"GSE8397_sample_classification_audit.csv"),row.names=FALSE)

if(!all(expected_regions %in% ph$region)) stop('One or more expected brain regions not found')

# Use original GPL96 annotations for probe -> gene symbol mapping.
gpl <- GEOquery::getGEO('GPL96',destdir='GSE8397')
a <- GEOquery::Table(gpl)
sym_choices <- c('Gene Symbol','Gene.symbol','Gene.Symbol','GENE_SYMBOL','Symbol','gene_symbol')
sym <- intersect(sym_choices,colnames(a))
if(!length(sym)) stop('Gene-symbol annotation column not found. Columns: ',paste(colnames(a),collapse=', '))
anno <- data.frame(probe=as.character(a$ID),symbol=as.character(a[[sym[1]]]))
anno <- anno[!is.na(anno$probe) & !is.na(anno$symbol) & nzchar(anno$symbol),]

custom <- c('MAPK9','PLCB4','MAPK8','MAPK10','MAPK14','IL1R1','PLCG1','IL1B','NGF',
            'PLCB2','PLCB3','MAPK13','PLCB1','MAPK11','MAPK12','IL1RAP','TRPV1','PLCG2')
pathways <- list(TRPV1_related_hsa04750_18genes=custom)
results <- list()
rank_exports <- list()
for(reg in expected_regions){
 i <- which(ph$region==reg & ph$disease %in% c('PD','Control'))
 g <- factor(ph$disease[i], levels=c('Control','PD'))
 cat('\n',reg,': '); print(table(g))
 if(any(table(g)<2)) stop('Insufficient samples in ',reg)
 design <- stats::model.matrix(~0+g)
 colnames(design) <- c('Control','PD')
 fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr[,i,drop=FALSE],design),
                                            limma::makeContrasts(PD-Control,levels=design)))
 tt <- limma::topTable(fit,coef=1,number=Inf,sort.by='none')
 tt$probe <- rownames(tt)
 m <- dplyr::left_join(tt,anno,by='probe')
 # Split ambiguous HG-U133A annotations (e.g., SHPK /// TRPV1), consistently with legacy workflow.
 m$symbol <- trimws(m$symbol)
 m <- m[!is.na(m$symbol) & nzchar(m$symbol) & is.finite(m$t),]
 m <- tidyr::separate_rows(m,symbol,sep='\\s*///\\s*')
 m <- m[nzchar(m$symbol),]
 # One probe per gene: strongest absolute moderated t, document this decision.
 m <- m[order(-abs(m$t),m$probe),]
 m <- m[!duplicated(m$symbol),]
 ranks <- stats::setNames(m$t,m$symbol)
 ranks <- sort(ranks,decreasing=TRUE)
 hit <- intersect(custom,names(ranks))
 cat('Mapped custom genes (',length(hit),'/18): ',paste(hit,collapse=', '),'\n',sep='')
 if(length(hit)<5) stop('Less than 5 mapped custom genes in ',reg)
 set.seed(8397)
 fg <- fgsea::fgseaMultilevel(pathways=pathways,stats=ranks,minSize=5,maxSize=500,eps=0)
 if(nrow(fg)!=1) stop('Unexpected fgsea result count in ',reg)
 results[[reg]] <- data.frame(region=reg,PD_n=sum(g=='PD'),Control_n=sum(g=='Control'),
 NES=fg$NES,pval=fg$pval,fgsea_single_set_padj=fg$padj,ES=fg$ES,
 mapped_genes=fg$size,leadingEdge=paste(fg$leadingEdge[[1]],collapse=','),
 measured_hsa04750_genes=paste(hit,collapse=','))
 rank_exports[[reg]] <- data.frame(region=reg,gene=m$symbol,probe=m$probe,t=m$t,logFC=m$logFC,P.Value=m$P.Value)
}
res <- dplyr::bind_rows(results)
res$FDR_BH_three_regions <- p.adjust(res$pval,method='BH')
print(res[,c('region','PD_n','Control_n','mapped_genes','NES','pval','FDR_BH_three_regions')],row.names=FALSE)
utils::write.csv(res,file.path(outdir,'GSE8397_hsa04750_fgsea_NES_FDR_3regions.csv'),row.names=FALSE)
utils::write.csv(dplyr::bind_rows(rank_exports),file.path(outdir,'GSE8397_all_genes_ranked_t.csv'),row.names=FALSE)

res$region <- factor(res$region,levels=expected_regions)
res$label <- sprintf('NES %.2f  |  FDR %.3g',res$NES,res$FDR_BH_three_regions)
p <- ggplot2::ggplot(res,ggplot2::aes(x=region,y=NES,fill=NES>0))+
 ggplot2::geom_hline(yintercept=0,color='grey55')+
 ggplot2::geom_col(width=0.6)+
 ggplot2::geom_text(ggplot2::aes(label=label),vjust=ifelse(res$NES>=0,-0.5,1.5),fontface='bold',size=4)+
 ggplot2::scale_fill_manual(values=c('TRUE'='#FA8072','FALSE'='#648FFF'),guide='none')+
 ggplot2::coord_cartesian(ylim=range(c(0,res$NES))*1.4,clip='off')+
 ggplot2::labs(title='TRPV1-related hsa04750 (18 genes)',x=NULL,y='Normalized enrichment score (NES)')+
 ggplot2::theme_classic(base_size=13)+ggplot2::theme(axis.text.x=ggplot2::element_text(face='bold',angle=15,hjust=1),plot.title=ggplot2::element_text(face='bold'))
ggplot2::ggsave(file.path(outdir,'GSE8397_hsa04750_fgsea_3regions.pdf'),p,width=10,height=5.5)
ggplot2::ggsave(file.path(outdir,'GSE8397_hsa04750_fgsea_3regions.svg'),p,width=10,height=5.5)
cat('\nDONE: ',normalizePath(outdir),'\n')
