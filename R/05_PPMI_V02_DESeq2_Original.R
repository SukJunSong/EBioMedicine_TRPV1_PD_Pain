# Input/output paths are configured via environment variables: PD_TRPV1_ROOT, PD_TRPV1_CLINICAL_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)
if (!nzchar(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""))) stop("Set PD_TRPV1_CLINICAL_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 05_PPMI_V02_DESeq2_Original.R
# Requires original project datasets; full R execution not validated here.
# V02 DESeq2: model formula type FIX, 2026-10-27
# Run in a freshly restarted R session, using source() (do not execute selected lines).
# This revision runs all stages in one function so an upstream failure cannot
# leave missing objects and continue to downstream analyses.

PPMI_run_V02 <- function() {
  # PPMI V02 | (1) original DESeq2 reproduction -> (2) medication sensitivity
  # Input: kallisto/Salmon .sf files at <PD_TRPV1_ROOT>/PPMI/genes V02
  # NOT a TPM regression. Group is NP1PAIN AT V02 (HC separate).
  # IMPORTANT: the original code replaced DESeq2 padj with manual BH; both are saved.
  # R packages (install beforehand): DESeq2, readxl, readr, dplyr, ggplot2, ggrepel
  # Run entire script after setting PD_TRPV1_ROOT and PD_TRPV1_CLINICAL_DIR.
  # If currently at Browse[n]>, type Q then restart R; do NOT source from Browse[].
  suppressPackageStartupMessages({
    library(DESeq2); library(readxl); library(readr); library(dplyr)
    library(ggplot2); library(ggrepel)
  })
  root <- Sys.getenv("PD_TRPV1_ROOT", unset="")
  ppmi <- file.path(root, 'PPMI')
  choose_file <- function(paths, label) {
    hit <- paths[file.exists(paths)]
    if (!length(hit)) stop('Missing ',label,'; checked:\n',paste(paths,collapse='\n'))
    hit[[1]]
  }
  meta_path <- choose_file(c(file.path(ppmi,'metaData V02.xlsx'),file.path(root,'metaData V02.xlsx')), 'metaData V02.xlsx')
  pain_path <- choose_file(c(file.path(ppmi,'painData V02.xlsx'),file.path(root,'painData V02.xlsx')), 'painData V02.xlsx')
  sf_dir <- choose_file(c(file.path(ppmi,'genes V02'),file.path(root,'genes V02')), 'genes V02 folder')
  conmed_path <- choose_file(c(file.path(root,'Concomitant_Medication_Log_26Oct2026.csv'),file.path(ppmi,'Concomitant_Medication_Log_26Oct2026.csv')), 'CONMED CSV')
  ledd_path <- choose_file(c(file.path(root,'LEDD_Concomitant_Medication_Log_26Oct2026.csv'),file.path(ppmi,'LEDD_Concomitant_Medication_Log_26Oct2026.csv')), 'LEDD CSV')
  clinical_path <- choose_file(c(file.path(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""), 'MDS-UPDRS_Part_I_Patient_Questionnaire_25Feb2025.csv'),file.path(root,'MDS-UPDRS_Part_I_Patient_Questionnaire_25Feb2025.csv')), 'MDS-UPDRS Part I patient V02 visit dates')
  out <- file.path(root,'V02_DESeq2_Original_vs_Medication')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  # Earlier runs may have created partial CSVs; do not interpret them as fresh results.
  stale <- c('03_ORIGINAL_hsa04750_18_genes.csv', '09_TRPV1_original_vs_medication_summary.csv')
  # Previous output files are preserved. Check timestamps and 00_FATAL_ERROR.txt after run.
  # for (f in stale) { p <- file.path(out,f); if(file.exists(p)) file.remove(p) }
  save_csv <- function(x, filename) readr::write_csv(as.data.frame(x),file.path(out,filename),na='')
  key <- function(x) trimws(as.character(x))
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  # --- metadata: fail rather than make up missing visit classification ---
  meta <- as.data.frame(readxl::read_excel(meta_path))
  pain <- as.data.frame(readxl::read_excel(pain_path))
  stopifnot(all(c('PATNO','DIAGNOSIS') %in% names(meta)), 'PATNO' %in% names(pain))
  meta$PATNO <- key(meta$PATNO); pain$PATNO <- key(pain$PATNO)
  if(anyDuplicated(meta$PATNO)) stop('Duplicate PATNO in V02 metadata; fix before analysis')
  if(anyDuplicated(pain$PATNO)) stop('Duplicate PATNO in V02 painData; fix before analysis')
  # NP1PAIN may be in pain file or metadata. Explicitly prefer pain V02 values.
  if ('NP1PAIN' %in% names(pain)) {
    pain_keep <- pain[,c('PATNO','NP1PAIN')]; names(pain_keep)[2] <- 'NP1PAIN_V02'
  } else stop('painData V02.xlsx needs NP1PAIN column')
  meta2 <- dplyr::left_join(meta, pain_keep, by='PATNO')
  meta2$DIAGNOSIS <- toupper(trimws(as.character(meta2$DIAGNOSIS)))
  meta2$Group <- ifelse(meta2$DIAGNOSIS %in% c('CONTROL','HC','HEALTHY CONTROL'), 'Control',
    ifelse(meta2$DIAGNOSIS=='PD' & !is.na(num(meta2$NP1PAIN_V02)),
      ifelse(num(meta2$NP1PAIN_V02)==0,'No pain PD',ifelse(num(meta2$NP1PAIN_V02)==1,'Low pain PD',
        ifelse(num(meta2$NP1PAIN_V02)>=2,'High pain PD', NA_character_))), NA_character_))
  # Original code included ALL HC regardless of missing HC pain; PD must have V02 pain.
  # --- SF: count matrix, same Name IDs across samples; explicitly check duplicates ---
  sf_paths <- list.files(sf_dir,pattern='\\.sf$',full.names=TRUE,ignore.case=TRUE)
  if(!length(sf_paths)) stop('No .sf files in: ',sf_dir)
  ids <- sub('^PPMI-Phase[12]-IR3\\.([0-9]+)\\..*$', '\\1', basename(sf_paths))
  bad <- !grepl('^PPMI-Phase[12]-IR3\\.[0-9]+\\.',basename(sf_paths))
  if(any(bad)) stop('Unable to extract PATNO from filenames: ',paste(head(basename(sf_paths[bad]),10),collapse=', '))
  manifest <- data.frame(PATNO=ids,File=sf_paths)
  save_csv(manifest,'00_all_V02_sf_manifest.csv')
  manifest <- manifest[manifest$PATNO %in% meta2$PATNO & !is.na(meta2$Group[match(manifest$PATNO,meta2$PATNO)]),]
  if(anyDuplicated(manifest$PATNO)) stop('Multiple V02 .sf files for PATNO; inspect 00_all_V02_sf_manifest.csv before selecting a sample')
  if(nrow(manifest)<15) stop('Too few matched samples; check metadata and .sf directory')
  message('Reading ',nrow(manifest),' sf files (can take a few minutes)...')
  read_sf <- function(path) {
    z <- readr::read_tsv(path, col_types=readr::cols(.default='c'),show_col_types=FALSE,progress=FALSE)
    if(!all(c('Name','NumReads') %in% names(z))) stop('Bad SF file: ',path)
    if(anyDuplicated(z$Name)) stop('Duplicate gene Name in: ',path)
    val <- num(z$NumReads)
    if(any(!is.finite(val))) stop('Nonfinite read counts: ',path)
    data.frame(Name=as.character(z$Name), value=round(val))
  }
  first <- read_sf(manifest$File[[1]]); gene_ids <- first$Name
  counts <- matrix(NA_real_,nrow=length(gene_ids),ncol=nrow(manifest),dimnames=list(gene_ids,manifest$PATNO))
  counts[,1] <- first$value
  if(nrow(manifest)>1) for(i in 2:nrow(manifest)) {
    z <- read_sf(manifest$File[[i]])
    idx <- match(gene_ids,z$Name)
    if(anyNA(idx) || nrow(z)!=length(gene_ids)) stop('Gene index mismatch in ',manifest$File[[i]],'; cannot silently substitute 0')
    counts[,i] <- z$value[idx]
  }
  if(any(counts<0) || any(!is.finite(counts))) stop('Invalid count matrix')
  # DESeq2 count input must be integers (NumReads from quant.sf is rounded, as in original).
  storage.mode(counts) <- 'integer'
  md <- meta2[match(colnames(counts),meta2$PATNO),,drop=FALSE]
  stopifnot(identical(md$PATNO,colnames(counts)))
  md$group <- factor(md$Group,levels=c('Control','No pain PD','Low pain PD','High pain PD'))
  rownames(md) <- md$PATNO
  save_csv(md %>% dplyr::select(PATNO,DIAGNOSIS,NP1PAIN_V02,Group),'01_V02_original_sample_cohort.csv')
  save_csv(as.data.frame(table(Group=md$group)),'02_original_group_N.csv')
  message('Original group N: '); print(table(md$group))
  if(any(table(md$group)<2)) stop('At least one group too small/missing')
  # -- fit one reusable function --
  fit_deseq <- function(counts_sub,coldata,design_formula,tag) {
    # FIX: accept formula text (e.g. '~ group') AND actual formula objects.
    # Passing text to model.matrix() causes x$terms / atomic vectors error.
    if (is.character(design_formula) && length(design_formula)==1L) {
      design_formula <- stats::as.formula(design_formula, env=baseenv())
    }
    if (!inherits(design_formula, 'formula')) stop(tag, ': design_formula must be a formula')
    if(!identical(colnames(counts_sub),rownames(coldata))) stop(tag, ': PATNO and colData order differ')
    if(ncol(counts_sub)<12) stop(tag, ': too few samples')
    # Supply only model variables; Excel metadata can contain duplicate/unsupported columns.
    needed <- all.vars(design_formula)
    if(!all(needed %in% names(coldata))) stop(tag, ': missing design columns: ', paste(setdiff(needed,names(coldata)),collapse=', '))
    minimal <- as.data.frame(coldata[,needed,drop=FALSE])
    rownames(minimal) <- rownames(coldata)
    if(anyNA(minimal)) stop(tag, ': missing covariate values')
    if('group' %in% names(minimal)) minimal$group <- factor(as.character(minimal$group), levels=c('Control','No pain PD','Low pain PD','High pain PD'))
    X <- model.matrix(design_formula,data=minimal)
    if(qr(X)$rank<ncol(X)) stop(tag, ': design matrix rank deficient')
    message('START ',tag, ' | n samples=',ncol(counts_sub),' | genes=',nrow(counts_sub))
    tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts_sub,colData=minimal,design=design_formula)
      # Save a diagnostic checkpoint before fitting.
      saveRDS(list(n=ncol(dds),levels=levels(minimal$group),design=as.character(design_formula),
                   group_n=table(minimal$group)),file.path(out,paste0(tag,'_preflight.rds')))
      dds <- DESeq2::DESeq(dds,quiet=TRUE)
      saveRDS(dds,file.path(out,paste0(tag,'_dds.rds')))
      message('SUCCESS ',tag)
      dds
    }, error=function(e){
      msg <- paste0('FAILED ',tag,': ',conditionMessage(e))
      writeLines(msg,file.path(out,paste0(tag,'_ERROR.txt')))
      message(msg)
      stop(msg,call.=FALSE)
    })
  }

  extract_de <- function(dds,cmp_tag,group_a,group_b='Control') {
    r <- DESeq2::results(dds,contrast=c('group',group_a,group_b),alpha=0.05)
    z <- as.data.frame(r)
    z$Ensembl_versioned <- rownames(z)
    z$Ensembl <- sub('\\.[0-9]+$','',z$Ensembl_versioned)
    z$padj_manual_BH_original <- p.adjust(z$pvalue,method='BH')
    z$comparison <- cmp_tag
    z <- z[,c('comparison','Ensembl','Ensembl_versioned','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj','padj_manual_BH_original')]
    save_csv(z,paste0(cmp_tag,'_all_genes.csv'))
    z
  }
  original <- fit_deseq(counts,md,'~ group','A_original')
  comparisons <- c('No pain PD','Low pain PD','High pain PD')
  original_res <- tryCatch(dplyr::bind_rows(lapply(comparisons, function(g) extract_de(original,
     paste0('A_',gsub(' ','_',g),'_vs_HC'),g))),error=function(e){
       writeLines(conditionMessage(e),file.path(out,'A_original_results_ERROR.txt'))
       stop('Original DESeq2 result extraction failed: ',conditionMessage(e),call.=FALSE)
     })
  saveRDS(original_res,file.path(out,'A_original_results.rds'))
  target <- c('ENSG00000050748','ENSG00000101333','ENSG00000107643','ENSG00000109339',
   'ENSG00000112062','ENSG00000115594','ENSG00000124181','ENSG00000125538',
   'ENSG00000134259','ENSG00000137841','ENSG00000149782','ENSG00000156711',
   'ENSG00000182621','ENSG00000185386','ENSG00000188130','ENSG00000196083',
   'ENSG00000196689','ENSG00000197943')
  # TRPV1 is ENSG00000196689, version handled only for matching.
  save_csv(original_res %>% dplyr::filter(Ensembl %in% target),'03_ORIGINAL_hsa04750_18_genes.csv')
  volcano <- function(z,cmp) {
    z <- z[is.finite(z$log2FoldChange)&is.finite(z$padj),]
    z$regulation <- ifelse(z$padj<0.05 & z$log2FoldChange>0,'Up',ifelse(z$padj<0.05 & z$log2FoldChange<0,'Down','NS'))
    z$label <- ifelse(z$Ensembl=='ENSG00000196689','TRPV1',ifelse(z$Ensembl %in% target,z$Ensembl,''))
    p <- ggplot(z,aes(log2FoldChange,-log10(pmax(padj,1e-300)),color=regulation))+
      geom_point(alpha=0.45,size=1)+
      geom_point(data=z[z$Ensembl %in% target,],shape=21,fill='white',color='black',size=1.9,stroke=0.5)+
      ggrepel::geom_text_repel(data=z[z$label=='TRPV1',],aes(label=label),color='black',size=4,max.overlaps=Inf)+
      scale_color_manual(values=c(Up='#FA8072',Down='#648FFF',NS='grey75'))+
      geom_hline(yintercept=-log10(0.05),linetype=2,color='grey35')+
      labs(title=paste('V02',cmp),x='log2 fold change',y='-log10(DESeq2 adjusted p-value)')+
      theme_classic(base_size=12)+theme(plot.title=element_text(hjust=0.5))
    ggsave(file.path(out,paste0(cmp,'_volcano.png')),p,width=12,height=7,dpi=350)
    ggsave(file.path(out,paste0(cmp,'_volcano.pdf')),p,width=12,height=7)
  }
  for(g in comparisons) volcano(original_res[original_res$comparison==paste0('A_',gsub(' ','_',g),'_vs_HC'),],paste0('A_',gsub(' ','_',g),'_vs_HC'))
  # ---- V02 medication dates: same month is POSSIBLE exposure, not verified pre-visit. ----
  vis <- readr::read_csv(clinical_path,col_types=readr::cols(.default='c'),show_col_types=FALSE)
  if(!all(c('PATNO','EVENT_ID','INFODT') %in% names(vis))) stop('V02 visit file needs PATNO EVENT_ID INFODT')
  month <- function(x) {
    x <- trimws(as.character(x)); x[x %in% c('','NA','N/A','UNK','Unknown')] <- NA_character_
    s <- ifelse(grepl('^[0-9]{1,2}/[0-9]{4}$',x),paste0('01/',x),x)
    dd <- suppressWarnings(as.Date(s,'%d/%m/%Y'))
    # ISO dates from export if present
    iso <- suppressWarnings(as.Date(s,'%Y-%m-%d'))
    dd[is.na(dd)] <- iso[is.na(dd)]
    as.Date(ifelse(is.na(dd),NA_character_,format(dd,'%Y-%m-01')))
  }
  vis <- vis[vis$EVENT_ID=='V02' & key(vis$PATNO) %in% md$PATNO,,drop=FALSE]
  vis$PATNO <- key(vis$PATNO); vis$visit_month <- month(vis$INFODT)
  if(any(!is.na(vis$INFODT)&is.na(vis$visit_month))) warning('Some V02 INFODT values not parsed; inspect QC')
  visit_dates <- vis %>% group_by(PATNO) %>% summarise(V02_month=if(length(unique(na.omit(visit_month)))==1) unique(na.omit(visit_month)) else as.Date(NA),n_visit_month=n_distinct(visit_month,na.rm=TRUE),.groups='drop')
  # Months with conflicting visit records are marked missing; never choose one silently.
  med_cm <- read_csv(conmed_path,col_types=cols(.default='c'),show_col_types=FALSE)
  med_ld <- read_csv(ledd_path,col_types=cols(.default='c'),show_col_types=FALSE)
  cm <- med_cm %>% transmute(PATNO=key(PATNO),drug=CMTRT,indication=CMINDC_TEXT,source='CONMED',start=month(STARTDT),stop=month(STOPDT),REC_ID)
  ld <- med_ld %>% transmute(PATNO=key(PATNO),drug=LEDTRT,indication=NA_character_,source='LEDD',start=month(STARTDT),stop=month(STOPDT),REC_ID)
  rx <- paste(c('ACETAMINOPHEN','PARACETAMOL','TYLENOL','IBUPROFEN','NAPROXEN','DICLOFENAC','CELECOXIB','MELOXICAM','ASPIRIN','TRAMADOL','CODEINE','HYDROCODONE','OXYCODONE','MORPHINE','FENTANYL','BUPRENORPHINE','TAPENTADOL','METHADONE','GABAPENTIN','PREGABALIN','DULOXETINE','AMITRIPTYLINE','NORTRIPTYLINE','CARBAMAZEPINE','LIDOCAINE','CAPSAICIN','KETOROLAC','SUMATRIPTAN'),collapse='|')
  painrx <- 'PAIN|HEADACHE|MIGRAINE|NEURALGIA|ARTHRALGIA|SCIATICA|BACKACHE|BACK ACHE'
  pdrx <- 'LEVODOPA|CARBIDOPA|SINEMET|RYTARY|STALEVO|ENTACAPONE|OPICAPONE|RASAGILINE|SELEGILINE|SAFINAMIDE|PRAMIPEXOLE|ROPINIROLE|ROTIGOTINE|APOMORPHINE|AMANTADINE|BENZTROPINE'
  med <- bind_rows(cm,ld) %>% filter(PATNO %in% md$PATNO) %>% left_join(visit_dates,by='PATNO') %>% mutate(
    drug_up=toupper(coalesce(drug,'')), indic_up=toupper(coalesce(indication,'')),
    analgesic_candidate=source=='CONMED' & (grepl(rx,drug_up)|grepl(painrx,indic_up)),
    PD_drug_candidate=source=='LEDD'|grepl(pdrx,drug_up),
    active_inclusive=!is.na(start)&!is.na(V02_month)&start<=V02_month&(is.na(stop)|stop>=V02_month),
    active_strict=!is.na(start)&!is.na(V02_month)&start<V02_month&(is.na(stop)|stop>=V02_month),
    same_month_start=!is.na(start)&!is.na(V02_month)&start==V02_month)
  save_csv(med %>% filter(analgesic_candidate|PD_drug_candidate),'04_V02_medication_drug_audit.csv')
  flags <- med %>% group_by(PATNO) %>% summarise(
    analgesic_possible=as.integer(any(analgesic_candidate & active_inclusive)),
    analgesic_strict=as.integer(any(analgesic_candidate & active_strict)),
    PD_drug_possible=as.integer(any(PD_drug_candidate & active_inclusive)),
    PD_drug_strict=as.integer(any(PD_drug_candidate & active_strict)),
    .groups='drop')
  visit_info <- data.frame(PATNO=md$PATNO,Group=md$Group) %>% left_join(visit_dates,by='PATNO') %>% left_join(flags,by='PATNO')
  # No medication rows means no *recorded* use, only when a visit date exists.
  for(k in c('analgesic_possible','analgesic_strict','PD_drug_possible','PD_drug_strict')) {
    visit_info[[k]][is.na(visit_info[[k]]) & !is.na(visit_info$V02_month)] <- 0L
  }
  save_csv(visit_info,'05_V02_medication_patient_flags_ALL_groups.csv')
  save_csv(visit_info %>% group_by(Group) %>% summarise(n=n(),n_V02_month=sum(!is.na(V02_month)),n_analgesic_possible=sum(analgesic_possible==1,na.rm=TRUE),n_PD_drug_possible=sum(PD_drug_possible==1,na.rm=TRUE),.groups='drop'),'06_medication_QC_group_N.csv')
  md2 <- as.data.frame(md)
  md2 <- dplyr::left_join(md2,visit_info %>% dplyr::select(PATNO,V02_month,analgesic_possible,analgesic_strict,PD_drug_possible,PD_drug_strict),by='PATNO')
  rownames(md2) <- md2$PATNO
  # Full cohort medication-adjusted model (HC included) only for rows with known V02 visit month.
  # Missing HC exposure cannot be assumed unexposed.
  medres <- list()
  for(mode in c('possible','strict')) {
    cols <- c(paste0('analgesic_',mode),paste0('PD_drug_',mode))
    tmp <- md2[!is.na(md2$V02_month) & complete.cases(md2[,cols,drop=FALSE]),,drop=FALSE]
    for(cn in cols) tmp[[cn]] <- as.integer(tmp[[cn]])
    if(nrow(tmp)<15 || any(table(tmp$group)<2)) {
      message('Skipping medication model (',mode,'): no adequate V02 visit-date coverage across all groups. See QC.'); next
    }
    # Drop constant covariates, but never the group variable.
    active <- cols[vapply(tmp[,cols,drop=FALSE],function(x) length(unique(x))>1,logical(1))]
    f <- as.formula(paste('~',paste(c(active,'group'),collapse=' + ')))
    mm <- model.matrix(f,data=tmp)
    if(qr(mm)$rank<ncol(mm)) { message('Skipping rank-deficient ',mode,' model; see 05 QC'); next }
    tag <- paste0('B_medication_',mode)
    fit <- fit_deseq(counts[,tmp$PATNO,drop=FALSE],tmp,f,tag)
    for(g in comparisons) medres[[length(medres)+1]] <- extract_de(fit,paste0(tag,'_',gsub(' ','_',g),'_vs_HC'),g)
    save_csv(tmp %>% dplyr::select(PATNO,Group,all_of(cols)),paste0(tag,'_N_and_exposures.csv'))
  }
  if(length(medres)) {
    z <- bind_rows(medres)
    save_csv(z %>% filter(Ensembl %in% target),'07_MEDICATION_adjusted_hsa04750_18_genes.csv')
  } else {
    writeLines('Medication-adjusted HC vs PD model could not be estimated: inspect V02 date/medication data for BOTH HC and PD. No missing HC values were set to zero.',file.path(out,'07_MEDICATION_model_not_run_README.txt'))
  }
  # Stronger sensitivity: among V02 date-verified patients, exclude recorded analgesic use;
  # do not claim that absence of a record proves no exposure.
  clean <- md2[!is.na(md2$V02_month)&!is.na(md2$analgesic_possible)&md2$analgesic_possible==0,,drop=FALSE]
  if(nrow(clean)>=15 && all(table(clean$group)>=2)) {
    fit_clean <- fit_deseq(counts[,clean$PATNO,drop=FALSE],clean,'~ group','C_exclude_analgesic_candidates')
    cleanres <- bind_rows(lapply(comparisons,function(g) extract_de(fit_clean,paste0('C_analgesic_excluded_',gsub(' ','_',g),'_vs_HC'),g)))
    save_csv(cleanres %>% filter(Ensembl %in% target),'08_EXCLUDED_analgesic_hsa04750_18_genes.csv')
    save_csv(as.data.frame(table(Group=clean$group)),'08_EXCLUDED_analgesic_group_N.csv')
  } else message('Skipping analgesic-excluded model because groups insufficient; no inference from small High Pain group')
  # TRPV1 summary. These estimates can only be compared after checking N and clinical visit dates.
  collected <- list(original_res)
  if(length(medres)) collected <- c(collected,list(bind_rows(medres)))
  if(exists('cleanres',inherits=FALSE)) collected <- c(collected,list(cleanres))
  allr <- bind_rows(collected)
  save_csv(allr %>% filter(Ensembl=='ENSG00000196689') %>%
    dplyr::select(comparison,baseMean,log2FoldChange,lfcSE,pvalue,padj,padj_manual_BH_original),
    '09_TRPV1_original_vs_medication_summary.csv')
  writeLines(c('Original: V02 raw NumReads DESeq2, HC vs V02 pain groups, design=~group.',
   'Original code had manual BH on pvalue, which is included separately from DESeq2 padj.',
   'Medication candidate lists are screenings and need clinical adjudication.',
   'Same-month starts count as POSSIBLE exposure; strict uses dates before visit month.',
   'Missing HC visit date/exposure is NOT coded as no use.',
   'A medication-adjusted fit can differ in sample size; compare n before attributing results to drugs.',
   'Raw counts rounded from Kallisto NumReads as in original code; preferred tximport workflow is a future methodological improvement.'),
   file.path(out,'README_interpretation.txt'))
  message('FINISHED. Output at: ',out)

  invisible(TRUE)
}

PPMI_V02_RUN_STATUS <- tryCatch(
  PPMI_run_V02(),
  error = function(e) {
    diagnostic_dir <- file.path(root, 'V02_DESeq2_Original_vs_Medication')
    dir.create(diagnostic_dir, recursive=TRUE, showWarnings=FALSE)
    err <- c(paste0("ERROR: ", conditionMessage(e)),
             paste0("CALL: ", paste(deparse(conditionCall(e)), collapse=" ")),
             "The analysis STOPPED. Later output files must not be interpreted as complete.",
             "Restart R and rerun the entire source() to reproduce the failure.")
    writeLines(err, file.path(diagnostic_dir, "00_FATAL_ERROR.txt"))
    message(paste(err, collapse="\n"))
    FALSE
  }
)
if (!isTRUE(PPMI_V02_RUN_STATUS)) message("FAILED. Inspect 00_FATAL_ERROR.txt (not FINISHED).")
