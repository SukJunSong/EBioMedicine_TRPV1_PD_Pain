# Input/output paths are configured via environment variables: PD_TRPV1_ROOT, PD_TRPV1_CLINICAL_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)
if (!nzchar(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""))) stop("Set PD_TRPV1_CLINICAL_DIR before running this script (see README.md).", call. = FALSE)

# 03 - PPMI BL/V02 paired Salmon gene-level matrices
# R >= 4.1. Source this entire file. No clinical subjects are silently dropped.
# Clinical age/sex are optional for matrix creation; missing values are audited.

run_ppmi_03 <- function() {
  required_pkgs <- c('readr','readxl','dplyr','tidyr','tibble','AnnotationDbi','org.Hs.eg.db')
  unavailable <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(unavailable)) stop('Install missing packages before running: ', paste(unavailable, collapse=', '), call.=FALSE)

  base_dir <- file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'PPMI')
  output_dir <- file.path(base_dir, 'BL_V02_paired_TRPV1_analysis')
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  input_file <- file.path(base_dir, 'BL_V02_common_PATNO_with_Diagnosis_and_PainGroup.csv')
  if (!file.exists(input_file)) stop('Run 02 first; missing: ', input_file, call.=FALSE)

  clean_id <- function(x) trimws(as.character(x))
  first_existing <- function(paths) {
    hit <- paths[file.exists(paths)]
    if (length(hit)) hit[[1]] else NA_character_
  }
  clean_number <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
  norm_path <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | !nzchar(x)] <- NA_character_
    x
  }
  target_genes <- c('NGF','PLCB4','MAPK10','MAPK12','TRPV1','MAPK13',
                    'IL1R1','MAPK11','IL1RAP','PLCB3','PLCG2','PLCB1',
                    'IL1B','MAPK14','MAPK9','PLCB2','MAPK8','PLCG1')

  meta <- readr::read_csv(input_file, show_col_types=FALSE, col_types=readr::cols(.default=readr::col_character()))
  required <- c('PATNO','DIAGNOSIS','Final_Group','BL_file_path','V02_file_path')
  miss <- setdiff(required, names(meta))
  if (length(miss)) stop('02 output missing columns: ',paste(miss,collapse=', '),call.=FALSE)
  meta$PATNO <- clean_id(meta$PATNO)
  if (anyNA(meta$PATNO) || any(!nzchar(meta$PATNO)) || anyDuplicated(meta$PATNO)) {
    stop('RNA PATNO is missing or duplicated in 02 output.',call.=FALSE)
  }
  meta$BL_file_path <- norm_path(meta$BL_file_path)
  meta$V02_file_path <- norm_path(meta$V02_file_path)
  meta$CLINICAL_PATNO <- meta$PATNO

  # Optional only when patient identity equivalence has been verified manually.
  map_file <- file.path(base_dir,'PATNO_clinical_mapping_verified.csv')
  if (file.exists(map_file)) {
    mp <- readr::read_csv(map_file,show_col_types=FALSE,col_types=readr::cols(.default=readr::col_character()))
    if (!all(c('RNA_PATNO','CLINICAL_PATNO') %in% names(mp))) {
      stop('Verified mapping must have RNA_PATNO, CLINICAL_PATNO.',call.=FALSE)
    }
    mp$RNA_PATNO <- clean_id(mp$RNA_PATNO)
    mp$CLINICAL_PATNO <- clean_id(mp$CLINICAL_PATNO)
    if (anyNA(mp$RNA_PATNO) || anyNA(mp$CLINICAL_PATNO) ||
        any(!nzchar(mp$RNA_PATNO)) || any(!nzchar(mp$CLINICAL_PATNO)) ||
        anyDuplicated(mp$RNA_PATNO) || anyDuplicated(mp$CLINICAL_PATNO) ||
        any(!mp$RNA_PATNO %in% meta$PATNO)) stop('Invalid/nonunique clinical mapping.',call.=FALSE)
    ix <- match(meta$PATNO,mp$RNA_PATNO)
    selected <- !is.na(ix)
    meta$CLINICAL_PATNO[selected] <- mp$CLINICAL_PATNO[ix[selected]]
    if (anyDuplicated(meta$CLINICAL_PATNO)) stop('Mapping creates duplicate clinical IDs.',call.=FALSE)
  }

  # Sex: use V02 first, then BL as fallback. Never guess sex.
  sex_path <- first_existing(c(file.path(base_dir,'metaData V02.xlsx'),
                               file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'metaData V02.xlsx'),
                               file.path(base_dir,'metaData BL.xlsx')))
  sex_tbl <- data.frame(CLINICAL_PATNO=character(),GENDER=character())
  if (!is.na(sex_path)) {
    sx <- as.data.frame(readxl::read_excel(sex_path))
    if (all(c('PATNO','GENDER') %in% names(sx))) {
      sx$CLINICAL_PATNO <- clean_id(sx$PATNO)
      sx$GENDER <- trimws(as.character(sx$GENDER))
      sx$GENDER[is.na(sx$GENDER) | !nzchar(sx$GENDER)] <- NA_character_
      sx <- sx[!is.na(sx$CLINICAL_PATNO) & nzchar(sx$CLINICAL_PATNO),c('CLINICAL_PATNO','GENDER')]
      sex_tbl <- dplyr::summarise(dplyr::group_by(sx,CLINICAL_PATNO),
          GENDER=if (dplyr::n_distinct(GENDER[!is.na(GENDER)])==1L) GENDER[!is.na(GENDER)][1] else NA_character_,
          n_values=dplyr::n_distinct(GENDER[!is.na(GENDER)]),.groups='drop')
      if (any(sex_tbl$n_values>1L)) stop('Conflicting GENDER values for a clinical PATNO; inspect metadata.',call.=FALSE)
      sex_tbl$n_values <- NULL
    } else warning('Sex metadata has no PATNO/GENDER; recording GENDER=NA.')
  } else warning('Sex metadata not found; recording GENDER=NA.')

  # Enrollment age: available when supplied; no age imputation or silently dropping patients.
  age_path <- first_existing(c(file.path(base_dir,'Participant_Status_19May2026.csv'),
                               file.path(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""), 'Participant_Status_19May2026.csv'),
                               file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'Participant_Status_19May2026.csv')))
  age_tbl <- data.frame(CLINICAL_PATNO=character(),ENROLL_AGE=numeric())
  if (!is.na(age_path)) {
    ag <- readr::read_csv(age_path,show_col_types=FALSE,col_types=readr::cols(.default=readr::col_character()))
    if (all(c('PATNO','ENROLL_AGE') %in% names(ag))) {
      ag$CLINICAL_PATNO <- clean_id(ag$PATNO)
      ag$ENROLL_AGE <- clean_number(ag$ENROLL_AGE)
      ag <- ag[!is.na(ag$CLINICAL_PATNO) & nzchar(ag$CLINICAL_PATNO),c('CLINICAL_PATNO','ENROLL_AGE')]
      age_tbl <- dplyr::summarise(dplyr::group_by(ag,CLINICAL_PATNO),
          ENROLL_AGE=if (dplyr::n_distinct(ENROLL_AGE[!is.na(ENROLL_AGE)])==1L) ENROLL_AGE[!is.na(ENROLL_AGE)][1] else NA_real_,
          n_values=dplyr::n_distinct(ENROLL_AGE[!is.na(ENROLL_AGE)]),.groups='drop')
      if (any(age_tbl$n_values>1L)) stop('Conflicting ENROLL_AGE values for a clinical PATNO.',call.=FALSE)
      age_tbl$n_values <- NULL
    } else warning('Participant_Status has no PATNO/ENROLL_AGE; recording age=NA.')
  } else warning('Participant_Status not found; recording age=NA.')

  meta <- dplyr::left_join(meta,sex_tbl,by='CLINICAL_PATNO')
  meta <- dplyr::left_join(meta,age_tbl,by='CLINICAL_PATNO')
  if (anyDuplicated(meta$PATNO)) stop('Clinical join duplicated RNA PATNO.',call.=FALSE)
  meta$GENDER <- as.character(meta$GENDER)
  meta$ENROLL_AGE <- as.numeric(meta$ENROLL_AGE)
  qc <- dplyr::transmute(meta,RNA_PATNO=PATNO,CLINICAL_PATNO,Final_Group,
                         GENDER,ENROLL_AGE,missing_sex=is.na(GENDER),missing_age=is.na(ENROLL_AGE),
                         candidate_without_suffix=sub('[A-Za-z]+$','',PATNO),
                         candidate_in_sex=sub('[A-Za-z]+$','',PATNO) %in% sex_tbl$CLINICAL_PATNO,
                         candidate_in_age=sub('[A-Za-z]+$','',PATNO) %in% age_tbl$CLINICAL_PATNO)
  readr::write_csv(qc,file.path(output_dir,'03_Clinical_Matching_QC.csv'))
  readr::write_csv(meta,file.path(base_dir,'BL_V02_common_PATNO_with_Diagnosis_PainGroup_Age_Gender.csv'))
  cat('Clinical QC: missing sex=',sum(qc$missing_sex),', missing age=',sum(qc$missing_age),' of ',nrow(qc),' subjects\n',sep='')
  if (any(qc$missing_sex | qc$missing_age)) {
    warning('Some age/sex values are missing; expression matrices will still be generated. Do NOT use incomplete clinical data for age/sex-adjusted models.')
  }

  # Paired RNA sample construction: all RNA PATNO preserved.
  sample_info <- dplyr::bind_rows(
    dplyr::transmute(meta,sample_id=paste0(PATNO,'_BL'),PATNO,Visit='BL',
                     Diagnosis_Group=DIAGNOSIS,Final_Group,Pain_Group=Final_Group,
                     ENROLL_AGE,GENDER,file_path=BL_file_path),
    dplyr::transmute(meta,sample_id=paste0(PATNO,'_V02'),PATNO,Visit='V02',
                     Diagnosis_Group=DIAGNOSIS,Final_Group,Pain_Group=Final_Group,
                     ENROLL_AGE,GENDER,file_path=V02_file_path))
  if (anyDuplicated(sample_info$sample_id)) stop('Duplicate sample_id.',call.=FALSE)
  if (anyNA(sample_info$file_path)) stop('BL/V02 file_path missing in 02 output.',call.=FALSE)
  missing_files <- !file.exists(sample_info$file_path)
  if (any(missing_files)) {
    readr::write_csv(sample_info[missing_files,,drop=FALSE],file.path(output_dir,'03_Missing_SF_Files.csv'))
    stop('Salmon .sf input files missing; list written to 03_Missing_SF_Files.csv. No samples removed.',call.=FALSE)
  }
  sample_info$Visit <- factor(sample_info$Visit,levels=c('BL','V02'))
  sample_info$Final_Group <- factor(sample_info$Final_Group,levels=c('Control','No Pain PD','Low Pain PD','High Pain PD'))
  sample_info$PATNO <- factor(sample_info$PATNO)
  sample_info <- dplyr::arrange(sample_info,PATNO,Visit)
  # Leave labels as-is, including Unknown or PD_Pain_NA; do not exclude.
  if (anyNA(sample_info$Final_Group)) warning('Some Final_Group values are outside 4 main groups; samples remain in matrices; see 02 output.')

  # Build per-sample wide matrices without stacking all 40k-gene files into one huge dataframe.
  cat('Reading ',nrow(sample_info),' Salmon files...\n',sep='')
  counts <- vector('list',nrow(sample_info)); tpms <- vector('list',nrow(sample_info))
  for (i in seq_len(nrow(sample_info))) {
    f <- sample_info$file_path[i]
    dat <- readr::read_tsv(f,show_col_types=FALSE,
                           col_types=readr::cols(.default=readr::col_skip(),Name=readr::col_character(),
                                                 TPM=readr::col_double(),NumReads=readr::col_double()))
    needed <- c('Name','TPM','NumReads')
    if (!all(needed %in% names(dat))) stop('Missing Name/TPM/NumReads in ',f,call.=FALSE)
    dat$Ensembl_ID <- sub('\\.[0-9]+$','',dat$Name)
    dat <- dat[!is.na(dat$Ensembl_ID) & nzchar(dat$Ensembl_ID),,drop=FALSE]
    dat <- dplyr::summarise(dplyr::group_by(dat,Ensembl_ID),
                            count=sum(NumReads,na.rm=TRUE),TPM=mean(TPM,na.rm=TRUE),.groups='drop')
    counts[[i]] <- dat[,c('Ensembl_ID','count')]
    tpms[[i]] <- dat[,c('Ensembl_ID','TPM')]
    names(counts[[i]])[2] <- sample_info$sample_id[i]
    names(tpms[[i]])[2] <- sample_info$sample_id[i]
    if (i %% 25L == 0L || i == nrow(sample_info)) cat('  read ',i,'/',nrow(sample_info),'\n',sep='')
  }
  # Union of gene IDs, assembled once (faster and less memory than repeated joins).
  ids <- unique(unlist(lapply(counts, `[[`, 'Ensembl_ID'), use.names=FALSE))
  n_genes <- length(ids)
  n_samples <- nrow(sample_info)
  count_mat <- matrix(0, nrow=n_genes, ncol=n_samples,
                      dimnames=list(ids,sample_info$sample_id))
  tpm_mat <- matrix(0, nrow=n_genes, ncol=n_samples,
                    dimnames=list(ids,sample_info$sample_id))
  for (i in seq_len(n_samples)) {
    k <- match(counts[[i]]$Ensembl_ID,ids)
    count_mat[k,i] <- counts[[i]][[2L]]
    k2 <- match(tpms[[i]]$Ensembl_ID,ids)
    tpm_mat[k2,i] <- tpms[[i]][[2L]]
  }
  count_mat <- round(count_mat)   # Preserve original pipeline's rounding.
  count_wide <- data.frame(Ensembl_ID=ids,count_mat,check.names=FALSE)
  tpm_wide <- data.frame(Ensembl_ID=ids,tpm_mat,check.names=FALSE)

  ids <- unique(as.character(count_wide$Ensembl_ID))
  mapped <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                  keys=ids,keytype='ENSEMBL',column='SYMBOL',multiVals='first')
  symbols <- as.character(mapped[ids]); symbols[is.na(symbols) | !nzchar(symbols)] <- ids[is.na(symbols) | !nzchar(symbols)]
  gene_annotation <- data.frame(Ensembl_ID=ids,GeneSymbol=symbols,stringsAsFactors=FALSE)
  gene_annotation$GeneName <- NA_character_
  targets <- gene_annotation[gene_annotation$GeneSymbol %in% target_genes,,drop=FALSE]
  targets$detected_in_matrix <- targets$Ensembl_ID %in% count_wide$Ensembl_ID
  missing_targets <- setdiff(target_genes,targets$GeneSymbol)
  if (length(missing_targets)) warning('Unmapped target gene symbol(s): ',paste(missing_targets,collapse=', '))

  readr::write_csv(sample_info,file.path(output_dir,'sample_info_BL_V02_paired.csv'))
  readr::write_csv(count_wide,file.path(output_dir,'count_matrix_BL_V02_paired.csv'))
  readr::write_csv(tpm_wide,file.path(output_dir,'TPM_matrix_BL_V02_paired.csv'))
  readr::write_csv(gene_annotation,file.path(output_dir,'gene_annotation_ensembl_to_symbol.csv'))
  readr::write_csv(targets,file.path(output_dir,'target_genes_detection_check.csv'))
  cat('DONE: paired matrices saved to ',output_dir,'\n',sep='')
  cat('Subjects: ',nrow(meta),' | samples: ',nrow(sample_info),' | genes: ',length(ids),'\n',sep='')
  invisible(list(sample_info=sample_info,clinical_qc=qc))
}

run_ppmi_03()
