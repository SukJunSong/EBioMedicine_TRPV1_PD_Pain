# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 01_PPMI_Match_BL_V02_Samples_FIXED(2).R
# Requires original project datasets; full R execution not validated here.
# ============================================================
# 01. PPMI BL-V02 common participant identification
# Self-contained. Run WHOLE file, not selected lines.
# ============================================================

run_match_bl_v02 <- function() {
  base_dir <- file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'PPMI')
  bl_dir  <- file.path(base_dir, 'genes BL')
  v02_dir <- file.path(base_dir, 'genes V02')
  out_bl_dir  <- file.path(base_dir, 'genes BL_matched_PATNO')
  out_v02_dir <- file.path(base_dir, 'genes V02_matched_PATNO')
  output_csv  <- file.path(base_dir, 'BL_V02_common_PATNO_file_list.csv')

  # Keep original copy behavior; set FALSE if only the CSV is needed.
  copy_matched_sf <- FALSE

  for (p in c(base_dir, bl_dir, v02_dir)) {
    if (!dir.exists(p)) stop('Required folder not found: ', p, call. = FALSE)
  }

  list_sf <- function(folder, visit) {
    files <- list.files(folder, pattern = '\\.genes\\.sf$', full.names = TRUE,
                        ignore.case = TRUE)
    if (length(files) == 0L) {
      stop('No *.genes.sf files in ', folder,
           '. Check folder and file extensions.', call. = FALSE)
    }
    filenames <- basename(files)
    regex <- paste0('^PPMI-Phase[0-9]+-IR3\\.([0-9]+[A-Za-z]*)\\.', visit, '\\..*\\.genes\\.sf$')
    valid <- grepl(regex, filenames, ignore.case = TRUE)
    if (any(!valid)) {
      stop('Cannot extract PATNO from ', visit, ' file(s): ',
           paste(head(filenames[!valid], 10L), collapse = ', '), call. = FALSE)
    }
    patno <- sub(regex, '\\1', filenames, ignore.case = TRUE)
    if (anyDuplicated(patno)) {
      stop('Duplicate PATNO in ', visit, ': ',
           paste(unique(patno[duplicated(patno)]), collapse = ', '),
           '. Resolve duplicated samples first.', call. = FALSE)
    }
    data.frame(PATNO = patno, file_path = files,
               filename = filenames, stringsAsFactors = FALSE)
  }

  bl <- list_sf(bl_dir, 'BL')
  v02 <- list_sf(v02_dir, 'V02')
  names(bl)[2:3]  <- c('BL_file_path', 'BL_file')
  names(v02)[2:3] <- c('V02_file_path', 'V02_file')

  matched_df <- merge(bl, v02, by = 'PATNO', sort = TRUE)
  cat('BL samples:', nrow(bl), '\n')
  cat('V02 samples:', nrow(v02), '\n')
  cat('Common PATNO:', nrow(matched_df), '\n')
  cat('BL IDs with letter suffix:', sum(grepl('[A-Za-z]$', bl$PATNO)), '\n')
  cat('V02 IDs with letter suffix:', sum(grepl('[A-Za-z]$', v02$PATNO)), '\n')
  if (nrow(matched_df) == 0L) {
    stop('No common PATNO. Check sample names and visit folders.', call. = FALSE)
  }

  if (copy_matched_sf) {
    dir.create(out_bl_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_v02_dir, recursive = TRUE, showWarnings = FALSE)
    success_bl <- file.copy(matched_df$BL_file_path,
                            file.path(out_bl_dir, matched_df$BL_file),
                            overwrite = TRUE)
    success_v02 <- file.copy(matched_df$V02_file_path,
                             file.path(out_v02_dir, matched_df$V02_file),
                             overwrite = TRUE)
    if (!all(success_bl) || !all(success_v02)) {
      stop('Failed to copy some *.sf files. Check disk space and write permissions.',
           call. = FALSE)
    }
  }

  utils::write.csv(matched_df, output_csv, row.names = FALSE)
  cat('Saved:', output_csv, '\n')
  invisible(matched_df)
}

# Runs without relying on any objects from a previous script/session.
run_match_bl_v02()
