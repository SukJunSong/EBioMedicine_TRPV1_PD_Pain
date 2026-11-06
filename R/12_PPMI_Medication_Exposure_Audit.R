# Input/output paths are configured via environment variables: PD_TRPV1_ROOT, PD_TRPV1_CLINICAL_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)
if (!nzchar(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""))) stop("Set PD_TRPV1_CLINICAL_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 12_PPMI_Medication_Exposure_Audit.R
# Requires original project datasets; full R execution not validated here.
# PPMI drug exposure audit for existing PD RNA-seq cohort
# Inclusion rule: medication initiated in the SAME MONTH as a visit is included as
# possible exposure, not proven use before that visit. Both strict and inclusive flags saved.
# R >= 4.1. Required packages: readr, readxl, dplyr, tidyr, stringr, lubridate
pkgs <- c('readr','readxl','dplyr','tidyr','stringr','lubridate')
miss <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)) stop('Install: install.packages(c(',paste(sprintf('"%s"',miss),collapse=','),'))')
suppressPackageStartupMessages({library(readr);library(readxl);library(dplyr);library(tidyr);library(stringr);library(lubridate)})
root <- Sys.getenv("PD_TRPV1_ROOT", unset="")
clinical_dir <- Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset="")
rna_dir <- file.path(root,'PPMI')
out <- file.path(root,'Medication_BL_V02_Audit')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fil <- list(conmed=file.path(root,'Concomitant_Medication_Log_26Oct2026.csv'),
            ledd=file.path(root,'LEDD_Concomitant_Medication_Log_26Oct2026.csv'),
            pain=file.path(clinical_dir,'MDS-UPDRS_Part_I_Patient_Questionnaire_25Feb2025.csv'),
            rna=file.path(rna_dir,'hsa04750pd_cleaned.xlsx'))
for(nm in names(fil)) if(!file.exists(fil[[nm]])) stop('Missing: ',fil[[nm]])
writeout <- function(x,file) write_csv(x,file.path(out,file),na='')
# Input dates are primarily MM/YYYY, and visits are also MM/YYYY.
# Convert to first of month, do not infer unknown day.
month_date <- function(x) {
  s <- str_trim(as.character(x)); s[s %in% c('','NA','N/A','Unknown','UNK')] <- NA_character_
  suppressWarnings(as.Date(ifelse(str_detect(s,'^\\d{1,2}/\\d{4}$'),paste0('01/',s),
     ifelse(str_detect(s,'^\\d{1,2}/\\d{1,2}/\\d{4}$'),s,NA_character_)),format='%d/%m/%Y'))
}
clean_id <- function(x) str_trim(as.character(x))
# Original cohort: patient IDs in PD TRPV1 sheet (188 subjects in uploaded version)
pd <- read_excel(fil$rna,sheet='TRPV1',col_types='text')
ids <- pd %>% transmute(PATNO=clean_id(PATNO)) %>% distinct()
writeout(ids,'01_PD_RNA_PATNO.csv')
# Visit dates in the MDS-UPDRS questionnaire; source dates explicitly traceable.
vis <- read_csv(fil$pain,show_col_types=FALSE,col_types=cols(.default=col_character())) %>%
  transmute(PATNO=clean_id(PATNO), EVENT_ID=str_trim(EVENT_ID),
            visit_date_raw=INFODT, visit_month=month_date(INFODT),
            pain= suppressWarnings(as.numeric(NP1PAIN))) %>%
  semi_join(ids,by='PATNO') %>% filter(EVENT_ID %in% c('BL','V02'))
# Unexpected duplicate/discordant dates should be audited rather than silently choosing.
writeout(vis %>% count(PATNO,EVENT_ID) %>% filter(n>1),'QC_duplicate_visit_rows.csv')
visit_summary <- vis %>% group_by(PATNO,EVENT_ID) %>%
  summarise(visit_month=if(n_distinct(visit_month,na.rm=TRUE)==1) first(visit_month[!is.na(visit_month)]) else as.Date(NA),
            visit_date_n=n_distinct(visit_month,na.rm=TRUE),
            pain=if(n_distinct(pain,na.rm=TRUE)==1 && any(!is.na(pain))) first(pain[!is.na(pain)]) else NA_real_,
            .groups='drop') %>%
  pivot_wider(names_from=EVENT_ID,values_from=c(visit_month,visit_date_n,pain),names_sep='_')
people <- ids %>% left_join(visit_summary,by='PATNO') %>%
 mutate(BL_group=case_when(is.na(pain_BL)~'Unknown', pain_BL==0~'No Pain',pain_BL==1~'Low Pain',pain_BL>=2~'High Pain',TRUE~'Unknown'))
writeout(people,'02_PD_visit_dates_and_baseline_pain.csv')

# Medication records. Retain raw strings and distinguish named drugs and indication.
cm <- read_csv(fil$conmed,show_col_types=FALSE,col_types=cols(.default=col_character())) %>%
 transmute(PATNO=clean_id(PATNO),source='CONMED',drug=CMTRT,indication=CMINDC_TEXT,
           start_raw=STARTDT,stop_raw=STOPDT,ongoing=ONGOING,event_log=EVENT_ID,REC_ID) %>%
 semi_join(ids,by='PATNO')
ld <- read_csv(fil$ledd,show_col_types=FALSE,col_types=cols(.default=col_character())) %>%
 transmute(PATNO=clean_id(PATNO),source='LEDD',drug=LEDTRT,indication=NA_character_,
           start_raw=STARTDT,stop_raw=STOPDT,ongoing=NA_character_,event_log=EVENT_ID,REC_ID) %>%
 semi_join(ids,by='PATNO')
# candidate patterns are screening rules, NOT proof of the clinical indication.
analgesic_regex <- paste(c(
 'ACETAMINOPHEN','PARACETAMOL','TYLENOL','IBUPROFEN','ADVIL','MOTRIN',
 'NAPROXEN','ALEVE','DICLOFENAC','VOLTAREN','CELECOXIB','CELEBREX','MELOXICAM',
 'PIROXICAM','ETODOLAC','INDOMETHACIN','KETOROLAC','NABUMETONE','SULINDAC',
 'ASPIRIN','SALICYLATE','TRAMADOL','ULTRAM','CODEINE','HYDROCODONE','OXYCODONE',
 'MORPHINE','FENTANYL','BUPRENORPHINE','TAPENTADOL','METHADONE','HYDROMORPHONE',
 'GABAPENTIN','NEURONTIN','PREGABALIN','LYRICA','DULOXETINE','CYMBALTA',
 'AMITRIPTYLINE','NORTRIPTYLINE','CARBAMAZEPINE','OXCARBAZEPINE','LIDOCAINE',
 'CAPSAICIN','SUMATRIPTAN','RIZATRIPTAN','BUTALBITAL'
),collapse='|')
pain_indication_regex <- 'PAIN|HEADACHE|MIGRAINE|NEURALGIA|ARTHRALGIA|SCIATICA|BACK ACHE|BACKACHE'
# LEDD table represents PD drug exposure; CONMED entries included as supplementary checks.
pd_drug_regex <- 'LEVODOPA|CARBIDOPA|SINEMET|RYTARY|STALEVO|ENTACAPONE|COMTAN|OPICAPONE|RASAGILINE|AZILECT|SELEGILINE|SAFINAMIDE|PRAMIPEXOLE|MIRAPEX|ROPINIROLE|REQUIP|ROTIGOTINE|NEUPRO|APOMORPHINE|AMANTADINE|GOCOVRI|TRIHEXYPHENIDYL|BENZTROPINE|BROMOCRIPTINE|PERGOLIDE'

med <- bind_rows(cm,ld) %>% mutate(
 drug_upper=str_to_upper(coalesce(drug,'')),indication_upper=str_to_upper(coalesce(indication,'')),
 start_month=month_date(start_raw), stop_month=month_date(stop_raw),
 analgesic_candidate=(source=='CONMED' & (str_detect(drug_upper,analgesic_regex) | str_detect(indication_upper,pain_indication_regex))),
 pd_med_candidate=(source=='LEDD' | str_detect(drug_upper,pd_drug_regex))) %>%
 left_join(people %>% select(PATNO,BL_group,visit_month_BL,visit_month_V02),by='PATNO')
# Flag if a medication interval overlaps each visit month, and timing relative to that visit.
# Unknown end date is treated as potentially ongoing (not a certainty).
flag_at <- function(start,stop,visit) {
 !is.na(start) & !is.na(visit) & start<=visit & (is.na(stop)|stop>=visit)
}
med <- med %>% mutate(
 start_before_BL = !is.na(start_month) & !is.na(visit_month_BL) & start_month < visit_month_BL,
 start_same_BL = !is.na(start_month) & !is.na(visit_month_BL) & start_month == visit_month_BL,
 start_before_V02 = !is.na(start_month) & !is.na(visit_month_V02) & start_month < visit_month_V02,
 start_same_V02 = !is.na(start_month) & !is.na(visit_month_V02) & start_month == visit_month_V02,
 active_BL_possible = flag_at(start_month,stop_month,visit_month_BL),
 active_V02_possible = flag_at(start_month,stop_month,visit_month_V02),
 active_BL_started_strict = active_BL_possible & start_before_BL,
 active_V02_started_strict = active_V02_possible & start_before_V02,
 BL_to_V02_start_inclusive = !is.na(start_month) & !is.na(visit_month_BL) & !is.na(visit_month_V02) &
       start_month >= visit_month_BL & start_month <= visit_month_V02,
 medication_before_V02_inclusive = !is.na(start_month) & !is.na(visit_month_V02) & start_month <= visit_month_V02 &
       (is.na(stop_month) | stop_month >= visit_month_BL),
 dates_ambiguous = is.na(start_month) | is.na(visit_month_BL) | is.na(visit_month_V02) | is.na(stop_month)
)
writeout(med %>% filter(pd_med_candidate),'03_PD_drug_records_all.csv')
writeout(med %>% filter(analgesic_candidate),'04_analgesic_CANDIDATES_review.csv')
writeout(med %>% filter(start_same_BL | start_same_V02),'05_same_month_start_records.csv')
# Review the list of patients already taking PD meds before BL; preserve actual drug names.
writeout(med %>% filter(pd_med_candidate,active_BL_started_strict) %>%
 select(PATNO,BL_group,source,drug,start_raw,stop_raw,visit_month_BL,visit_month_V02,
        active_BL_started_strict,active_V02_possible,REC_ID),
 '06_PD_drugs_ACTIVE_at_BL_started_before_BL.csv')
# Patient summary: ANY matched drug; include same-month in inclusive group, separately reported.
flags <- med %>% group_by(PATNO) %>% summarise(
 pd_active_BL_strict=any(pd_med_candidate & active_BL_started_strict),
 pd_active_BL_inclusive=any(pd_med_candidate & active_BL_possible),
 pd_started_same_BL=any(pd_med_candidate & start_same_BL),
 pd_active_V02_strict=any(pd_med_candidate & active_V02_started_strict),
 pd_active_V02_inclusive=any(pd_med_candidate & active_V02_possible),
 pd_started_same_V02=any(pd_med_candidate & start_same_V02),
 pd_started_BL_to_V02_inclusive=any(pd_med_candidate & BL_to_V02_start_inclusive),
 pd_any_before_or_same_V02=any(pd_med_candidate & medication_before_V02_inclusive),
 analgesic_active_BL_inclusive=any(analgesic_candidate & active_BL_possible),
 analgesic_active_V02_inclusive=any(analgesic_candidate & active_V02_possible),
 analgesic_any_before_or_same_V02=any(analgesic_candidate & medication_before_V02_inclusive),
 analgesic_started_same_V02=any(analgesic_candidate & start_same_V02),
 .groups='drop')
flags <- people %>% left_join(flags,by='PATNO') %>%
 mutate(across(where(is.logical),~replace_na(.x,FALSE)))
writeout(flags,'07_patient_flags_188_inclusive_and_strict.csv')
writeout(flags %>% filter(pd_active_BL_inclusive | pd_started_same_BL),'08_PATNO_PD_drug_at_BL_including_same_month.csv')
writeout(flags %>% filter(pd_any_before_or_same_V02),'09_PATNO_PD_drug_before_or_same_V02.csv')
writeout(flags %>% filter(analgesic_any_before_or_same_V02),'10_PATNO_analgesic_candidates_before_or_same_V02.csv')
summary <- flags %>% group_by(BL_group) %>% summarise(
 n=n(),n_missing_visit_BL=sum(is.na(visit_month_BL)),n_missing_visit_V02=sum(is.na(visit_month_V02)),
 pd_active_BL_strict=sum(pd_active_BL_strict),pd_active_BL_inclusive=sum(pd_active_BL_inclusive),
 pd_same_BL=sum(pd_started_same_BL),pd_before_or_same_V02=sum(pd_any_before_or_same_V02),
 pd_same_V02=sum(pd_started_same_V02),
 analgesic_before_or_same_V02=sum(analgesic_any_before_or_same_V02),
 analgesic_same_V02=sum(analgesic_started_same_V02),.groups='drop')
writeout(bind_rows(summary,summary %>% ungroup() %>% summarise(BL_group='ALL',across(where(is.numeric),sum))),
         '11_group_summary_counts.csv')
writeLines(c(
 'Same-month starts are INCLUDED as potential medication exposures, not proven before assessment.',
 'Visit/medication dates are often MM/YYYY; unknown day cannot establish event order.',
 'BL-exposure checks interval overlap: start <= BL month and stop missing or >= BL month.',
 'PD drugs include LEDD rows and heuristic CONMED PD-drug names. Inspect detailed records.',
 'Analgesic candidates are intentionally broad and REQUIRE medication indication review.',
 'Missing STOPDT is treated as possibly ongoing, not definitely ongoing.',
 'Records with missing STARTDT cannot be classified as exposed; review manually.',
 'For sensitivity analyses, report strict and inclusive scenarios separately.',
 'The 188-ID set is obtained from the PD TRPV1 sheet, not assumed to represent the earlier 255-subject cohort.'
),file.path(out,'README_important.txt'))
cat('\nSaved to: ',out,'\n',sep='')
