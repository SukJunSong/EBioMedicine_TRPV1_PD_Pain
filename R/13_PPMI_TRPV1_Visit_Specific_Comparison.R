# Input/output paths are configured via environment variables: PD_TRPV1_ROOT, PD_TRPV1_CLINICAL_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)
if (!nzchar(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""))) stop("Set PD_TRPV1_CLINICAL_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 13_PPMI_TRPV1_Visit_Specific_Comparison.R
# Requires original project datasets; full R execution not validated here.
# Reproduce visit-specific TRPV1 comparison: HC vs PD No/Low/High Pain.
# Sensitivity analysis uses same visit medication flags, never BL pain group at V02.
# IMPORTANT: This RNA workbook has PD 188 and HC 130, not historic 255/170.
required <- c('readxl','readr','dplyr','tidyr','broom')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install missing packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(readxl);library(readr);library(dplyr);library(tidyr);library(broom)})
root <- Sys.getenv("PD_TRPV1_ROOT", unset="")
rna <- file.path(root,'PPMI')
clin <- Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset="")
out <- file.path(root,'Medication_TRPV1_CORRECTED')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
find_file <- function(n) {p <- file.path(rna,n);if(!file.exists(p))p<-file.path(root,n);if(!file.exists(p))stop('Missing: ',n);p}
num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
id <- function(x)trimws(as.character(x))
read_sheet <- function(fname,sheet,cohort) {
  x<-read_excel(find_file(fname),sheet=sheet,col_types='text');names(x)<-trimws(names(x))
  x%>%transmute(PATNO=id(PATNO),cohort=cohort,BL=num(Baseline),V02=num(`6 months`))
}
pd_gene<-read_sheet('hsa04750pd_cleaned.xlsx','TRPV1','PD')
hc_gene<-read_sheet('hsa04750ctrl_cleaned.xlsx','TRPV1','HC')
pd_pain<-read_sheet('hsa04750pd_cleaned.xlsx','pain','PD')%>%select(PATNO,pain_BL=BL,pain_V02=V02)
if(anyDuplicated(pd_gene$PATNO)||anyDuplicated(hc_gene$PATNO)||anyDuplicated(pd_pain$PATNO))stop('Duplicate RNA PATNO')
meta<-read_excel(find_file('metaData BL.xlsx'),col_types='text')%>%
  transmute(PATNO=id(PATNO),sex=trimws(as.character(GENDER)))%>%
  group_by(PATNO)%>%summarise(sex=if(n_distinct(sex[!is.na(sex)&sex!=''])==1)first(sex[!is.na(sex)&sex!='']) else NA_character_,.groups='drop')
age_file<-file.path(clin,'Age_at_visit_20Mar2025.csv')
if(!file.exists(age_file))stop('Missing ',age_file)
ages<-read_csv(age_file,col_types=cols(.default=col_character()),show_col_types=FALSE)%>%
  transmute(PATNO=id(PATNO),visit=trimws(EVENT_ID),age=num(AGE_AT_VISIT))%>%
  filter(visit%in%c('BL','V02'))%>%group_by(PATNO,visit)%>%
  summarise(age=if(n_distinct(age[!is.na(age)])==1)first(age[!is.na(age)])else NA_real_,.groups='drop')
flagpath<-file.path(root,'Medication_BL_V02_Audit','07_patient_flags_188_inclusive_and_strict.csv')
if(!file.exists(flagpath))stop('Run medication audit first: ',flagpath)
flags<-read_csv(flagpath,col_types=cols(.default=col_character()),show_col_types=FALSE)%>%
  transmute(PATNO=id(PATNO),
    analgesic_BL=toupper(analgesic_active_BL_inclusive)=='TRUE',
    analgesic_V02=toupper(analgesic_active_V02_inclusive)=='TRUE',
    pd_drug_BL=toupper(pd_active_BL_inclusive)=='TRUE',
    pd_drug_V02=toupper(pd_active_V02_inclusive)=='TRUE')
pd_base<-pd_gene%>%left_join(pd_pain,by='PATNO')%>%left_join(flags,by='PATNO')
allbase<-bind_rows(pd_base,hc_gene)
long<-bind_rows(lapply(c('BL','V02'),function(v){
  z<-allbase%>%transmute(PATNO,cohort,visit=v,TRPV1=.data[[v]],pain=if(v=='BL')pain_BL else pain_V02,
    analgesic=if(v=='BL')analgesic_BL else analgesic_V02,
    pd_drug=if(v=='BL')pd_drug_BL else pd_drug_V02)
  z
}))%>%left_join(ages,by=c('PATNO','visit'))%>%left_join(meta,by='PATNO')%>%
  mutate(group=case_when(cohort=='HC'~'HC',is.na(pain)~NA_character_,pain==0~'No Pain',pain==1~'Low Pain',pain>=2~'High Pain'),
    group=factor(group,levels=c('HC','No Pain','Low Pain','High Pain')),
    sex=if_else(sex%in%c('Male','Female'),sex,NA_character_))
write_csv(long,file.path(out,'01_patient_visit_comparison_data.csv'))
qc<-long%>%group_by(visit,group)%>%summarise(n=n(),n_gene=sum(!is.na(TRPV1)),n_age_sex=sum(!is.na(TRPV1)&!is.na(age)&!is.na(sex)),
 mean_TRPV1=mean(TRPV1,na.rm=TRUE),n_analgesic_possible=sum(analgesic%in%TRUE,na.rm=TRUE),.groups='drop')
write_csv(qc,file.path(out,'02_QC_n_by_VISIT_specific_group.csv'))
# HC drug data aren't available from the PD-only flag audit. Drug-adjusted 4-group models
# are deliberately NOT run with unknown HC exposure coded as zero.
models<-list()
for(v in c('BL','V02')){
 d<-long%>%filter(visit==v,!is.na(TRPV1),!is.na(group))
 for(spec in c('raw','age_sex')){
  dd<-if(spec=='age_sex') filter(d,!is.na(age),!is.na(sex)) else d
  if(spec=='age_sex')dd<-dd%>%mutate(sex=factor(sex))
  f<-if(spec=='raw') TRPV1~group else TRPV1~group+age+sex
  fit<-lm(f,data=dd)
  res<-broom::tidy(fit,conf.int=TRUE)%>%mutate(visit=v,model=spec,n=nobs(fit))
  models[[paste(v,spec)]]<-res
 }
}
model_tab<-bind_rows(models)%>%group_by(visit,model)%>%
 mutate(FDR_BH_3_vs_HC=if_else(grepl('^group',term),p.adjust(ifelse(grepl('^group',term),p.value,NA_real_),method='BH'),NA_real_))%>%ungroup()
write_csv(model_tab,file.path(out,'03_HC_reference_group_models.csv'))
# Separate PD-only medication sensitivity (same visit pain group)
pdmodels<-list()
for(v in c('BL','V02')){
 d<-long%>%filter(visit==v,cohort=='PD',!is.na(group),!is.na(TRPV1),!is.na(age),!is.na(sex))%>%
  mutate(group=factor(as.character(group),levels=c('No Pain','Low Pain','High Pain')),analgesic=as.logical(analgesic),pd_drug=as.logical(pd_drug),sex=factor(sex))
 for(spec in c('age_sex','plus_analgesic','plus_analgesic_PDdrug')){
  formula<-switch(spec,age_sex=TRPV1~group+age+sex,
      plus_analgesic=TRPV1~group+age+sex+analgesic,
      plus_analgesic_PDdrug=TRPV1~group+age+sex+analgesic+pd_drug)
  fit<-lm(formula,data=d)
  pdmodels[[paste(v,spec)]]<-broom::tidy(fit,conf.int=TRUE)%>%mutate(visit=v,model=spec,n=nobs(fit))
 }
}
write_csv(bind_rows(pdmodels),file.path(out,'04_PD_ONLY_visit_pain_group_medication_models.csv'))
writeLines(c('Critical: V02 pain groups use V02 pain scores, NOT BL pain groups.',
 '03: HC reference group models vs PD No/Low/High. 04: PD-only within-group sensitivity; DIFFERENT COMPARISON.',
 'Historical published significant result may use a different n=255/170 cohort, normalization or model.',
 'Cannot adjust HC/PD combined model for medication until HC medication exposures are defined.',
 'Same-month starts counted as POSSIBLE drug exposure, not confirmed ingestion before blood draw.'),file.path(out,'README.txt'))
message('Saved to ',out)
