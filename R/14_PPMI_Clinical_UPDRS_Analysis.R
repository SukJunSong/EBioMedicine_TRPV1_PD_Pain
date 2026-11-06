# Input/output paths are configured via environment variables: PD_TRPV1_ROOT, PD_TRPV1_CLINICAL_DIR.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)
if (!nzchar(Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset=""))) stop("Set PD_TRPV1_CLINICAL_DIR before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 14_PPMI_Clinical_UPDRS_Analysis.R
# Requires original project datasets; full R execution not validated here.
# eBioMedicine revision — PPMI, improved age/sex-adjusted clinical and TRPV1 analyses
# Changes from prior version:
# 1) Part III duplicate visits are NOT all dropped; one row per PATNO/EVENT_ID is chosen, preferring OFF state.
# 2) TRPV1 correlation figures are generated.
# 3) Baseline pain colours swapped: No Pain = blue, Low Pain = green, High Pain = red.

packages <- c('readxl','dplyr','tidyr','ggplot2','readr','purrr','tibble','lme4','lmerTest','broom','emmeans')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install packages first: install.packages(c(', paste(sprintf('"%s"', missing), collapse=','), '))')
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr); library(purrr); library(tibble)
})

clinical_dir <- Sys.getenv("PD_TRPV1_CLINICAL_DIR", unset="")
rna_dir <- file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), 'PPMI')
out_dir <- file.path(rna_dir,'Revision_Analysis_R_adjusted_v2')
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

f <- list(
  p1   = file.path(clinical_dir,'MDS-UPDRS_Part_I_21Mar2025.csv'),
  p1p  = file.path(clinical_dir,'MDS-UPDRS_Part_I_Patient_Questionnaire_25Feb2025.csv'),
  p2   = file.path(clinical_dir,'MDS_UPDRS_Part_II__Patient_Questionnaire_21Mar2025.csv'),
  p3   = file.path(clinical_dir,'MDS-UPDRS_Part_III_21Mar2025.csv'),
  age  = file.path(clinical_dir,'Age_at_visit_20Mar2025.csv'),
  meta = file.path(rna_dir,'metaData BL.xlsx'),
  pd   = file.path(rna_dir,'hsa04750pd_cleaned.xlsx'),
  ctrl = file.path(rna_dir,'hsa04750ctrl_cleaned.xlsx')
)
if(any(!file.exists(unlist(f)))) stop('Missing inputs:\n', paste(unlist(f)[!file.exists(unlist(f))], collapse='\n'))

visit_map  <- c('Baseline'='BL','6 months'='V02','12 months'='V04','24 months'='V08','36 months'='V12')
months_map <- c('BL'=0,'V02'=6,'V04'=12,'V08'=24,'V12'=36)
key <- c('PATNO','EVENT_ID')
pain_colours <- c('No Pain'='#619CFF', 'Low Pain'='#00BA38', 'High Pain'='#F8766D')

clean_id <- function(x) trimws(as.character(x))
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
write_csv2 <- function(dat, name) readr::write_csv(dat, file.path(out_dir, name), na='')

#------------------------------
# RNA expression sheets
#------------------------------
read_gene <- function(file, sheet, cohort) {
  x <- readxl::read_excel(file, sheet=sheet, col_types='text')
  names(x) <- trimws(names(x))
  needed <- c('PATNO', names(visit_map))
  if(!all(needed %in% names(x))) stop('Missing expected columns in sheet ', sheet)
  x %>%
    mutate(PATNO = clean_id(PATNO)) %>%
    pivot_longer(cols = all_of(names(visit_map)), names_to='visit_label', values_to='raw_value') %>%
    mutate(EVENT_ID = unname(visit_map[visit_label]),
           months   = unname(as.numeric(months_map[EVENT_ID])),
           value    = safe_num(trimws(raw_value)),
           cohort   = cohort)
}

expr_raw <- bind_rows(read_gene(f$pd, 'TRPV1', 'PD'), read_gene(f$ctrl, 'TRPV1', 'HC'))
write_csv2(expr_raw %>% filter(!is.na(raw_value), nzchar(trimws(raw_value)), is.na(value)) %>% select(PATNO, EVENT_ID, cohort, raw_value), 'QC_TRPV1_nonnumeric.csv')
expr <- expr_raw %>% select(PATNO, EVENT_ID, months, cohort, TRPV1=value)
if(anyDuplicated(expr[c('PATNO','EVENT_ID','cohort')])) stop('RNA expression has duplicate patient/visit/cohort rows')

pain_gene <- bind_rows(read_gene(f$pd, 'pain', 'PD'), read_gene(f$ctrl, 'pain', 'HC')) %>%
  select(PATNO, EVENT_ID, cohort, pain_gene=value)

#------------------------------
# Age / sex metadata
#------------------------------
a_raw <- read_csv(f$age, show_col_types=FALSE,
                  col_types=cols(PATNO=col_character(), EVENT_ID=col_character(), AGE_AT_VISIT=col_double())) %>%
  mutate(PATNO=clean_id(PATNO), EVENT_ID=trimws(EVENT_ID))
write_csv2(a_raw %>% count(PATNO, EVENT_ID, name='n_rows') %>% filter(n_rows > 1), 'QC_Age_duplicate_visits.csv')
age <- a_raw %>%
  group_by(PATNO, EVENT_ID) %>%
  summarise(age_visit = if(n_distinct(AGE_AT_VISIT, na.rm=TRUE)==1) first(AGE_AT_VISIT) else NA_real_,
            age_conflict = n_distinct(AGE_AT_VISIT, na.rm=TRUE) > 1,
            .groups='drop')

meta_raw <- readxl::read_excel(f$meta, sheet=1) %>%
  transmute(PATNO = clean_id(PATNO),
            sex_raw = trimws(as.character(GENDER)),
            metadata_diagnosis = as.character(DIAGNOSIS))
write_csv2(meta_raw %>% count(PATNO) %>% filter(n > 1), 'QC_Metadata_duplicate_PATNO.csv')
meta <- meta_raw %>%
  group_by(PATNO) %>%
  summarise(sex = if(n_distinct(sex_raw[!is.na(sex_raw)]) == 1) first(sex_raw[!is.na(sex_raw)]) else NA_character_,
            metadata_diagnosis = first(metadata_diagnosis), .groups='drop') %>%
  mutate(sex = if_else(sex %in% c('Male','Female'), sex, NA_character_),
         sex = factor(sex, levels=c('Female','Male')))

#------------------------------
# Clinical data
#------------------------------
p1  <- read_csv(f$p1,  show_col_types=FALSE) %>% transmute(PATNO=clean_id(PATNO), EVENT_ID, NP1RTOT=safe_num(NP1RTOT))
p1p <- read_csv(f$p1p, show_col_types=FALSE) %>% transmute(PATNO=clean_id(PATNO), EVENT_ID, NP1PAIN=safe_num(NP1PAIN), NP1PTOT=safe_num(NP1PTOT))
p2  <- read_csv(f$p2,  show_col_types=FALSE) %>% transmute(PATNO=clean_id(PATNO), EVENT_ID, PartII=safe_num(NP2PTOT))
for(z in list(p1, p1p, p2)) if(anyDuplicated(z[key])) stop('Unexpected duplicate keys in Part I / Part II file')

p3_raw <- read_csv(f$p3, show_col_types=FALSE) %>%
  transmute(PATNO=clean_id(PATNO), EVENT_ID,
            PartIII=safe_num(NP3TOT),
            PDSTATE=trimws(as.character(PDSTATE)),
            PDTRTMNT=safe_num(PDTRTMNT),
            PDMEDYN=safe_num(PDMEDYN),
            OFFEXAM=safe_num(OFFEXAM),
            ONEXAM=safe_num(ONEXAM),
            INFODT=as.character(INFODT))

# Choose one Part III row per patient/visit, preferring OFF state.
choose_part3_row <- function(df) {
  if(nrow(df) == 1) {
    df$choice_reason <- 'single_row'
    return(df)
  }
  off <- df %>% filter(toupper(PDSTATE) == 'OFF' | OFFEXAM == 1)
  if(nrow(off) >= 1) {
    out <- off %>% slice(1)
    out$choice_reason <- 'preferred_OFF'
    return(out)
  }
  on <- df %>% filter(toupper(PDSTATE) == 'ON' | ONEXAM == 1)
  if(nrow(on) >= 1) {
    out <- on %>% slice(1)
    out$choice_reason <- 'fallback_ON'
    return(out)
  }
  out <- df %>% slice(1)
  out$choice_reason <- 'fallback_first'
  out
}

p3_dup <- p3_raw %>% count(PATNO, EVENT_ID, name='n_rows') %>% filter(n_rows > 1)
write_csv2(p3_raw %>% semi_join(p3_dup, by=key), 'QC_PartIII_duplicate_rows_all.csv')

p3 <- p3_raw %>%
  group_by(PATNO, EVENT_ID) %>%
  group_modify(~ choose_part3_row(.x)) %>%
  ungroup()
write_csv2(p3 %>% semi_join(p3_dup, by=key) %>% select(PATNO, EVENT_ID, PartIII, PDSTATE, PDTRTMNT, PDMEDYN, OFFEXAM, ONEXAM, choice_reason), 'QC_PartIII_duplicate_rows_selected.csv')
if(anyDuplicated(p3[key])) stop('Part III still has duplicate PATNO/EVENT_ID after selection')

clinical <- p1p %>%
  full_join(p1, by=key) %>%
  full_join(p2, by=key) %>%
  full_join(p3 %>% select(-choice_reason), by=key) %>%
  left_join(age, by=key) %>%
  left_join(meta, by='PATNO') %>%
  mutate(PartI = if_else(!is.na(NP1RTOT) & !is.na(NP1PTOT), NP1RTOT + NP1PTOT, NA_real_),
         PartI_noPain = PartI - NP1PAIN,
         Total = if_else(!is.na(PartI) & !is.na(PartII) & !is.na(PartIII), PartI + PartII + PartIII, NA_real_),
         Total_noPain = Total - NP1PAIN)
if(anyDuplicated(clinical[key])) stop('Clinical merge produced duplicate patient/visit rows')

#------------------------------
# Merge analysis dataset
#------------------------------
analysis <- expr %>%
  left_join(clinical, by=key) %>%
  left_join(pain_gene, by=c(key, 'cohort')) %>%
  mutate(PainGroup = case_when(
    is.na(NP1PAIN) ~ NA_character_,
    NP1PAIN == 0   ~ 'No Pain',
    NP1PAIN == 1   ~ 'Low Pain',
    NP1PAIN >= 2   ~ 'High Pain'
  ),
  PainGroup = factor(PainGroup, levels=c('No Pain','Low Pain','High Pain')))
write_csv2(analysis, '00_merged_all.csv')
write_csv2(analysis %>% filter(cohort=='PD') %>% distinct(PATNO, cohort, sex, metadata_diagnosis) %>% filter(is.na(sex)), 'QC_PD_subjects_missing_sex.csv')
write_csv2(analysis %>% filter(!is.na(pain_gene), !is.na(NP1PAIN), pain_gene != NP1PAIN) %>% select(PATNO, EVENT_ID, cohort, pain_gene, NP1PAIN), 'QC_gene_vs_clinical_pain_disagreement.csv')

# Overall completeness
qc_visit <- analysis %>%
  group_by(cohort, EVENT_ID, months) %>%
  summarise(n_unique = n_distinct(PATNO),
            n_trpv1 = sum(!is.na(TRPV1)),
            n_pain  = sum(!is.na(NP1PAIN)),
            n_sex   = sum(!is.na(sex)),
            n_age   = sum(!is.na(age_visit)),
            n_age_sex_complete = sum(!is.na(age_visit) & !is.na(sex)),
            n_PartIII = sum(!is.na(PartIII)), .groups='drop')
write_csv2(qc_visit, 'QC_completeness_by_visit.csv')

# Baseline pain groups for PD cohort
base <- analysis %>% filter(cohort=='PD', EVENT_ID=='BL') %>%
  transmute(PATNO,
            baseline_pain = NP1PAIN,
            age_bl = age_visit,
            sex_bl = sex,
            BL_group = case_when(
              is.na(NP1PAIN) ~ NA_character_,
              NP1PAIN == 0 ~ 'No Pain',
              NP1PAIN == 1 ~ 'Low Pain',
              NP1PAIN >= 2 ~ 'High Pain'
            )) %>%
  mutate(BL_group = factor(BL_group, levels=c('No Pain','Low Pain','High Pain')))

write_csv2(base %>% count(BL_group, name='n_PD_baseline_group'), 'QC_PD_baseline_pain_group_counts.csv')

# Group x visit counts for each outcome
vars <- c('PartI_noPain','PartII','PartIII','Total_noPain')
qc_group_counts <- analysis %>%
  filter(cohort=='PD') %>%
  left_join(base %>% select(PATNO, BL_group), by='PATNO') %>%
  pivot_longer(cols=all_of(vars), names_to='outcome', values_to='score') %>%
  group_by(outcome, EVENT_ID, BL_group) %>%
  summarise(n = n_distinct(PATNO[!is.na(score)]), .groups='drop') %>%
  arrange(outcome, EVENT_ID, BL_group)
write_csv2(qc_group_counts, 'QC_PD_group_counts_by_visit_outcome.csv')

#------------------------------
# 01. Spearman correlations
#------------------------------
spcor <- function(dat, x, y) {
  d <- dat %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  if(nrow(d) < 8 || n_distinct(d[[x]]) < 2 || n_distinct(d[[y]]) < 2) {
    return(tibble(n=nrow(d), rho=NA_real_, p=NA_real_))
  }
  ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method='spearman', exact=FALSE))
  tibble(n=nrow(d), rho=unname(ct$estimate), p=ct$p.value)
}

cor_results <- map_dfr(c('pain','TRPV1'), function(predictor) {
  x <- if(predictor=='pain') 'NP1PAIN' else 'TRPV1'
  targets <- if(predictor=='TRPV1') c('NP1PAIN', vars) else vars
  map_dfr(names(months_map), function(event) {
    map_dfr(targets, function(v) {
      bind_cols(tibble(predictor=predictor, EVENT_ID=event, outcome=v),
                spcor(filter(analysis, cohort=='PD', EVENT_ID==event), x, v))
    })
  })
}) %>% group_by(predictor, EVENT_ID) %>% mutate(p_BH = p.adjust(p, 'BH')) %>% ungroup()
write_csv2(cor_results, '01_within_visit_Spearman.csv')

#------------------------------
# 02. Age/sex-adjusted cross-sectional linear models
#------------------------------
fit_cross <- function(dat, outcome, predictor) {
  x <- if(predictor=='pain') 'NP1PAIN' else 'TRPV1'
  d <- dat %>% filter(is.finite(.data[[outcome]]), is.finite(.data[[x]]), is.finite(age_visit), !is.na(sex))
  if(nrow(d) < 15 || n_distinct(d[[x]]) < 2 || n_distinct(d$sex) < 2) {
    return(tibble(n=nrow(d), beta=NA_real_, se=NA_real_, ci_low=NA_real_, ci_high=NA_real_, p=NA_real_))
  }
  model <- lm(reformulate(c(x, 'age_visit', 'sex'), response=outcome), data=d)
  co <- summary(model)$coefficients
  if(!x %in% rownames(co)) return(tibble(n=nrow(d), beta=NA_real_, se=NA_real_, ci_low=NA_real_, ci_high=NA_real_, p=NA_real_))
  ci <- confint(model, parm=x)
  tibble(n=nrow(d), beta=unname(co[x,1]), se=unname(co[x,2]), ci_low=ci[1], ci_high=ci[2], p=unname(co[x,4]))
}

cross_results <- map_dfr(c('pain','TRPV1'), function(predictor) {
  targets <- if(predictor=='TRPV1') c('NP1PAIN', vars) else vars
  map_dfr(names(months_map), function(event) {
    map_dfr(targets, function(v) {
      bind_cols(tibble(predictor=predictor, EVENT_ID=event, outcome=v),
                fit_cross(filter(analysis, cohort=='PD', EVENT_ID==event), v, predictor))
    })
  })
}) %>% group_by(predictor, EVENT_ID) %>% mutate(p_BH = p.adjust(p, 'BH')) %>% ungroup()
write_csv2(cross_results, '02_age_sex_adjusted_linear_models.csv')

#------------------------------
# 03. Longitudinal prognosis by baseline pain group
#------------------------------
long <- analysis %>%
  filter(cohort=='PD') %>%
  select(-PainGroup, -sex) %>%
  left_join(base, by='PATNO') %>%
  filter(!is.na(BL_group)) %>%
  mutate(visit = factor(EVENT_ID, levels=names(months_map)),
         PATNO = factor(PATNO))
write_csv2(long, '03_longitudinal_dataset.csv')

fit_lmm <- function(v) {
  d <- long %>% filter(!is.na(.data[[v]]), !is.na(age_bl), !is.na(sex_bl), !is.na(visit)) %>% droplevels()
  if(nrow(d) < 30 || n_distinct(d$PATNO) < 15 || n_distinct(d$BL_group) < 2 || n_distinct(d$visit) < 2 || n_distinct(d$sex_bl) < 2) {
    return(tibble(outcome=v, term='insufficient_data', estimate=NA_real_, std_error=NA_real_, p=NA_real_, n=nrow(d), n_people=n_distinct(d$PATNO)))
  }
  form <- as.formula(paste0(v, ' ~ visit * BL_group + age_bl + sex_bl + (1|PATNO)'))
  fit <- tryCatch(lmerTest::lmer(form, data=d, REML=FALSE), error=function(e) NULL)
  if(is.null(fit)) {
    return(tibble(outcome=v, term='model_failed', estimate=NA_real_, std_error=NA_real_, p=NA_real_, n=nrow(d), n_people=n_distinct(d$PATNO)))
  }
  co <- as.data.frame(coef(summary(fit))); co$term <- rownames(co)
  write_csv2(tibble(outcome=v, term=co$term, estimate=co$Estimate, std_error=co[['Std. Error']], p=co[['Pr(>|t|)']]), paste0('03_LMM_', v, '.csv'))
  test <- as.data.frame(anova(fit, type=3)); test$term <- rownames(test)
  write_csv2(bind_cols(tibble(outcome=v), tibble::as_tibble(test)), paste0('03_LMM_typeIII_', v, '.csv'))
  tibble(outcome=v, term=co$term, estimate=co$Estimate, std_error=co[['Std. Error']], p=co[['Pr(>|t|)']], n=nrow(d), n_people=n_distinct(d$PATNO))
}
models <- bind_rows(lapply(vars, fit_lmm))
write_csv2(models, '03_all_LMM_coefficients.csv')

# Descriptive trajectories
summ <- long %>%
  pivot_longer(cols=all_of(vars), names_to='outcome', values_to='score') %>%
  filter(!is.na(score)) %>%
  group_by(outcome, BL_group, months) %>%
  summarise(n=n_distinct(PATNO), mean=mean(score), sd=sd(score), .groups='drop') %>%
  mutate(se=sd/sqrt(n),
         ci_low=if_else(n>1, mean - qt(.975, n-1)*se, NA_real_),
         ci_high=if_else(n>1, mean + qt(.975, n-1)*se, NA_real_))
write_csv2(summ, '03_descriptive_trajectories_95CI.csv')

for(v in vars) {
  p <- ggplot(filter(summ, outcome==v), aes(months, mean, colour=BL_group, group=BL_group)) +
    geom_line(linewidth=.8, na.rm=TRUE) +
    geom_point(size=2.7, na.rm=TRUE) +
    geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.7, na.rm=TRUE) +
    scale_color_manual(values=pain_colours, drop=FALSE) +
    scale_x_continuous(breaks=c(0,6,12,24,36)) +
    labs(title=paste('Baseline pain and', v),
         subtitle='Descriptive means and 95% CIs (not model-adjusted)',
         x='Months', y='MDS-UPDRS score', colour='Baseline pain') +
    theme_classic(base_size=12)
  ggsave(file.path(out_dir, paste0('03_trajectory_', v, '.png')), p, width=7, height=4.8, dpi=300)
}

#------------------------------
# 04. TRPV1 correlation figures
#------------------------------
make_label_df <- function(event) {
  cor_results %>%
    filter(predictor=='TRPV1', EVENT_ID==event) %>%
    mutate(label = paste0('n=', n, '\n', 'rho=', ifelse(is.na(rho), 'NA', sprintf('%.2f', rho)), '\n', 'BH p=', ifelse(is.na(p_BH), 'NA', formatC(p_BH, format='e', digits=2)))) %>%
    select(outcome, label)
}

for(event in names(months_map)) {
  plot_dat <- analysis %>%
    filter(cohort=='PD', EVENT_ID==event) %>%
    transmute(TRPV1, NP1PAIN, PartI_noPain, PartII, PartIII, Total_noPain) %>%
    pivot_longer(cols=c(NP1PAIN, PartI_noPain, PartII, PartIII, Total_noPain), names_to='outcome', values_to='score') %>%
    filter(is.finite(TRPV1), is.finite(score))
  lab_dat <- make_label_df(event)
  p <- ggplot(plot_dat, aes(TRPV1, score)) +
    geom_point(alpha=.65, size=1.8) +
    geom_smooth(method='lm', se=FALSE, linewidth=.6, colour='black') +
    facet_wrap(~ outcome, scales='free_y') +
    geom_text(data=lab_dat, aes(x=-Inf, y=Inf, label=label), inherit.aes=FALSE, hjust=-0.1, vjust=1.1, size=3) +
    labs(title=paste0('TRPV1 correlations at ', event), x='TRPV1 expression', y='Clinical score') +
    theme_classic(base_size=11)
  ggsave(file.path(out_dir, paste0('04_TRPV1_correlations_', event, '.png')), p, width=10, height=6.5, dpi=300)
}

# Combined overview plot for TRPV1 vs pain/UPDRS coefficients by visit
coef_plot_dat <- cross_results %>%
  filter(predictor=='TRPV1') %>%
  mutate(EVENT_ID=factor(EVENT_ID, levels=names(months_map)),
         outcome=factor(outcome, levels=c('NP1PAIN', vars)))
if(nrow(coef_plot_dat) > 0) {
  pcoef <- ggplot(coef_plot_dat, aes(EVENT_ID, beta, ymin=ci_low, ymax=ci_high)) +
    geom_hline(yintercept=0, linetype=2, colour='grey50') +
    geom_pointrange() +
    facet_wrap(~ outcome, scales='free_y') +
    labs(title='Age/sex-adjusted association of TRPV1 with pain and MDS-UPDRS measures',
         x='Visit', y='Beta (95% CI)') +
    theme_classic(base_size=11)
  ggsave(file.path(out_dir, '04_TRPV1_adjusted_coefficients_overview.png'), pcoef, width=10, height=6.5, dpi=300)
}

writeLines(c(
  'Main fix versus the prior script: Part III duplicate ON/OFF visits are no longer all excluded.',
  'One Part III row per PATNO/EVENT_ID is selected, preferring OFF state; see QC_PartIII_duplicate_rows_selected.csv.',
  'Because of that change, PartIII and Total_noPain trajectory plots should be re-interpreted using the new outputs only.',
  'Trajectory plots are descriptive means with 95% CIs, not model-adjusted estimates.',
  'TRPV1 correlation plots are now generated for each visit.',
  'Pain-group colours: No Pain = blue, Low Pain = green, High Pain = red.',
  'The uploaded RNA sheets contain PD n=188 and HC n=130, not the original manuscript PD n=255 and HC n=170.',
  'Analgesic medication has not been adjusted here.',
  'Age/sex-adjusted cross-sectional models use complete cases only.',
  'If publication-ready adjusted trajectory plots are needed, use emmeans/marginal means from the LMM rather than the descriptive plots.'
), file.path(out_dir, 'README_v2_notes.txt'))

cat('\nFinished. Files saved to: ', out_dir, '\n', sep='')


# =============================================================
# 05. Publication-sized adjusted beta heatmaps
# =============================================================
# Age/sex-adjusted beta, BH-corrected p stars. These are NOT Spearman rho.
# One standardized colour scale per panel around 0; beta units differ
# across clinical outcomes, so cell magnitudes are not comparable across rows.
star_bh <- function(p) dplyr::case_when(
  is.na(p) ~ '', p < .001 ~ '***', p < .01 ~ '**', p < .05 ~ '*', TRUE ~ ''
)
make_beta_heatmap <- function(predictor_key, name, subtitle) {
  labels <- c(NP1PAIN='Pain score',
              PartI_noPain='MDS-UPDRS Part I (pain excluded)',
              PartII='MDS-UPDRS Part II', PartIII='MDS-UPDRS Part III',
              Total_noPain='MDS-UPDRS Total (pain excluded)')
  visits <- c(BL='Baseline', V02='6 months', V04='12 months', V08='24 months', V12='36 months')
  dat <- cross_results %>% filter(predictor == predictor_key) %>%
    mutate(outcome_label=unname(labels[outcome]),
           visit_label=unname(visits[EVENT_ID]),
           cell=ifelse(is.finite(beta), paste0(sprintf('%.2f', beta), star_bh(p_BH)), 'NA'),
           outcome_label=factor(outcome_label, levels=rev(unname(labels[if(predictor_key=='TRPV1') c('NP1PAIN',vars) else vars]))),
           visit_label=factor(visit_label, levels=unname(visits)))
  limit <- max(abs(dat$beta), na.rm=TRUE)
  if(!is.finite(limit) || limit==0) limit <- 1
  p <- ggplot(dat, aes(visit_label, outcome_label, fill=beta)) +
    geom_tile(colour='white', linewidth=1) +
    geom_text(aes(label=cell), size=4.3, fontface='bold') +
    scale_fill_gradient2(low='#648FFF', mid='white', high='#FA8072',
                         midpoint=0, limits=c(-limit,limit), na.value='grey90',
                         name='Adjusted beta') +
    coord_cartesian(clip='off') +
    labs(title=name, subtitle=subtitle, x=NULL, y=NULL,
         caption='Cell labels: age-/sex-adjusted beta. * BH p < .05; ** p < .01; *** p < .001.') +
    theme_minimal(base_size=14) +
    theme(panel.grid=element_blank(),
          axis.text.x=element_text(angle=35, hjust=1, size=12, colour='black'),
          axis.text.y=element_text(face='bold',size=12, colour='black'),
          plot.title=element_text(size=18, face='bold', hjust=.5, margin=margin(b=9)),
          plot.subtitle=element_text(size=11, hjust=.5, margin=margin(b=17)),
          plot.caption=element_text(size=10, hjust=.5, margin=margin(t=13)),
          plot.margin=margin(t=24,r=45,b=28,l=35),
          legend.title=element_text(size=11), legend.position='right')
  # Wide enough for long titles; additional PDF is vector output.
  ggsave(file.path(out_dir, paste0('05_',predictor_key,'_adjusted_beta_heatmap_WIDE.png')),
         p, width=12.5, height=7, units='in', dpi=350, limitsize=FALSE, bg='white')
  ggsave(file.path(out_dir, paste0('05_',predictor_key,'_adjusted_beta_heatmap_WIDE.pdf')),
         p, width=12.5, height=7, units='in', limitsize=FALSE, bg='white')
  invisible(p)
}
make_beta_heatmap('pain','Pain associations with MDS-UPDRS measures',
                  'Age- and sex-adjusted linear regression')
make_beta_heatmap('TRPV1','TRPV1 associations with pain and MDS-UPDRS',
                  'Age- and sex-adjusted linear regression')

# =============================================================
# 06. Explain exactly why analytic n decreases
# =============================================================
# Original manuscript cohort (255 PD) not present in this input; the
# available RNA file contains 188 PD participants.
qc_waterfall <- analysis %>% filter(cohort=='PD') %>%
  group_by(EVENT_ID, months) %>%
  summarise(n_RNA_subjects=n_distinct(PATNO),
            n_TRPV1=sum(!is.na(TRPV1)),
            n_clinical_pain=sum(!is.na(NP1PAIN)),
            n_PartI=sum(!is.na(PartI_noPain)),
            n_PartII=sum(!is.na(PartII)),
            n_PartIII=sum(!is.na(PartIII)),
            n_Total=sum(!is.na(Total_noPain)),
            n_with_age_sex=sum(!is.na(age_visit) & !is.na(sex)),
            n_TRPV1_PartIII_eligible=sum(!is.na(TRPV1) & !is.na(PartIII) &
                                           !is.na(age_visit) & !is.na(sex)),
            .groups='drop')
write_csv2(qc_waterfall,'06_QC_sample_size_waterfall_by_visit.csv')
qc_missing <- analysis %>% filter(cohort=='PD') %>%
  transmute(PATNO, EVENT_ID, TRPV1_available=!is.na(TRPV1),
            pain_available=!is.na(NP1PAIN), age_available=!is.na(age_visit),
            sex_available=!is.na(sex), PartI_available=!is.na(PartI_noPain),
            PartII_available=!is.na(PartII), PartIII_available=!is.na(PartIII),
            Total_available=!is.na(Total_noPain))
write_csv2(qc_missing,'06_QC_availability_by_patient_visit.csv')
cat('Wide adjusted beta heatmaps and n-QC files saved to ',out_dir,'\n')
