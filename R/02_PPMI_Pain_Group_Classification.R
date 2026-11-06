# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 02_PPMI_Pain_Group_Classification.R
# Requires original project datasets; full R execution not validated here.
# =========================
# 0. Package
# =========================

library(readxl)
library(dplyr)
library(readr)
library(stringr)


# =========================
# 1. 경로 설정
# =========================

base_dir <- file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), "PPMI")

common_file <- file.path(base_dir, "BL_V02_common_PATNO_file_list.csv")
metadata_file <- file.path(base_dir, "metaData V02.xlsx")
paindata_file <- file.path(base_dir, "painData V02.xlsx")


# =========================
# 2. 파일 불러오기
# =========================

common_df <- read_csv(common_file, show_col_types = FALSE)

meta_df <- read_excel(metadata_file)

pain_df <- read_excel(paindata_file)


# =========================
# 3. column 이름 확인
# =========================
# 한 번 확인용
# print(colnames(common_df))
# print(colnames(meta_df))
# print(colnames(pain_df))


# =========================
# 4. PATNO 형식 통일
# =========================

common_df <- common_df %>%
  mutate(PATNO = as.character(PATNO))

meta_df <- meta_df %>%
  mutate(PATNO = as.character(PATNO))

pain_df <- pain_df %>%
  mutate(PATNO = as.character(PATNO))


meta_group <- meta_df %>%
  select(PATNO, DIAGNOSIS) %>%
  distinct() %>%
  mutate(
    DIAGNOSIS = str_trim(as.character(DIAGNOSIS)),
    Diagnosis_Group = case_when(
      DIAGNOSIS == "PD" ~ "PD",
      DIAGNOSIS == "Control" ~ "Control",
      TRUE ~ NA_character_
    )
  )
# =========================
# 6. painData에서 NP1PAIN 정보 정리
# =========================

pain_group <- pain_df %>%
  select(PATNO, NP1PAIN) %>%
  distinct() %>%
  mutate(
    NP1PAIN = as.numeric(NP1PAIN),
    Pain_Group = case_when(
      NP1PAIN == 0 ~ "No Pain PD",
      NP1PAIN == 1 ~ "Low Pain PD",
      NP1PAIN >= 2 ~ "High Pain PD",
      TRUE ~ NA_character_
    )
  )


# Prevent many-to-many joins silently multiplying subjects.
for (label in c("meta_group", "pain_group")) {
  z <- get(label)
  if (anyDuplicated(z$PATNO)) {
    stop(label, " contains multiple records per PATNO. Resolve conflicting records before continuing.")
  }
}

# =========================
# 7. 공통 PATNO 파일 리스트에 meta / pain 정보 붙이기
# =========================

final_df <- common_df %>%
  left_join(meta_group, by = "PATNO") %>%
  left_join(pain_group, by = "PATNO") %>%
  mutate(
    Final_Group = case_when(
      Diagnosis_Group == "Control" ~ "Control",
      Diagnosis_Group == "PD" & Pain_Group == "No Pain PD" ~ "No Pain PD",
      Diagnosis_Group == "PD" & Pain_Group == "Low Pain PD" ~ "Low Pain PD",
      Diagnosis_Group == "PD" & Pain_Group == "High Pain PD" ~ "High Pain PD",
      Diagnosis_Group == "PD" & is.na(Pain_Group) ~ "PD_Pain_NA",
      TRUE ~ "Unknown"
    )
  )


# =========================
# 8. 그룹별 숫자 확인
# =========================

print(table(final_df$Diagnosis_Group, useNA = "ifany"))
print(table(final_df$Final_Group, useNA = "ifany"))


# =========================
# 9. 최종 파일 저장
# =========================

write_csv(
  final_df,
  file.path(base_dir, "BL_V02_common_PATNO_with_Diagnosis_and_PainGroup.csv")
)


# =========================
# 10. 그룹별로 따로 저장
# =========================

control_df <- final_df %>%
  filter(Final_Group == "Control")

no_pain_pd_df <- final_df %>%
  filter(Final_Group == "No Pain PD")

low_pain_pd_df <- final_df %>%
  filter(Final_Group == "Low Pain PD")

high_pain_pd_df <- final_df %>%
  filter(Final_Group == "High Pain PD")

write_csv(control_df, file.path(base_dir, "BL_V02_Control_common_PATNO.csv"))
write_csv(no_pain_pd_df, file.path(base_dir, "BL_V02_No_Pain_PD_common_PATNO.csv"))
write_csv(low_pain_pd_df, file.path(base_dir, "BL_V02_Low_Pain_PD_common_PATNO.csv"))
write_csv(high_pain_pd_df, file.path(base_dir, "BL_V02_High_Pain_PD_common_PATNO.csv"))

cat("완료: BL-V02 공통 PATNO 기준으로 Control / No Pain PD / Low Pain PD / High Pain PD 분류 저장했습니다.\n")
