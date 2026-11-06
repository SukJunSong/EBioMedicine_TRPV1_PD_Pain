# Raw input files and paths (input names preserved)

Raw dataset **basenames have NOT been changed** to avoid breaking the PPMI/GEO matching and joins. Only local machine root paths have been replaced with environment variables. The repository contains no raw datasets.

Set in a private `00_LOCAL_PATHS.R`:

- `PD_TRPV1_ROOT`: directory that previously contained `PPMI/` and `V02_DESeq2_Original_vs_Medication/`, `V02_DESeq2_AgeSex_Adjusted/`.
- `PD_TRPV1_CLINICAL_DIR`: clinical CSV directory containing `Age_at_visit_20Mar2025.csv` and MDS-UPDRS CSV files.
- `PD_TRPV1_GSE8397_DIR`: output/input directory for GSE8397 GSEA and the optional series matrix.
- `HOSPITAL_CORRELATION_XLSX`: restricted hospital correlation workbook.

Examples of **exact input basenames** required by original code:

| Input | Description | Scripts |
|---|---|---|
| `PPMI/genes BL/*.genes.sf` | Salmon BL gene quantification | 01, 03 |
| `PPMI/genes V02/*.genes.sf` | Salmon V02 gene quantification | 01, 03, 05 |
| `PPMI/metaData V02.xlsx` | Participant V02 metadata | 02, 03, 05, 06 |
| `PPMI/painData V02.xlsx` | MDS-UPDRS item 1.9 | 02, 05 |
| `Participant_Status_19May2026.csv` | Enrollment age (searches PPMI folder, clinical folder, project root) | 03 |
| `Age_at_visit_20Mar2025.csv` | Age at visit | 06, 13, 14 |
| `Concomitant_Medication_Log_26Oct2026.csv` | Medication | 05, 12 |
| `LEDD_Concomitant_Medication_Log_26Oct2026.csv` | Medication | 05, 12 |
| `MDS-UPDRS_Part_I_Patient_Questionnaire_25Feb2025.csv` | Visit pain assessment | 05, 12, 14 |
| `PPMI/hsa04750pd_cleaned.xlsx` and `PPMI/hsa04750ctrl_cleaned.xlsx` | Workbook-based analysis separate from 425-person paired analysis | 12, 13, 14 |
| `GSE8397-GPL96_series_matrix.txt` | GEO microarray expression matrix | 16 |
| `Hospital correlation.xlsx` (actual path via env variable) | Non-public PBMC dataset | 20 |

**Important:** PPMI and hospital participant-level data should NOT be committed to a public repository. Even generated result CSVs with patient identifiers should be private.
