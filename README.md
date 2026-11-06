# TRPV1-associated inflammatory signaling and pain in Parkinson's disease

R source scripts for the PPMI longitudinal and V02 RNA-seq analyses, medication sensitivity analyses, public GEO GSE8397 brain transcriptomics, and restricted hospital PBMC correlation. This archive contains **code only**, not raw human data.

## Set input paths once (before running scripts)

Make a private copy of `00_LOCAL_PATHS.example.R` named `00_LOCAL_PATHS.R`, edit the four placeholder path values, then in RStudio run:

```r
source("00_LOCAL_PATHS.R")
source("R/01_PPMI_Match_BL_V02_Samples.R")
```

The private `00_LOCAL_PATHS.R` is ignored by Git. The original raw dataset filenames are intentionally unchanged; input directory roots are controlled through environment variables.

## Order and inclusion

### Analysis scripts by category

| Analysis category | Script numbers | Description |
|---|---|---|
| PPMI preprocessing | 01–03 | Participant matching, pain-group classification, and expression matrix generation |
| PPMI longitudinal analysis | 07–08 | TRPV1-related gene score, mixed-effects models, and longitudinal expression |
| PPMI differential expression and medication analysis | 05–06, 09–15 | DESeq2 analysis, volcano plots, heatmaps, medication sensitivity, and clinical analyses |
| GSE8397 brain transcriptomics | 16–19 | Differential expression, GO/GSEA, and enrichment visualization |
| Hospital PBMC analysis | 20 | Correlation analysis of TRP channel expression and clinical parameters |

**Script 04 was omitted** because its matrix-merging functionality was consolidated into Script 03. The original numbering was retained for traceability.

For the paired longitudinal analysis, run `01 → 02 → 03 → 07 → 08`. The 03 script is the **COMPLETE** version that successfully created 817-subject/1634-sample expression matrices, not the older strict age/sex-matching version. The longitudinal models should use the documented classified 425-subject subset and verify group counts.

For V02 DESeq2 and figures, run `05 → 06 → 09 → 10 → 11` (05 fits group-only, 06 adds age/sex, 10 is age/sex-adjusted Fig 2D, 11 is group-only medication-exclusion sensitivity for Supplementary Figure S2).

For additional clinical/medication studies run `12 → 13 → 14 → 15` as applicable; 13 requires 12's flags, and 15 requires 05 and 07's outputs. These are **not the same cohort** as the paired 425-subject series.

For public postmortem brain analysis run `16 → 17 → 18 → 19`. For hospital PBMC data run `20` after setting the private workbook path.

## Statistics and intervals

- 07 and 08 descriptive group trajectories use arithmetic mean with two-sided Student-t **95% confidence intervals**, not mean ± SE.
- 07 uses a mixed-effects model with visit × group, enrollment age, gender and participant random intercept.
- 05 uses DESeq2 group-only; 06 uses age + sex + group; 10 reads 06; 11 uses group-only refits after medication exclusions. Do **not** label 11 as age/sex adjusted.
- 14 clinical trajectories are descriptive mean with 95% confidence intervals; its other panels use age/sex-adjusted regression or LMM where specified.
- 10/11 heatmap cells show log2 fold change and BH-FDR significance, **not confidence intervals**.
- 20 uses Pearson correlation with raw p-value stars; this is **not BH corrected**, so verify alignment with Fig. 1 legend.

## Verification status

These scripts are copied from the user-supplied R analysis versions, with **only the local absolute paths replaced by environment settings and filenames cleaned**. The 03 source was taken from the working COMPLETE version. Analysis/statistical calculations are unchanged. Run and compare the output n, estimates, adjusted p-values and actual figures with the final submission before public release. Full R execution and data-based reproduction could not be performed in the current sandbox, which lacks the private PPMI/hospital datasets and the R runtime.

See `INPUT_FILE_MAP.md`, `STATISTICS_AUDIT.md`, and `FILE_NAME_MAP.md`.


