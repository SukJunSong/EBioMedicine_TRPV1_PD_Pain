# Audit of statistical methods and confidence intervals

This audit was conducted by reading the uploaded R source code (static inspection). Dataset-level reproduction has NOT been checked. No existing statistical test, coefficient, summary, or significance threshold was changed in the published-code copies.

| Script | Observed method | CI treatment |
|---|---|---|
| 01-03 | PPMI file matching, pain grouping and matrices | Not applicable |
| 05 | DESeq2 group-only model; outputs both default DESeq2 `padj` and manual BH adjusted p | Not applicable |
| 06 | DESeq2 `~ age + sex + group` | Not applicable |
| 07 | Gene-wise log2(TPM+1) z-score composite, visit/group LMM with age and gender and subject intercept | Pointwise Student-t 95% CI around raw group/visit means |
| 08 | Individual gene longitudinal log2(TPM+1) trajectories | Pointwise Student-t 95% CI around raw group/visit means |
| 09 | Volcano from 06 age-sex adjusted results | Not applicable |
| 10 | 06-adjusted log2FC/BH-FDR heatmap | Not applicable |
| 11 | Original 05 DESeq2 group-only and group-only refits for medication sensitivity | Not applicable |
| 12 | Analgesic/dopaminergic exposure categorization | Not applicable |
| 13 | Visit-specific HC vs PD linear models and medication sensitivity | Inspect output n and FDR family |
| 14 | Age-/sex-adjusted clinical regressions and mixed models; raw descriptive trajectories | Descriptive mean ± Student-t 95% CI |
| 15 | Original DESeq2 direct contrast and separate age/sex-adjusted score regression | Regression confidence interval (95%) |
| 16 | GSE8397 limma PD-vs-control microarray | Not applicable |
| 17 | Region-specific GO fgsea | Not applicable; BH-FDR `padj` |
| 18-19 | GO descriptions, GSEA NES dotplots and `padj < 0.05` marker | Not applicable |
| 20 | Pearson pairwise-complete correlation | No 95% CI; raw p significance stars |

**Unresolved checks before public release**:

1. Confirm that the final manuscript Fig. 1 correlation legend uses **raw p** rather than FDR-adjusted p, or update one of them as a separate scientific correction. The code has NOT been altered.
2. Confirm the difference between paired BL/V02 425 classified samples and cross-sectional V02 DESeq2 451 subjects is appropriately described.
3. Confirm GSE8397 panel annotation/plot matches final Fig. 3 and any GSEA statistical claims.
4. Medication heatmap S2 (11) refits a **group-only** model; it is not adjusted for age/sex. Keep Methods/caption consistent.
5. R package version and environment files (`sessionInfo()`) should be recorded after real execution.
