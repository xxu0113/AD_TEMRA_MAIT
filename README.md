# AD_TEMRA_MAIT

I read tony wyss's paper named Clonally expanded CD8 T cells patrol the cerebrospinal fluid in Alzheimer’s disease.

They think there is increased CD8⁺ TEMRA (T effector memory CD45RA⁺) cells in AD CSF and it correlated with bad cognitive score.

I am wondering if there is increased MAIT signature in their dataset (They have scRNA-seq from 18 CSF samples and 13 blood samples from AD, MCI, and HC people.)


## Source
- Paper: https://www.nature.com/articles/s41586-019-1895-7#Sec2
- CSF dataset: [GSE134577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134577)
- Blood dataset: [GSE134576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134576)

## Structure
- `scripts/`
  - `01_temra_seurat_MAITscore.R` — CSF scRNA-seq analysis with QC, clustering, and MAIT module scores.
  - `02_blood_seurat_MAITmarkers.R` — Blood scRNA-seq analysis with QC, clustering, and MAIT marker expression.
- `reports/`
  - `CSF_Blood_MAIT.qmd` — Quarto report integrating CSF + blood results.
  - `CSF_Blood_MAIT.html/pdf` — Rendered reports integrating CSF + blood results.
- `outputs/` large rds *(ignored in Git)*, but keep figures
- `data/` *(ignored in Git)* Raw GEO downloads (GSE134576 / GSE134577).


