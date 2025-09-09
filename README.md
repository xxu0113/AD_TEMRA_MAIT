# AD_TEMRA_MAIT

I read tony wyss's paper named Clonally expanded CD8 T cells patrol the cerebrospinal fluid in Alzheimer’s disease.

They think there is increased CD8⁺ TEMRA (T effector memory CD45RA⁺) cells in AD CSF and it correlated with bad cognitive score.

I am wondering if there is increased MAIT signature in AD their dataset for either blood or CSF. (They have human scRNA-seq from 18 CSF samples and 13 blood samples from AD, MCI, and HC people.)

I think there may be increased MAIT cell signature (cluster 5 for the blood) in AD, but I didn't see MAIT cell signature in CSF except KLRB1 expression in a lot of cells. This led me wonder maybe in AD, the effect of MAIT cells may be more through periphery? 

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
- `outputs/` large rds ignored
- `data/` Raw GEO downloads (GSE134576 / GSE134577) ignored.


## Results for blood analysis
I think in blood, in AD, there are 2 sample/patients with higher mait cell in terms of count.
<img width="1142" height="1278" alt="image" src="https://github.com/user-attachments/assets/6d8e704d-3361-438d-84c4-c2f11743cf58" />

I think in blood, in AD, there are 2 sample/patients with higher mait cell in terms of MAIT propotion over all blood cells.
<img width="1142" height="1278" alt="image" src="https://github.com/user-attachments/assets/6a28eb93-316f-4cd4-88dd-6b66bd89f44d" />




