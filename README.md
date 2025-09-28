
# Airway RNA-seq Case Study — DESeq2 + Random Forest (ICR-ready)

![Run pipeline](https://github.com/<your-username>/icr-airway-case-study/actions/workflows/run-rnaseq-pipeline.yml/badge.svg)

**Goal:** End-to-end bulk RNA-seq analysis on the Bioconductor `airway` dataset with reproducible code, clear figures, and a tiny ML demo — an artefact to link in CVs (ICR: RNA-seq, R/Python/Bash/Git, HPC-friendly).

## Data
- Bioconductor `airway` dataset (human airway smooth muscle; dexamethasone treated vs untreated).

## Methods (high-level)
1) Differential expression with **DESeq2** (QC/filter → DE table).
2) Visualisations: **Volcano** and **PCA** (VST).
3) Toy ML: **Random Forest** predicts sample condition from top-variable VST genes.
   > Note: sample size is small; this is a demonstration of workflow competence.

## Outputs
- `figures/volcano.png`, `figures/pca.png`, `figures/roc.png`
- `results/DESeq2_results.csv`, `results/metrics.txt`, `results/sessionInfo.txt`

## Reproduce locally
```bash
R -q -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
R -q -e "BiocManager::install(c('DESeq2','airway'), ask=FALSE, update=TRUE)"
R -q -e "install.packages(c('ggplot2','randomForest','pROC','matrixStats'), repos='https://cloud.r-project.org')"
Rscript scripts/airway_rnaseq_pipeline.R
```

## Run in CI (GitHub Actions)
A workflow is included at `.github/workflows/run-rnaseq-pipeline.yml` using the Bioconductor Docker image for reliable installs. Enable it in the **Actions** tab.

## Key Figures
*(appear after running the pipeline)*
![Volcano](figures/volcano.png)
![PCA](figures/pca.png)
![ROC](figures/roc.png)

## Case Study (1-page)
➡️ See [CASE_STUDY.md](CASE_STUDY.md)

## Citation
Please cite `airway` and `DESeq2`. See `CITATION.cff`.

## Troubleshooting
- **Error: there is no package called 'DESeq2'**  
  Install Bioconductor packages first:
  `R -q -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('DESeq2','airway'), ask=FALSE, update=TRUE)"`
- On Windows, install **Rtools** matching your R version if compilation is needed.
