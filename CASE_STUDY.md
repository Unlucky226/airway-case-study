
# Bulk RNA-seq Case Study — DESeq2 + Random Forest (Airway)

**Objective.** Demonstrate an end-to-end bulk RNA-seq pipeline with reproducible code, clear figures, and a small ML demo suitable for fast-moving translational teams (ICR fit).

**Data.** Bioconductor `airway` (human airway smooth muscle; dexamethasone treated vs untreated).

**Methods.**
- Differential expression: `DESeq2` (QC/filter → DE table)
- Visualisations: Volcano (log2FC vs -log10 padj) and PCA (VST)
- Toy classifier: Random Forest predicts sample condition from top-variable VST genes (demonstration only; tiny n)

**Results.**
- PCA shows treatment-driven separation.
- DE table saved to `results/DESeq2_results.csv`.
- RF demo reports accuracy/AUC (see `results/metrics.txt`) and `figures/roc.png` if both classes present in test split.

**Reproducibility.**
- Single script: `scripts/airway_rnaseq_pipeline.R`
- CI workflow (GitHub Actions) runs the pipeline in a clean environment
- `sessionInfo` logged to `results/sessionInfo.txt`

**Why it matches the ICR spec.**
- RNA-seq (bulk), R/Bioconductor, Linux/HPC-ready, Git, CI/reproducibility, clear figures for team discussion.
