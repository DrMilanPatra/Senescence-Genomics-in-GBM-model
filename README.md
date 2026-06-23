Multi-omics analysis of senescence in human glioblastoma models

Overview

This repository contains R and Python scripts used for the analysis of multi-modal genomic datasets investigating therapy-induced senescence in glioblastoma (GBM) models. The analyses integrate bulk and single-cell transcriptomics, chromatin accessibility, chromatin immunoprecipitation, spatial transcriptomics, and chromatin interaction data to characterize transcriptional and epigenetic programs associated with senescence in GBM.
Each script is designed to perform an independent analytical task and can be executed separately.

Data availability

Raw and processed datasets are not included in this repository and are available through public repositories (e.g., NCBI GEO) as described in the associated manuscript.


Data types analyzed

•	Bulk RNA-seq (proliferating vs senescent GBM models)
•	Single-cell RNA-seq (scRNA-seq GBM)
•	Single-cell ATAC-seq (scATAC-seq GBM)
•	Bulk ATAC-seq (proliferating vs senescent GBM models)
•	ChIP-seq: H3K27ac, H3k4me1, CTCF
•	Spatial transcriptomics for GBM samples
•	HiChIP chromatin interaction data for GBM samples


Software requirements


All required packages and dependencies are specified within individual scripts (R and Python).

R environment
•	R ≥ 4.0
•	Common packages include DESeq2, Seurat, Signac, edgeR, ggplot2, and Bioconductor dependencies

Python environment
•	Python ≥ 3.8
•	Common packages include scanpy, pandas, numpy, scipy, matplotlib, and related genomics tools


Analysis script description (summary):

Population level Survival analysis models:

survival_analysis_cox.r

scRNA-seq analysis of senescence: data processing, NEBULA differential analysis, cell state analysis, ligand-receptor interaction (L-R), and other analysis:

processing and doublet discrimination: scRNA_data_processing.r, run_scrublet.py
Signature high vs low cells: signature_high_vs_low_Nebula_differential_analysis.r
p16+ vs p16- cells: p16pos_vs_p16neg_Nebula_differential_analysis.r
Mutation effect included in NEBULA model: Mutation_fixed_effect_signature_high_vs_low_Nebula_differential_analysis.r
L-R interaction: run_cellphonedb_signature_statistical_3.py, run_cellphonedb_cellstate_statistical_3.py, cellphonedb_mixed models.r
cell state: Cell state assignments.r, DE_CellState_Wilcox.r

spatial transcriptomics analysis of senescence: spatial_GBM_statistics.r, 



Usage:

Rscript R_scripts/differential_expression.R
python python_scripts/scRNA_processing.py
Users must update input file paths within each script prior to execution.


Notes
•	No raw data, intermediate files, or final figures are included in this repository.
•	No unified pipeline framework is used; scripts are modular and independent.


Citation


If you use this repository, please cite the associated publication.


Contact

For questions or collaboration, please contact: [mln.patra@gmail.com]
