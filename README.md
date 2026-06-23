Multi-omics analysis of senescence in glioblastoma models

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


Repository structure

•	R_scripts/ : R-based analysis scripts for data processing, differential expression, chromatin analysis, and statistical models.
•	python_scripts/ : Python-based computational analysis for scrublet, CellPhonedb
•	config/ : Optional configuration files


Each script is fully standalone and can be executed independently.


Software requirements


All required packages and dependencies are specified within individual scripts (R and Python).

R environment
•	R ≥ 4.0
•	Common packages include DESeq2, Seurat, Signac, edgeR, ggplot2, and Bioconductor dependencies

Python environment
•	Python ≥ 3.8
•	Common packages include scanpy, pandas, numpy, scipy, matplotlib, and related genomics tools


Analysis workflow (summary)


Analyses include:
•	Quality control and preprocessing of sequencing data
•	Transcriptomic profiling of senescent and proliferating GBM models
•	Chromatin accessibility analysis (ATAC-seq and scATAC-seq)
•	Chromatin state analysis (ChIP-seq)
•	Chromatin interaction analysis (HiChIP)
•	Integration of multi-omics datasets to characterize senescence-associated regulatory programs in GBM models


Usage


Each script is intended to be executed independently.
Example:
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
