# Processing of Sample 306-4 and 306-3 scCITE-seq Data

Single-cell RNA sequencing data for samples **306-4** and **306-3** were processed using **Cell Ranger** (10x Genomics) and downstream **R-based workflows**: [https://github.com/cheungf/306](https://github.com/cheungf/306).

Libraries were sequenced across 24 lanes spanning multiple flowcells to enhance read depth. GEX, ADT, and VDJ libraries were combined using a unified `multi_config.csv`, allowing `cellranger multi` to merge data across lanes while preserving barcode and library integrity.

Post-processing included aggregation with `cellranger aggr` to normalize sequencing depth across lanes. Each aggregated dataset was then processed using sample-specific R scripts (`306_4_release_v01.R` or `306_3_release_v01.R`), implementing quality control, donor demultiplexing, ADT normalization, and clustering using the Seurat framework.

Gene expression and ADT data were imported using `Read10X_h5()` and assembled into Seurat objects with metadata fields for tissue, lane, donor, and flowcell. Donor assignments were determined using **DEMUXALOT2**, based on SNP profiles generated from matched bulk RNA-seq. Cells were classified as singlets, doublets, or ambiguous.

QC metrics—including mitochondrial, ribosomal, and hemoglobin gene percentages—were calculated and visualized. Low-RNA singlets were used as background for ADT denoising, while high-quality cells were retained for analysis. ADT data were normalized using **DSB** and reintegrated.

Seurat objects were normalized with **SCTransform**, merged by tissue, reduced with PCA, clustered using the Louvain algorithm, and embedded with UMAP. Cell types were annotated using reference-based classifiers.

For quality validation and metadata consistency, two custom R packages were developed to assist with barcode overlap confirmation and Cell Ranger metric summarization. Final Seurat objects were used for all downstream analyses.

**Full details, including code and configuration files, are available at** [https://github.com/cheungf/306](https://github.com/cheungf/306).
