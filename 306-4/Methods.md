### Processing of Sample 306-4 scRNA-seq Data

Single-cell RNA sequencing data for sample 306-4 were processed using Cell Ranger (10x Genomics) and downstream R-based analysis. The sample was prepared using a single 10X kit but sequenced across **24 lanes distributed across multiple flowcells and sequencing runs** to increase read depth. **Each lane's data was merged across flowcells**, and gene expression, ADT, and VDJ libraries were specified in a unified `multi_config.csv` used by `cellranger multi`. This approach allowed Cell Ranger to natively combine runs from the same library prep across different lanes while preserving barcode integrity and library feature types.

Each `cellranger multi` output dataset was then processed using a custom R script, [306_4_release_v01.R](https://github.com/cheungf/CHI-306-4), which implemented sample demultiplexing, antibody-derived tag (ADT) normalization, filtering, and clustering using the Seurat pipeline.

Gene expression and ADT data were loaded from 10X HDF5 files using `Read10X_h5()` and assembled into a Seurat object with both RNA and CITE assays. **Seurat objects were created separately for each tissue type**, and additional metadata columns were added to reflect lane identifiers, tissue origin, and donor information. Cell demultiplexing assignments from DEMUXALOT2 were integrated per lane, and metadata fields were appended to classify droplets as singlets (`SNG`), doublets (`DBL`), or ambiguous (`AMB`). Mitochondrial, ribosomal, and hemoglobin gene percentages were calculated for quality control, and violin plots were generated to visualize these metrics.

Cells were filtered into two groups: low RNA content singlets (used as background for DSB normalization) and high-quality singlets (used for downstream analysis). ADT normalization was performed using the **DSB** package, leveraging isotype controls and low-RNA background cells to denoise the CITE signal. Normalized protein expression matrices were then inserted back into the Seurat object.

Seurat objects were normalized with SCTransform and merged by tissue. These merged objects were then used for downstream dimensionality reduction via PCA, clustering (using the Louvain algorithm), and UMAP embedding. Cell type annotations were assigned using `monaco_ann1()` and `monaco_ann2()` functions, which apply reference-based classification. Summary CSVs and Seurat RDS files were saved for filtered positive cells, negative background cells, and the full unfiltered dataset.

To support additional quality control and multi-modal validation, two tools were used:

1. **barcode\_overlapper** ([https://github.com/foocheung/barcode\_overlapper](https://github.com/foocheung/barcode_overlapper))
   This R script compares barcode content between GEX and ADT FASTQ files. It uses the `ShortRead` package to stream and extract barcodes, computes pairwise overlap matrices, and visualizes barcode sharing as a heatmap using `ggplot2`. The tool helps detect sample mix-ups and barcode contamination in multi-modal single-cell experiments.
   **In this project, it was instrumental in debugging a lane-level mix-up between ADT and GEX libraries across flowcells.** The barcode overlap heatmap provided direct visual evidence of mismatched or misassigned lanes, enabling rapid correction of the `multi_config.csv` file and re-processing of affected runs.

2. **MetricRanger** ([https://github.com/foocheung/MetricRanger](https://github.com/foocheung/MetricRanger))
   This Golem-based R/Shiny app summarizes multiple Cell Ranger outputs into a unified interface. It extracts run metrics from web summary files, presents them in searchable and sortable tables, and provides embedded previews of each run’s HTML report. This allows rapid assessment of sequencing quality, complexity, and consistency across runs.

Due to Cell Ranger’s internal renumbering of input lanes during aggregation, original lane identities (e.g., `[1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20]`) were relabeled sequentially (e.g., `[1–12]`). This renumbering did not affect analysis outcomes since B, CD4 T, and CD8 T cell samples were processed independently. The YAML configuration files preserve the original lane mapping.

Additional metadata fields were appended to support donor-level analysis. This included a subject assignment field (`option4`) and a harmonized identity column in which all cells from twin donor T124 were reassigned to T125. Ambiguous (`AMB`) and doublet (`DBL`) classifications were unified accordingly. The entire workflow was repeated for three configurations (Sections A, B, and C), each defined by a separate YAML file.

All scripts are modular and reproducible and can be found here : https://github.com/cheungf/CHI-306-4. For questions or assistance with rerunning or adapting the workflow, please contact [foocheung@nih.gov](mailto:foocheung@nih.gov).


