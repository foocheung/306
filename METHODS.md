# Processing of Sample 306-4 and 306-3 scCITE-seq Data 

Single-cell RNA sequencing data for samples 306-4 and 306-3 were processed using Cell Ranger (10x Genomics) and downstream R-based analysis: [https://github.com/cheungf/306](https://github.com/cheungf/306). Samples were sequenced across 24 lanes distributed over multiple flowcells and sequencing runs to increase read depth. Each lane's data was merged across flowcells, and gene expression (GEX), antibody-derived tag (ADT), and VDJ libraries were specified in a unified `multi_config.csv` used by `cellranger multi`. This approach allowed Cell Ranger to natively combine runs from the same library prep across different lanes while preserving barcode integrity and library feature types.

Following `cellranger multi`, the output data from multiple libraries within each tissue were further aggregated using `cellranger aggr` to normalize sequencing depth across lanes. This ensured uniformity before downstream analyses.

Each `cellranger multi` and `aggr` output was processed using custom R scripts (`306_4_release_v01.R` or `306_3_release_v01.R`) implementing sample demultiplexing, ADT normalization, filtering, and clustering using the Seurat pipeline.

Gene expression and ADT data were loaded from 10X HDF5 files using `Read10X_h5()` and assembled into Seurat objects with both RNA and CITE assays. Objects were created per tissue, and metadata columns were added to capture lane identifiers, tissue origin, and donor information. Cell demultiplexing from bulk RNA-seq using `samtools mpileup` (Snakemake pipeline: [https://github.com/foocheung/Demultiplexing\_SNPSv2](https://github.com/foocheung/Demultiplexing_SNPSv2)) was used with DEMUXALOT2 to assign cells to donors. Metadata fields were appended to classify droplets as singlets (SNG), doublets (DBL), or ambiguous (AMB).

Mitochondrial, ribosomal, and hemoglobin gene percentages were calculated for QC, with violin plots generated to visualize these metrics. Cells were filtered into low RNA singlets (background for DSB normalization) and high-quality singlets (used for analysis). ADT normalization was performed using the DSB package with isotype controls and low-RNA cells to denoise protein signals. Normalized CITE matrices were reintegrated into Seurat objects.

Each object was normalized using `SCTransform` and merged by tissue. Merged datasets were reduced with PCA, clustered via the Louvain algorithm, and embedded using UMAP. Cell types were annotated using `monaco_ann1()` and `monaco_ann2()` based on reference classification. Summary CSVs and Seurat `.rds` files were saved for filtered positives, negative background cells, and unfiltered datasets.

Two in-house Shiny tools supported QC and validation:

* **barcode\_overlapper** https://github.com/foocheung/barcodeoverlapper: A Golem-based app that compares barcode overlap between libraries (e.g., GEX vs. ADT FASTQs), streaming barcodes using `ShortRead`, computing overlap matrices, and visualizing heatmaps via `ggplot2`. Used to confirm proper ADT–GEX lane alignment, detect barcode mismatches, and guide corrections to `multi_config.csv`.
* **MetricRanger** https://github.com/foocheung/MetricRanger: A Golem app that summarizes metrics from multiple `cellranger` outputs, extracts QC summaries, and previews HTML reports, enabling rapid assessment of sequencing consistency.

Due to Cell Ranger’s renumbering behavior during aggregation, original lane IDs (e.g., \[1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20]) were relabeled sequentially (e.g., \[1–12]). This relabeling did not impact results since B, CD4 T, and CD8 T cells were analyzed independently. YAML configuration files preserve original mappings.

Additional metadata fields supported donor-level harmonization. This included the `option4` subject field and a harmonized identity column where twin donor T124 was reassigned to T125. Classifications (SNG, DBL, AMB) were updated accordingly. The full pipeline was executed across three configurations (Sections A, B, and C), defined by separate YAML files.

All scripts are modular and reproducible. For assistance rerunning or customizing the workflow, please contact me.
