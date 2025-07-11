# Processing of Sample 306-4 and 306-3 scCITE-seq Data

Single-cell RNA sequencing data for samples **306-4** and **306-3** were processed using **Cell Ranger** (10x Genomics) and downstream **R-based workflows**. Each sample was prepared using a **single 10X Chromium kit**, but sequencing was performed across **24 independent lanes**, spanning multiple **flowcells** and **NovaSeq runs**, to maximize read depth and library complexity.

To ensure complete data integration across these distributed runs:

* **All flowcell FASTQs for a given lane were combined**, preserving lane identity.
* A single `multi_config.csv` was used per sample with `cellranger multi` to define all **GEX**, **ADT**, and **VDJ** libraries in a unified run.
* After `cellranger multi`, a **lane-level normalization step** was performed using `cellranger aggr`, which scaled libraries by mapped read depth to mitigate technical variation across lanes.

This workflow allowed Cell Ranger to **accurately integrate gene expression and protein tag data from the same physical library prep**, while maintaining barcode integrity and enabling cross-lane aggregation.

---

Each `cellranger multi` output was then passed into a sample-specific custom R script—either [306_4_release_v01.R](https://github.niaid.nih.gov/cheungf/306/blob/main/306-4/306_4_release_v01.R) or [306_3_release_v01.R](https://github.niaid.nih.gov/cheungf/306/blob/main/306-3/306_3_release_v01.R) — which executed the following major steps:

1. **Loading and Assembly of Seurat Objects**
   Using `Read10X_h5()`, gene expression and antibody-derived tag (ADT) matrices were imported from the HDF5 files into Seurat. Both `"RNA"` and `"CITE"` assays were initialized.

   * **Seurat objects were created separately for each tissue type** and labeled with metadata fields such as:

     * `Lane` (inferred from barcode)
     * Tissue origin
     * Flowcell and donor information

2. **Demultiplexing via Bulk RNA-Seq SNP Matching**
   Cell assignments were determined using [DEMUXALOT2](https://github.com/foocheung/Demultiplexing_SNPSv2), a tool that compares single-cell expression barcodes to bulk SNP profiles derived using `samtools mpileup`. For each lane:

   * Refined donor assignments were read from `assignments_refined.tsv`.
   * Cells were labeled as singlets (`SNG`), doublets (`DBL`), or ambiguous (`AMB`).
   * Barcode suffixes were appended to ensure per-lane uniqueness and prevent collisions during merging.

3. **Quality Control**
   For each cell:

   * **Mitochondrial gene content** (`percent.mito`) was calculated using `^MT-` feature patterns.
   * **Ribosomal and hemoglobin content** (`percent_ribo`, `percent_hb`) was computed using regex-matched gene symbols.
   * Pre-filtering QC plots were generated via `VlnPlot()` and saved as `*_pre-filtering.pdf`.

4. **Filtering and ADT Normalization via DSB**
   Cells were separated into:

   * **Low-RNA singlets** (`<500 features`) as **background** for ADT denoising
   * **High-quality singlets** (filtered by feature count, RNA count, and mito%) for downstream analysis
     ADT matrices were normalized using [**DSB**](https://cran.r-project.org/web/packages/dsb/) with isotype controls and inserted back into the Seurat object.

5. **Lane-Specific Normalization and Merging**

   * Seurat objects were split by `Lane` and normalized using `SCTransform()`.
   * After normalization, they were merged, and PCA was performed.
   * Dimensionality reduction and clustering were done via:

     * `RunPCA()`
     * `FindNeighbors()` and `FindClusters()` (Louvain)
     * `RunUMAP()` using RNA-based PCA coordinates

6. **Annotation and Output Generation**

   * Cell types were assigned using `monaco_ann1()` and `monaco_ann2()`, which map cells to reference immune types.
   * Final Seurat objects were saved as:

     * `*_pos_monaco.rds`: High-quality, filtered, annotated
     * `*_pos.rds`: High-quality, filtered
     * `*_neg.rds`: Background (low RNA)
     * `*.rds`: Full object before filtering
   * UMAPs split by lane were saved using `DimPlot_scCustom()` as `*_UMAP.pdf`.

---

### Quality Control Tools

To further support validation and reproducibility, two in-house R/Shiny tools were developed:

1. **[barcode\_overlapper](https://github.com/foocheung/barcodeoverlapper)**

   * Compares barcode overlap between FASTQ-derived libraries (e.g., GEX vs. ADT)
   * Streams barcodes using `ShortRead`, computes overlap matrices, and visualizes sharing via `ggplot2` heatmaps
   * **Used here to confirm per-lane alignment between GEX and ADT libraries** across flowcells
   * Helped identify and correct sample mismatches via updates to `multi_config.csv`

2. **[MetricRanger](https://github.com/foocheung/MetricRanger)**

   * Summarizes `cellranger` output across all runs
   * Aggregates HTML metrics into a searchable table and provides visual previews
   * Enables easy tracking of sequencing depth, library complexity, and sample outliers

---

### Additional Notes

* During `cellranger aggr`, **original lane numbers** (e.g., `[1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20]`) were automatically **renumbered sequentially** (e.g., `[1–12]`) by Cell Ranger. This did **not impact downstream results**, as tissue types were processed independently.
* YAML configuration files explicitly retain original lane assignments.
* A harmonized subject ID column (`option4`) was added for donor-level analyses. Twin donor `T124` was reassigned to `T125`, and cell classification fields (e.g., `AMB`, `DBL`) were unified.

This workflow was executed three times for different processing groups—**Sections A, B, and C**—each defined by a separate YAML configuration file.

### Reproducibility and Transparency
This pipeline is shared to support reproducibility of the dataset and pre-analysis steps. By providing the exact code, configurations, and processing logic, other researchers can:

* Reproduce the processed Seurat objects
* Audit technical decisions
* Apply or adapt the workflow to new datasets

This aligns with FAIR data principles and promotes long-term reusability and transparency.
