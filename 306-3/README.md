
## Directory Contents and Processing Scripts 306-3 Aggr

- All scripts for processing the Cell Ranger aggr output for sample 306-3 are located in this directory.
- One script, `306_3_release_v01.R`, and five YAML files were used.
- The script `306_3_release_v01.R` has already been used to process the data

- The output of `306_3_release_v01.R`, including plots and Seurat outputs from various levels, is located at:  
  `/hpcdata/scratch/cheungf/306_3_RELEASE/`
- The Cell Ranger aggr web summary files (ending with `web_summary.html`) are located in this directory.


**Note:** 
All B, CD4T, and CD8T lanes are processed separately. Please refer to the YAML files for verification.
Note only **1 lane** for **CD4T and CD8T** and **did not** get processed with cell ranger aggr (check yaml files config_3.yaml and config_2.yaml)
Only **B lanes** got processed with cell ranger aggr. There are **3** different versions B, again check the yaml files as below depending which lanes you want to include/exclude:

**config_1.yaml:**
<PRE>
data_path: "/hpcdata/chi/PROJECTS_Archive/2022_CHI_PROPOSALS/Manthiram_Covid-tonsil_CHI-306/NEW_DATA_9_sept_2023/230908_A00941_1368_AHKVGVDSX2/RUN_23_306_3/RUN/NEW2_3AGGR/outs/count/filtered_feature_bc_matrix.h5"
lanes: [1, 2, 3, 4, 5, 6]
output_prefix: "306_3_release/B"
</PRE>
**config_4.yaml:**
<PRE>
data_path: "/hpcdata/chi/PROJECTS_Archive/2022_CHI_PROPOSALS/Manthiram_Covid-tonsil_CHI-306/NEW_DATA_9_sept_2023/230908_A00941_1368_AHKVGVDSX2/RUN_23_306_3/RUN/NEW2_1AGGR/outs/count/filtered_feature_bc_matrix.h5"
lanes: [1, 2, 3, 6]
output_prefix: "306_3_release/B_1"
</PRE>
**config_5.yaml:**
<PRE>
data_path: "/hpcdata/chi/PROJECTS_Archive/2022_CHI_PROPOSALS/Manthiram_Covid-tonsil_CHI-306/NEW_DATA_9_sept_2023/230908_A00941_1368_AHKVGVDSX2/RUN_23_306_3/RUN/NEW2_2AGGR/outs/count/filtered_feature_bc_matrix.h5"
lanes: [1, 2, 3, 4, 6]
output_prefix: "306_3_release/B_2"
</PRE>


From useful email, this has been addressed as above:

<PRE>
"Downsampling with cellranger aggregate for B, CD4T, CD8T lanes separately.
Keep low RNA content cells for using as background for ADT normalization.
Perform ADT normalization with DSB package and create Seurat object with filtered cells.
    1.Conduct this for section A, B, C separately
Use demuxalot to demultiplex all lanes.
Add demuxalot results as metadata for mapping cells to subjects."
</PRE>


I also ran the data using Seurat5 and using the "Integrative analysis in Seurat v5". I can upload the data (rds) files, code and plots before and after intergration for 306-4 for CD4,CD8 and B cells and 306-3 for B cells across lanes within each run. I will post next week after our bioinfo meeting. From what I can see very small batch to batch effects across lanes within the same run within same tissue

**Refining or rerunning the code, should be straightforward. Questions or help, please send an email to me.**
