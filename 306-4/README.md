

## Directory Contents and Processing Overview

This directory contains all scripts and configuration files used to process the **Cell Ranger `aggr` output for sample 306-4**.

### Contents

- **Main script**: `306_4_release_v01.R`
- **Configuration files**: Three YAML files
- **Post-processing**: Small script in `T124_T125_merged/` to address item 5

### Processing Status

- `306_4_release_v01.R` has already been executed to complete items 1–4.
- Output files, including plots and Seurat objects from various levels, are available at:

```

/hpcdata/scratch/cheungf/CHI\_306\_4\_RELEASE\_v2/306\_4\_release/

```

- Cell Ranger `aggr` web summary files (`*.web_summary.html`) are also included in this directory.

### Lane Renumbering Note

Cell Ranger `aggr` does not preserve original lane numbers. Please refer to the YAML files for original lane mappings.

For example, the original B cell lanes:

```

\[1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20]

```

are renumbered as:

```

\[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

```

> This renumbering is not an issue, as B, CD4T, and CD8T samples were processed separately.

---

## Processing Plan (from internal communication)

The following strategy was implemented and is reflected in the code:

- Downsample using **Cell Ranger `aggr`** for B, CD4T, and CD8T lanes separately
- **Retain low RNA content cells** for use as background in ADT normalization
- Normalize ADT using the **DSB** package
- Create **Seurat objects with filtered cells**, processed independently for Sections A, B, and C
- Add metadata column (`option4`) to map cells to subjects
- Add a second metadata column to merge twin subjects:
  - Convert `AMB`, `DBL` of `T124` and `T125`
  - Map all `T124` entries to `T125`

---

## Questions or Support

Refining or rerunning the code should be straightforward.  
If you have any questions or need assistance, please contact **foocheung@nih.gov**.


