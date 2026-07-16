# Pseudotime Analysis Workflow (Monocle2)

This directory contains three R scripts that together perform **pseudotime trajectory analysis** using the Monocle2 package. The workflow extracts cluster‑specific count matrices, builds a CellDataSet, infers a developmental trajectory, and produces visualisations with downstream functional enrichment.

---

## Script Overview

| Script | Purpose |
|--------|---------|
| **`extract_cluster_counts.R`** | Extracts raw count matrices for individual clusters from a merged Seurat object and saves them as separate RDS files. |
| **`run_monocle_trajectory.R`** | Builds a Monocle CellDataSet from selected clusters, estimates size factors/dispersions, identifies differentially expressed genes, and constructs a trajectory (DDRTree). |
| **`monocle_analysis_visualization.R`** | Loads a pre‑computed CDS object, generates all downstream plots (pseudotime, cell type, time‑course, density, heatmap), performs differential expression along pseudotime, and runs GO enrichment on gene clusters. |

These scripts are meant to be run **in order**, though the visualisation script can be re‑run independently after a CDS object exists.

---

## Step‑by‑Step Workflow

### 1. Extract cluster counts (`extract_cluster_counts.R`)

**Purpose**  
To create a count matrix for each cluster (e.g., `CK3D_cluster18_count.rds`) that can be used as input for Monocle.

**Input**
- A metadata file (e.g., `CK3D_metadata_res1.txt`) containing cluster assignments (`harmony_clusters_res1`).
- Cell Ranger raw count matrices for each biological replicate (10x `filtered_feature_bc_matrix` folders).

**Output**
- One RDS file per cluster, named `CK3D_clusterX_count.rds`, containing a sparse matrix of raw UMI counts.

**Example command**
```bash
Rscript extract_cluster_counts.R
```

---

### 2. Build Monocle trajectory (`run_monocle_trajectory.R`)

**Purpose**  
To select two or more clusters of interest, combine their count matrices, and infer a pseudotime trajectory.

**Input**
- The cluster count RDS files generated in the previous step.
- A metadata file (same as above) to annotate cells with their cluster ID and a new `celltype` label.

**Output**
- A Monocle CellDataSet object (`.rds` file, e.g., `CK3D_c18c20_monocle.rds`).
- A table of differentially expressed genes between the selected cell types (`.txt` file).

**Example command**
```bash
Rscript run_monocle_trajectory.R
```

---

### 3. Visualise and interpret (`monocle_analysis_visualization.R`)

**Purpose**  
To load a CDS object and generate publication‑ready figures, perform pseudotime‑dependent differential expression, and functionally annotate gene clusters.

**Input**
- A CDS object (e.g., `TT_combined_Endo_5D7D8D9D11D_monocle.rds`).
- Optional: a DEG file (can be computed on the fly).
- The `org.Brapaoleracea.eg.db` annotation package for GO enrichment.

**Output**  
All plots are saved as SVG/PNG files in the `plots/` directory. The following figures are generated:

- **`monocle_time.svg`** – trajectory coloured by pseudotime.
- **`monocle_type.svg`** – trajectory coloured by simplified cell type (e.g., Endosperm_PCD vs active).
- **`monocle_facet.svg`** – trajectory split by time point (facets).
- **`monocle_density.svg`** – density distribution of pseudotime per time point.
- **`monocle_heatmap.svg`** – heatmap of genes that vary with pseudotime, with hierarchical clustering.
- **`monocle_GO_clusterX.txt`** – GO enrichment results for each gene cluster (X = 1–5).

**Example command**
```bash
Rscript monocle_analysis_visualization.R
```

---

## Example Figures

Below is a brief description of the key output images (placeholders; actual plots will appear in the output directory).

### Trajectory coloured by Pseudotime
<img width="341" height="195" alt="image" src="https://github.com/user-attachments/assets/670f0614-7d5e-42b0-a3a1-84eb60a8a166" />

Cells are shaded from pink (early) to darkred (late) along the inferred path.

### Trajectory coloured by Cell Type
<img width="264" height="200" alt="image" src="https://github.com/user-attachments/assets/aa45cd76-9f7f-40b2-a7a0-32f208ea2ad7" />

Different cell populations (e.g., PCD vs active endosperm) are shown in distinct colours.

### Faceted Trajectory by Time Point
<img width="689" height="195" alt="image" src="https://github.com/user-attachments/assets/3d5b2037-a690-4672-98cd-e7de32187cce" />

The same trajectory is shown separately for each developmental time point (TT5D to TT11D), highlighting sample distribution along the path.

### Pseudotime Density Distribution
<img width="341" height="195" alt="image" src="https://github.com/user-attachments/assets/888a3575-1735-4b58-be67-f745336c0213" />

Overlap of pseudotime values across time points, showing when different stages are most abundant.

### Heatmap of Pseudotime‑varying Genes
<img width="746" height="237" alt="image" src="https://github.com/user-attachments/assets/4b2d296f-63d4-4ce8-9813-eabccb587134" />

Expression of significant genes ordered by pseudotime; gene clusters are identified by hierarchical clustering.

### GO Enrichment per Gene Cluster
GO terms for each gene cluster are saved as tab‑separated text files (e.g., `monocle_GO_cluster1.txt`), which can be used for functional interpretation.

---

## Dependencies

All scripts require the following R packages:

```r
install.packages(c("monocle", "dplyr", "Seurat", "ggplot2", "svglite",
                   "ggsci", "ggrepel", "scales", "tidyr", "ggpubr",
                   "clusterProfiler"))
```

The annotation package `org.Brapaoleracea.eg.db` must be built and installed (see the `setup_GO_database` workflow).

---

## Notes

- File paths in the scripts are placeholders (e.g., `path/to/your/analysis_directory`). Update them to match your own directory structure.
- The cluster count extraction script assumes a specific naming pattern for Cell Ranger output folders; adjust the `data.dir` arguments accordingly.
- If you wish to run the pipeline on a different dataset (e.g., TT instead of CK), modify the cluster and file names inside the scripts accordingly.
