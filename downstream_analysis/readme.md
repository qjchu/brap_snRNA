# Downstream Analysis

This directory contains a comprehensive collection of R scripts for **downstream analysis of snRNA‑seq data** from *Brassica napus*. The workflows cover the major steps after initial clustering and integration, including cell‑type annotation, pseudotime trajectory inference, cell‑cell communication, and figure generation.

---

## Directory Structure

```
downstream_analysis/
├── cell_trajectory_analysis/      # Monocle2 pseudotime trajectory analysis
├── cell_type_annotation/          # Marker visualisation, GO enrichment, cross‑dataset correlation
├── cellchat_analysis/             # Cell‑cell communication inference with CellChat
├── clusters_correlation_plot/     # Network plots for cluster relationships
└── fig_scripts/                   # Scripts for main and supplementary figures
```

---

## Module Overview

| Module | Description | Key Scripts |
|--------|-------------|-------------|
| **cell_trajectory_analysis** | Pseudotime trajectory inference using Monocle2, from count extraction to visualisation and GO enrichment of trajectory‑associated genes. | `extract_cluster_counts.R`, `run_monocle_trajectory.R`, `monocle_analysis_visualization.R` |
| **cell_type_annotation** | Three complementary strategies for cell‑type identification: marker gene dot plots, GO enrichment, and cross‑dataset correlation with external *Arabidopsis* data. | `marker_gene_anno.R`, `diff_gene_GOres_anno.R`, `datasets_correlation_anno.R` |
| **cellchat_analysis** | Custom ligand‑receptor database for *B. napus* and CellChat analysis to infer cell‑cell communication networks across timepoints. | `run_cellchat_analysis.R`, `plot_cellchat_summary.R`, `plot_selected_LR_pathways.R` |
| **clusters_correlation_plot** | Network‑style plots visualising cluster relationships, with nodes representing clusters and edges representing correlation strengths. | `plot_cluster_links.R` (requires `geom-curve2.R`) |
| **fig_scripts** | Scripts that generate the main and supplementary figures for the publication, covering QC, marker dot plots, differential expression, and hormone gene expression. | `plot_QC_violin.R`, `plot_CK_marker_dotplot.R`, `compare_TT_vs_CK_celltypes.R`, `plot_hormone_gene_expression.R` |

---

## Module Details

### 1. Cell Type Annotation (`cell_type_annotation/`)

**Purpose**: Assign biological identities to clusters using three complementary strategies.

**Scripts**:
- **`marker_gene_anno.R`** – Generates UMAP plots, sample proportion bar charts, and extensive DotPlot panels for known marker genes (PlantscRNAdb, region‑specific, subregion‑specific, qPCR, FISH, cell cycle).
- **`diff_gene_GOres_anno.R`** – Performs GO enrichment (BP, CC, MF) for each cluster using `clusterProfiler` and creates bubble plots for the top terms.
- **`datasets_correlation_anno.R`** – Compares cluster‑averaged expression with external LCM (GSE12404) and RNA‑seq (GSE157145) datasets from *Arabidopsis*, with and without ComBat batch correction.

**Outputs**: UMAP SVG files, dot plots, GO enrichment tables (`*_cluster_marker_GOres.txt`), bubble plots, correlation heatmaps, dendrograms, and PCA biplots.

---

### 2. Cell Trajectory Analysis (`cell_trajectory_analysis/`)

**Purpose**: Infer developmental trajectories using Monocle2.

**Workflow**:
1. **`extract_cluster_counts.R`** – Extracts raw count matrices for individual clusters from a merged Seurat object and saves them as RDS files.
2. **`run_monocle_trajectory.R`** – Builds a Monocle CellDataSet from selected clusters, estimates size factors/dispersions, identifies differentially expressed genes, and constructs a DDRTree trajectory.
3. **`monocle_analysis_visualization.R`** – Loads a pre‑computed CDS object and generates downstream plots (pseudotime, cell type, time‑course, density, heatmap), performs differential expression along pseudotime, and runs GO enrichment on gene clusters.

**Outputs**: Monocle CDS objects (`.rds`), DEG tables, trajectory plots (SVG/PNG), and GO enrichment results.

---

### 3. Cell‑Cell Communication (`cellchat_analysis/`)

**Purpose**: Infer cell‑cell communication networks using CellChat with a customised ligand‑receptor database for *B. napus*.

**Workflow**:
1. **Prepare the database** – `LR_pair_atha_bchi.txt` is formatted into a CellChat‑compatible list.
2. **Run CellChat** – `run_cellchat_analysis.R` loads the Seurat object, assigns the custom database, identifies overexpressed genes/interactions, and computes communication probabilities. Exports `*_net_lr.txt`, `*_net_pathway.txt`, and a `.rds` object.
3. **Summarise** – `plot_cellchat_summary.R` processes all `net_lr.txt` files and generates scatter plots showing outgoing, incoming, and self‑signals per cluster, faceted by timepoint.
4. **Visualise selected pathways** – `plot_selected_LR_pathways.R` creates a dot plot for a curated set of ligand‑receptor interactions, showing probability (size) and significance (colour).

**Outputs**: Ligand‑receptor interaction tables, pathway summaries, CellChat objects (`.rds`), and summary scatter plots.

---

### 4. Cluster Correlation Plot (`clusters_correlation_plot/`)

**Purpose**: Visualise relationships between cell clusters as a network plot.

**Input**:
- **Node file** (`CK_clusters.xlsx`) – defines cluster positions (`newX`, `newY`), tissue group (`group`), cell type label (`celltype3`), and node size (`value3`).
- **Edge file** (`CK_clusters_cor.xlsx`) – defines connections between clusters (`startx`, `starty`, `endx`, `endy`), curvature, and correlation strength (`size2`).

**Script**: `plot_cluster_links.R` uses the custom `geom_curve2` geometry to draw curved edges with endpoints, producing a network plot where:
- Node size reflects cell proportion.
- Edge thickness reflects correlation strength.
- Colour indicates tissue identity (Endosperm, EP, SC).

**Output**: SVG/PNG network plot.

---

### 5. Figure Scripts (`fig_scripts/`)

**Purpose**: Generate the main and supplementary figures for the publication.

| Script | Figure(s) | Description |
|--------|-----------|-------------|
| `plot_QC_violin.R` | Supplementary QC | Flipped violin plots showing `nFeature_RNA` and `nCount_RNA` distributions per CK sample. |
| `plot_CK_marker_dotplot.R` | Main marker dot‑plot | Faceted dot plot of conserved marker genes across CK timepoints, grouped by cell type. |
| `compare_TT_vs_CK_celltypes.R` | Figure 4A | Compares matched cell types between TT (mutant) and CK (wild‑type) at 7D; visualises genes with `|log₂FC| > 2` as a column‑scatter plot with key genes labelled. |
| `plot_hormone_gene_expression.R` | Supplementary hormone plots | For selected cell types, plots expression of ABA, IAA, GA, and XTH gene sets over time (CK vs TT) as boxplot + line plots. |

**Outputs**: SVG/PNG figures saved to the appropriate output directories.

---

## Dependencies

All scripts require R (≥ 4.0) and the following packages (install as needed):

```r
install.packages(c(
  "dplyr", "Seurat", "ggplot2", "svglite", "ggsci", "ggrepel", 
  "tidyr", "stringr", "patchwork", "ggpubr", "scales", 
  "RColorBrewer", "reshape2", "readxl", "ggnewscale"
))
```

Additional packages for specific modules:

- **CellChat**: `CellChat` (available from GitHub: `devtools::install_github("sqjin/CellChat")`)
- **Monocle2**: `monocle` (Bioconductor: `BiocManager::install("monocle")`)
- **GO enrichment**: `clusterProfiler`, `org.Brapaoleracea.eg.db` (custom package)
- **Batch correction**: `sva`
- **Heatmaps**: `pheatmap`
- **PCA visualisation**: `factoextra`
