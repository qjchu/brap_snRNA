# brap_snRNA

**brap_snRNA** is a resource and toolkit for downstream analysis of single‑nucleus RNA‑sequencing (snRNA‑seq) data from **oilseed rape (*Brassica napus*)**. It provides analysis scripts, visualisation pipelines, and custom annotation databases tailored for *B. napus*, enabling a complete post‑alignment workflow from raw data loading to biological interpretation.

## Repository Structure

```
brap_snRNA/
├── raw_data_processing/              # Process Cell Ranger output into Seurat objects
│   ├── get_single_sample_rds.R       # Load a single sample
│   ├── get_multi_samples_rds.R       # Load and integrate multiple samples
│   └── readme.md                     # Detailed instructions
│
├── setup_GO_database/                # Build custom GO annotation R packages
│   ├── setup_brap_bole_GO_database.R # Build org.Brapaoleracea.eg.db from eggNOG‑mapper
│   └── readme.md                     # Step‑by‑step guide
│
└── downstream_analysis/              # Core downstream analysis modules
    ├── cell_trajectory_analysis/     # Monocle2 pseudotime trajectory inference
    ├── cell_type_annotation/         # Marker visualisation, GO enrichment, cross‑dataset correlation
    ├── cellchat_analysis/            # Cell‑cell communication with CellChat
    ├── clusters_correlation_plot/    # Network plots for cluster relationships
    └── fig_scripts/                  # Scripts for main and supplementary figures
```

---

## Modules Overview

### 1. Raw Data Processing (`raw_data_processing/`)

Scripts to convert Cell Ranger output into fully analysed Seurat objects.

| Script | Description |
|--------|-------------|
| `get_single_sample_rds.R` | Loads a single 10x sample, performs QC, normalisation, PCA, clustering (Leiden), UMAP/t‑SNE, and exports the Seurat object with metadata and marker genes. |
| `get_multi_samples_rds.R` | Loads multiple samples, merges them, applies integration methods (CCA, RPCA, Harmony), and exports the integrated object with all embeddings and cluster assignments. |

**Input**: Cell Ranger `filtered_feature_bc_matrix` directories.  
**Output**: Seurat object (`.rds`), metadata (`.txt`), cluster markers, and aggregated expression tables.

---

### 2. GO Database Setup (`setup_GO_database/`)

Builds a custom Gene Ontology annotation R package (`org.Brapaoleracea.eg.db`) for *B. napus* using eggNOG‑mapper annotations.

**Workflow**:
1. Generate `*.out.emapper.annotations` files using the [eggNOG‑mapper web server](http://eggnog5.embl.de/#/app/home).
2. Run `setup_brap_bole_GO_database.R` to parse annotations, extract GO terms, and build the package with `AnnotationForge::makeOrgPackage()`.
3. Install the package locally.

**Pre‑built packages** (Baidu Cloud) are also available for download.

---

### 3. Downstream Analysis (`downstream_analysis/`)

#### 3.1 Cell Type Annotation (`cell_type_annotation/`)

Three complementary strategies for cell‑type identification:

| Script | Description |
|--------|-------------|
| `marker_gene_anno.R` | Generates UMAP plots, sample proportion bar charts, and extensive DotPlot panels for known marker genes (PlantscRNAdb, region‑specific, subregion‑specific, qPCR, FISH, cell cycle). |
| `diff_gene_GOres_anno.R` | Performs GO enrichment (BP, CC, MF) for each cluster using `clusterProfiler` and creates bubble plots for the top terms. |
| `datasets_correlation_anno.R` | Compares cluster‑averaged expression with external LCM (GSE12404) and RNA‑seq (GSE157145) datasets from *Arabidopsis*, with and without ComBat batch correction. |

**Outputs**: UMAP SVG files, dot plots, GO enrichment tables, bubble plots, correlation heatmaps, and PCA biplots.

---

#### 3.2 Cell Trajectory Analysis (`cell_trajectory_analysis/`)

Pseudotime trajectory inference using **Monocle2**.

| Script | Description |
|--------|-------------|
| `extract_cluster_counts.R` | Extracts raw count matrices for individual clusters from a Seurat object and saves them as RDS files. |
| `run_monocle_trajectory.R` | Builds a Monocle CellDataSet from selected clusters, estimates size factors/dispersions, identifies DEGs, and constructs a DDRTree trajectory. |
| `monocle_analysis_visualization.R` | Loads a pre‑computed CDS object, generates trajectory plots (pseudotime, cell type, time‑course, density, heatmap), performs differential expression along pseudotime, and runs GO enrichment on gene clusters. |

**Outputs**: CDS objects (`.rds`), DEG tables, trajectory plots, and GO enrichment results.

---

#### 3.3 Cell‑Cell Communication (`cellchat_analysis/`)

Inference of cell‑cell communication networks using **CellChat** with a customised ligand‑receptor database for *B. napus*.

| File | Description |
|------|-------------|
| `LR_pair_atha_bchi.txt` | Ligand‑receptor pairs from *Arabidopsis* with *B. rapa* orthologs. |
| `run_cellchat_analysis.R` | Builds the custom database, runs CellChat on a sample, and exports `*_net_lr.txt`, `*_net_pathway.txt`, and a `.rds` object. |
| `plot_cellchat_summary.R` | Processes all `net_lr.txt` files and generates scatter plots showing outgoing, incoming, and self‑signals per cluster, faceted by timepoint. |
| `plot_selected_LR_pathways.R` | Creates a dot plot for a curated set of ligand‑receptor interactions, showing probability (size) and significance (colour). |

**Outputs**: Ligand‑receptor interaction tables, pathway summaries, CellChat objects, and summary scatter plots.

---

#### 3.4 Cluster Correlation Plot (`clusters_correlation_plot/`)

Network‑style plots visualising cluster relationships.

| File | Description |
|------|-------------|
| `geom-curve2.R` | Custom ggplot2 geometry for drawing curved lines with endpoint nodes. |
| `plot_cluster_links.R` | Reads node (cluster) and edge (correlation) data from Excel files and generates a network plot with nodes representing clusters (size ~ cell proportion) and curved edges representing correlations (thickness ~ strength). |

**Input**: Two Excel files – `CK_clusters.xlsx` (node positions, labels, sizes) and `CK_clusters_cor.xlsx` (edges, curvature, correlation strength).  
**Output**: SVG/PNG network plot.

---

#### 3.5 Figure Scripts (`fig_scripts/`)

Scripts that generate the main and supplementary figures for the publication.

| Script | Figure(s) | Description |
|--------|-----------|-------------|
| `plot_QC_violin.R` | Supplementary QC | Flipped violin plots showing `nFeature_RNA` and `nCount_RNA` distributions per CK sample. |
| `plot_CK_marker_dotplot.R` | Main marker dot‑plot | Faceted dot plot of conserved marker genes across CK timepoints, grouped by cell type. |
| `compare_TT_vs_CK_celltypes.R` | Figure 4A | Compares matched cell types between TT (mutant) and CK (wild‑type) at 7D; visualises genes with `|log₂FC| > 2` as a column‑scatter plot with key genes labelled. |
| `plot_hormone_gene_expression.R` | Supplementary hormone plots | For selected cell types, plots expression of ABA, IAA, GA, and XTH gene sets over time (CK vs TT) as boxplot + line plots. |

---

## Installation & Dependencies

### R Packages

All scripts require R (≥ 4.0). Install the core dependencies:

```r
install.packages(c(
  "dplyr", "Seurat", "ggplot2", "svglite", "ggsci", "ggrepel",
  "tidyr", "stringr", "patchwork", "ggpubr", "scales",
  "RColorBrewer", "reshape2", "readxl", "ggnewscale"
))
```

### Additional Packages for Specific Modules

- **CellChat**: `devtools::install_github("sqjin/CellChat")`
- **Monocle2**: `BiocManager::install("monocle")`
- **GO enrichment**: `clusterProfiler`, plus the custom `org.Brapaoleracea.eg.db` package
- **Batch correction**: `sva`
- **Heatmaps**: `pheatmap`
- **PCA visualisation**: `factoextra`

---

## Citation

If you use this repository in your research, please cite the original publication (details to be added later) and provide a link to this GitHub repository.

---

## Contact

For questions or suggestions, please open an issue on the [GitHub repository](https://github.com/qjchu/brap_snRNA).
