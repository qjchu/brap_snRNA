# Cell Type Annotation

This directory contains scripts for annotating cell types in our *Brassica napus* snRNA‑seq dataset. We adopt **three complementary strategies** to assign biological identities to clusters:

1. **Marker gene visualisation** – dot plots of known marker genes from published resources (PlantscRNAdb, region‑specific, subregion‑specific, qPCR, FISH, cell cycle).
2. **GO enrichment analysis** – functional enrichment of cluster‑specific marker genes to infer biological processes and cellular compartments.
3. **Cross‑dataset correlation** – compare our cluster‑averaged expression profiles with external LCM and RNA‑seq datasets from *Arabidopsis* to validate conservation of cell types.

Together, these approaches provide robust evidence for cell type assignment.

---

## Scripts Overview

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `marker_gene_anno.R` | Generates UMAP plots, sample proportion bar charts, and a series of DotPlot panels with known marker genes. | `*_metadata.txt`<br>`*_merge.rds` | UMAP SVG files, bar plot, dot plots (viewed interactively) |
| `diff_gene_GOres_anno.R` | Performs GO enrichment (BP, CC, MF) for each cluster using `clusterProfiler`; creates bubble plots for top terms. | `*_cluster_markers.txt`<br>`org.Brapa.eg.db` package | `*_cluster_marker_GOres.txt` (combined results)<br>Bubble plots for BP/CC/MF |
| `datasets_correlation_anno.R` | Compares cluster average expression with external expression data from PNAS LCM (GSE12404) and NaturePlant RNA‑seq (GSE157145), with and without ComBat batch correction. | `*_clusters_res1_avg.txt`<br>External expression files<br>ID mapping files (`atha_to_bchi_*.txt`) | Correlation heatmaps, hierarchical clustering dendrograms, PCA biplots |

---

## Step‑by‑Step Workflow

### 1. Visualise marker genes (`marker_gene_anno.R`)

This script assumes you have already performed integration (e.g., Harmony) and clustering. It loads the metadata and Seurat object, then:

- Plots UMAP coloured by cluster and by original sample (with labels).
- Plots a stacked bar chart showing the proportion of each sample per cluster.
- Generates numerous **DotPlot** panels for various marker gene sets:
  - **PlantscRNAdb** – embryo, seed coat, endosperm.
  - **Region‑specific** – pg, g, h stages (embryo, endosperm, seed coat, etc.).
  - **Subregion‑specific** – finer subdivisions (EP, MCE, PEN, CZE, CZSC, SC, SUS).
  - **qPCR‑validated** – short lists of experimentally confirmed markers.
  - **Cell cycle** – G0, G1, G1/S, S, G2, M phase markers.
  - **FISH‑validated** – markers from in situ hybridisation.

These dot plots help identify which clusters express known cell‑type markers, guiding initial annotation.

**Run**:
```r
source("marker_gene_anno.R")
```

---

### 2. GO enrichment analysis (`diff_gene_GOres_anno.R`)

After obtaining cluster markers (from `FindAllMarkers`), this script:

- Filters markers by adjusted p‑value < 0.01 and |log₂FC| > 1.
- For each cluster (1–22), runs `enrichGO` using the `org.Brapa.eg.db` package (for *B. rapa*).
- Saves all results to `TT9D_cluster_marker_GOres.txt`.
- Reads back the combined file and creates three separate bubble plots for **B**iological **P**rocess, **C**ellular **C**omponent, and **M**olecular **F**unction, showing the top 2 enriched terms per cluster.

The bubble plots show GeneRatio (x‑axis), term description (y‑axis), point size = gene count, and colour = adjusted p‑value. Cluster numbers are labelled on each point.

**Run**:
```r
source("diff_gene_GOres_anno.R")
```

**Dependencies**: `clusterProfiler`, `org.Brapa.eg.db`.

---

### 3. Cross‑dataset correlation (`datasets_correlation_anno.R`)

To assess whether our clusters correspond to known cell types from other species, we compare our cluster‑averaged expression with two external *Arabidopsis* datasets:

- **PNAS LCM** (GSE12404) – laser‑capture microdissection of seed compartments.
- **NaturePlant RNA‑seq** (GSE157145) – bulk RNA‑seq from seed regions.

The script:

- Converts *Arabidopsis* gene identifiers to *B. rapa* orthologs using pre‑computed mapping files (provided).
- Filters to common genes and merges the matrices.
- Scales the combined matrix and computes Spearman correlations, hierarchical clustering, and PCA.
- Repeats the analysis **after ComBat batch correction** to remove dataset‑specific technical variation.

Visual outputs include heatmaps (cross‑correlation between our clusters and external samples), dendrograms, and PCA plots with sample labels.

**Run**:
```r
source("datasets_correlation_anno.R")
```

**Dependencies**: `sva` (for ComBat), `factoextra`, `ggpubr`.

---

## Required Input Files

| File | Description |
|------|-------------|
| `*_metadata.txt` | Cell metadata with UMAP coordinates and cluster assignments. |
| `*_merge.rds` | Seurat object containing the integrated data. |
| `*_cluster_markers.txt` | Output from `FindAllMarkers` (used for GO). |
| `*_clusters_avg.txt` | Average expression per cluster (used for correlation). |
| `GSE12404_expr_matrix_normalized.txt` | External LCM expression matrix (PNAS). |
| `GSE157145_CPM_average_expression.txt` | External RNA‑seq CPM matrix (NaturePlant). |
| `atha_to_bchi_probe_80.txt` | Mapping from Arabidopsis probes to *B. rapa* genes (for LCM). |
| `atha_to_bchi_besthit_80.txt` | Mapping from Arabidopsis gene IDs to *B. rapa* genes (for RNA‑seq). |

All paths in the scripts may need to be adjusted to match your directory structure.

---

## Notes and Tips

- **Order of execution**: It is recommended to run `visualize_markers.R` first to get an initial idea of cell identities, then `GO_enrichment.R` to confirm functional signatures, and finally `cross_dataset_correlation.R` to validate assignments with external data.
- **Cluster numbering**: The scripts assume clusters are numbered 0–27 (or 1–22 for GO). Adjust loop ranges if your clustering resolution differs.
- **ID mapping**: The cross‑dataset correlation relies on accurate orthology mapping. The provided mapping files were generated with 80% sequence identity cutoffs; you may need to update them if using a different genome version.
- **Batch correction**: ComBat correction is applied within the script, but you can toggle the `if(TRUE)` blocks to run with or without correction for comparison.
- **Dot plot sizes**: For long gene lists, the text size is reduced (`size = 2`) to avoid over‑crowding. Increase it if your screen resolution allows.

---

## Dependencies

All scripts require the following R packages:

```r
install.packages(c("dplyr", "Seurat", "ggplot2", "patchwork", "svglite", 
                   "sva", "ggsci", "ggrepel", "tidyr", "stringr", 
                   "clusterProfiler", "org.Brapa.eg.db", "factoextra", 
                   "pheatmap", "ggpubr"))
```

`org.Brapa.eg.db` is a custom package built from eggNOG‑mapper annotations (see the `setup_GO_database` directory).
