# Figure Scripts

This directory contains R scripts used to generate the main and supplementary figures for the publication. Each script focuses on a specific visualisation task, from quality control to cell‑type annotation and differential comparisons between wild‑type (CK) and mutant (TT) samples.

---

## Script Overview

| Script | Figure(s) | Purpose |
|--------|-----------|---------|
| `plot_QC_violin.R` | Supplementary quality‑control violin plots | Generates flipped violin plots (boxplot + violin) showing the distribution of `nFeature_RNA` and `nCount_RNA` per sample for the CK dataset. |
| `plot_CK_marker_dotplot.R` | Main marker dot‑plot (CK time‑series) | Creates a faceted dot plot displaying the expression of conserved marker genes across clusters from all CK timepoints. Genes are grouped by cell type (EP, CZE, MCE, PEN, CZSC, SC‑ii, SC‑oi, SC‑SUS, etc.). |
| `compare_TT_vs_CK_celltypes.R` | Differential expression (TT vs CK) per cell type (Figure 4A) | Compares matched cell types between TT (mutant) and CK (wild‑type) at 7D. For each cell type, it calculates log₂ fold change and visualises genes with |logFC| > 2 as a column‑scatter plot, highlighting key genes. |
| `plot_hormone_gene_expression.R` | Hormone‑related gene expression trends (Supplementary) | For selected cell types (CZSC, SC‑SUS, EP, etc.), extracts expression of ABA, IAA, GA, and XTH gene sets, and plots their expression over time (CK vs TT) as boxplot + line plots. |

---

## Detailed Descriptions

### 1. `plot_QC_violin.R`

**Input**: Seurat object (RDS) of merged CK samples (e.g., `CK_10samples_batch3_metadata_res1.rds`).  
**Output**: SVG plots: `CK_vlnplot_nFeature_RNA.svg` and `CK_vlnplot_nCount_RNA.svg`.  
**Function**:

- Extracts the sample name from `orig.ident` (e.g., “CK3D_1” → “CK3D”).
- Generates a `VlnPlot` for `nFeature_RNA` and `nCount_RNA` (without points) using `Seurat`.
- Re‑plots the data with `ggplot2` as a flipped violin + boxplot, with custom colour palettes.
- Prints quantile summaries and cell counts per sample to the console.

**Example figure**:  
<img width="200" height="345" alt="image" src="https://github.com/user-attachments/assets/ab76ac52-b85f-43aa-97f1-5eb4660e3a25" />

*Flipped violin plots showing the distribution of detected genes (nFeature_RNA) and UMI counts (nCount_RNA) for each CK timepoint.*

---

### 2. `plot_CK_marker_dotplot.R`

**Input**:
- Marker gene tables (from `FindAllMarkers`) for each CK timepoint (annotated with Arabidopsis symbols, descriptions).
- Seurat objects (RDS) for each CK timepoint.  
**Output**: `CK_markers_dotplot_top3.svg` – a faceted dot plot.  
**Function**:

- For each cell type, intersects the top 100 marker genes across relevant clusters/timepoints to identify conserved markers.
- Selects 3 representative genes per cell type (30 genes total).
- Extracts DotPlot data from each Seurat object, combines them, and annotates each cluster with a simplified cell type label.
- Creates a faceted `ggplot` with:
  - X‑axis: marker genes (grouped by cell type).
  - Y‑axis: cluster IDs (labelled by cluster number + cell type).
  - Point size: percent expressed.
  - Point colour: scaled average expression.
  - Facets: cell type (rows) × timepoint (columns).
- Saves the figure.

**Example figure**:  
<img width="1174" height="499" alt="image" src="https://github.com/user-attachments/assets/12f90a89-1655-403f-b2ae-2318a57275d1" />

*Dot plot showing expression of 30 selected marker genes across CK clusters (Y‑axis) and timepoints (facets). Each point represents a gene‑cluster combination; colour indicates average expression, size indicates fraction of cells expressing the gene.*

---

### 3. `compare_TT_vs_CK_celltypes.R`

**Input**:
- Average expression matrices for TT7D and CK7D (`TT7D_combined_RNA_harmony_clusters_res1_5_avg.txt` and `CK7D_RNA_harmony_clusters_res1_avg.txt`).
- A gene annotation file (`besthit_anno.txt`) for symbol mapping.  
**Output**:
- `TT7D_VS_CK7D_celltype_fig4A_v2.svg/png` – the main figure.
- Supplementary tables (`TT7D_VS_CK7D_celltype_tableS7.txt`).  
**Function**:

- Defines a helper `celltype_DEG()` that, for a matched pair of cell types (one from TT, one from CK), computes log₂FC (TT/CK) for all genes expressed in both, keeping only those with |log₂FC| > 2.
- Runs this for 13 matched cell‑type pairs (Endosperm‑active, CZSC, CZSC‑phloem‑xylem, EP‑dividing‑cell, Endosperm‑PCD, SC‑oi, SC‑SUS, SC‑dividing‑cell, three different SC groups, SC‑ii, SC‑Endosperm).
- Combines results and creates a **horizontal column‑scatter plot**:
  - Background grey bars indicate the max/min logFC range per cell type.
  - Points represent individual genes, coloured by logFC (gradient from green to red).
  - Labels (symbols) are added for key genes of interest (e.g., `AGL62`, `ICE1`, `LEC1`, `TT2`, etc.) using `geom_text_repel`.
  - The Y‑axis lists cell types, and the X‑axis shows log₂FC.
- Output is saved in both SVG and PNG formats, and a full DEG table is exported.

**Example figure**:  
<img width="1596" height="533" alt="image" src="https://github.com/user-attachments/assets/a216d77f-bdfd-4ba5-9159-54bb594695bb" />

*Horizontal bar‑scatter plot comparing gene expression (log₂FC) between TT7D and CK7D for matched cell types. Each dot is a gene; red/green gradient indicates up‑regulation in TT vs CK. Key genes are labelled.*

---

### 4. `plot_hormone_gene_expression.R` (Supplementary hormone plots)

**Input**:
- Average expression matrices for all CK and TT timepoints.
- Pre‑defined Arabidopsis gene sets for ABA synthesis/degradation, IAA synthesis/degradation, GA synthesis/degradation, and XTH family.
- Orthology mapping file (`bchi_to_atha_besthit.txt`).  
**Output**: Multiple SVG files in `../revised_data/CK_TT_geneset/` (e.g., `CK_TT_CZSC_ABAsyn.svg`, `CK_TT_SUS_IAAsyn.svg`, etc.).  
**Function**:

- For a given cell type (e.g., CZSC, SC‑SUS, EP‑dividing‑cell, PEN, SC‑ii, SC‑oi, CZSC‑phloem‑xylem), the script selects the appropriate cluster columns from the expression matrices (both CK and TT).
- Uses `FindBchiHomo()` to map Arabidopsis gene IDs to *B. rapa* orthologs and extract their expression values.
- Melts the data and creates a **facet‑wrapped boxplot + line plot**:
  - X‑axis: individual timepoints (e.g., CK3D, CK5D, … TT11D).
  - Y‑axis: normalised expression.
  - Boxes show distribution across genes; points and connecting lines show the trajectory of each individual gene over time.
  - Facets separate CK and TT.
- Saves each gene‑set/cell‑type combination as a separate SVG file.

**Example figure**:  
<img width="230" height="193" alt="image" src="https://github.com/user-attachments/assets/671f639f-668c-4d1e-a70b-957aa539dbec" />

*Expression of ABA synthesis genes across CK (left) and TT (right) timepoints in the CZSC cell type. Boxes show distribution, lines show individual gene trends.*

---

## Dependencies

All scripts require the following core packages (among others, depending on the script):

```r
install.packages(c("dplyr", "Seurat", "ggplot2", "svglite", "ggsci",
                   "ggrepel", "tidyr", "stringr", "sva", "patchwork",
                   "ggpubr", "scales", "RColorBrewer", "reshape2"))
```

Additionally, `compare_TT_vs_CK_celltypes.R` uses `sva` (for ComBat batch correction) and `pheatmap` for optional correlation heatmaps.
