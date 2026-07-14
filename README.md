# brap_snRNA

**brap_snRNA** is a resource and toolkit for downstream analysis of single‑nucleus RNA‑sequencing (snRNA‑seq) data from **oilseed rape (*Brassica napus*)**. It provides analysis scripts, visualisation pipelines, and a custom gene annotation database (`org.Brapa.eg.db`) tailored for *B. napus*, enabling a complete post‑alignment workflow from data loading to biological interpretation.

---

## Repository Structure

```
brap_snRNA/
├── downstream_analysis/          # Core downstream analysis scripts
│   ├── cellchat/                 # Cell‑cell communication (CellChat)
│   ├── fig_scripts/              # Scripts to reproduce publication figures
│   ├── monocle2/                 # Pseudotime trajectory (Monocle2)
│   ├── datasets_correlation.R    # Cross‑dataset correlation analysis
│   ├── diff_gene_GO.R            # Differential expression and GO enrichment
│   └── marker_gene_anno.R        # Marker‑based cell type annotation
├── org.Brapa.eg.db/              # R package for B. napus gene annotation
│   ├── R/                        # Package source code
│   ├── inst/extdata/             # External annotation data
│   ├── man/                      # Documentation
│   ├── DESCRIPTION               # Package metadata
│   └── NAMESPACE                 # Namespace declarations
├── get_multi_samples_rds.R       # Load and merge multiple sample RDS files
├── get_single_sample_rds.R       # Load a single sample RDS file
└── setup_brap_bole_GO_database.R # Configure B. napus GO database
```

---

## Typical snRNA‑seq Analysis Workflow

The scripts in this repository are organised to follow a standard snRNA‑seq downstream analysis pipeline. The recommended order of use is:

1. **Data loading** – import your pre‑processed snRNA‑seq data (Seurat objects saved as RDS).
2. **Database setup** – configure the *B. napus* GO annotation and gene ID mapping.
3. **Cell clustering & visualisation** (basic Seurat steps, not included but expected before downstream).
4. **Marker gene annotation** – identify cell types using known markers.
5. **Differential expression & GO enrichment** – find DE genes and perform functional enrichment.
6. **Pseudotime analysis** – reconstruct developmental trajectories with Monocle2.
7. **Cell‑cell communication** – infer ligand‑receptor interactions with CellChat.
8. **Figure generation** – reproduce publication‑quality plots.

Below we describe each step and the corresponding scripts.

---

## Installation & Dependencies

### Requirements
- R (≥ 4.0)
- RStudio (recommended)

### Install Required R Packages
Before running the scripts, install the core dependencies:

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "monocle", "CellChat"))
```

### Install the *B. napus* Annotation Package
`org.Brapa.eg.db` provides Entrez‑to‑symbol mapping and GO annotations for *Brassica napus*. Install it from GitHub:

```r
# If devtools is not installed
install.packages("devtools")

devtools::install_github("qjchu/brap_snRNA/org.Brapa.eg.db")
```

---

## Step‑by‑Step Instructions

### Step 1 – Load Your Data
Use the provided scripts to import your snRNA‑seq data.

- **Single sample**: `get_single_sample_rds.R` loads one Seurat object (RDS).
- **Multiple samples**: `get_multi_samples_rds.R` loads and merges several samples into a combined Seurat object.

These scripts assume that your data have already been pre‑processed (quality control, normalisation, dimension reduction) and saved as Seurat objects.

---

### Step 2 – Set Up the GO Database
Before performing enrichment analyses, run the GO configuration script:

```r
source("setup_brap_bole_GO_database.R")
```

This initialises the *B. napus*‑specific GO annotation environment, which is required by `diff_gene_GO.R` and other functional analysis scripts.

---

### Step 3 – Cell Clustering and Dimensionality Reduction (Optional)
The repository does **not** include basic Seurat clustering scripts, but the standard workflow would involve:

```r
# Example – not provided, but expected prior to downstream steps
library(Seurat)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)
```

You can then use `DimPlot()` to visualise clusters.

---

### Step 4 – Cell Type Annotation Using Marker Genes
Run `downstream_analysis/marker_gene_anno.R` to assign cell identities based on known marker genes for *B. napus*. This script typically:

- Compares expression of established marker genes across clusters.
- Outputs a preliminary cell‑type annotation table that can be used in subsequent analyses.

---

### Step 5 – Differential Expression and GO Enrichment
Use `downstream_analysis/diff_gene_GO.R` to:

- Identify marker genes for each cluster (via `FindAllMarkers` or similar).
- Perform Gene Ontology enrichment on the DE gene lists.
- Visualise enriched GO terms.

The script relies on the GO database set up in Step 2.

---

### Step 6 – Pseudotime Trajectory Analysis (Monocle2)
The `downstream_analysis/monocle2/` directory contains scripts to:

- Convert Seurat objects into Monocle2 `CellDataSet` objects.
- Order cells along developmental trajectories.
- Identify genes that change along pseudotime.

Check the scripts inside that folder for usage details.

---

### Step 7 – Cell‑Cell Communication (CellChat)
The `downstream_analysis/cellchat/` folder provides scripts to:

- Construct a CellChat object from your Seurat data.
- Infer ligand‑receptor interactions between cell types.
- Visualise communication networks and signalling pathways.

---

### Step 8 – Reproduce Publication Figures
Use the scripts in `downstream_analysis/fig_scripts/` to generate the exact figures from the associated publication. For example, `fig5B_cellchat_plot_v3.R` reproduces a specific CellChat visualisation.

Additionally, `downstream_analysis/datasets_correlation.R` can be used to assess the correlation between independent datasets, which is useful for validating reproducibility.

---

## Additional Notes

- The provided scripts are designed for **B. napus** snRNA‑seq data. They may require adjustments if your data originate from a different species or tissue.
- The `org.Brapa.eg.db` package is built from a specific reference genome version; verify that your gene IDs match the annotation used.
- For large datasets, consider running analyses on a high‑performance computing cluster.

---

## Contact & Citation

If you use this repository in your research, please cite the original publication (details to be added) and provide a link to this GitHub repository.

For questions or suggestions, please open an issue on the [GitHub repository](https://github.com/qjchu/brap_snRNA).
