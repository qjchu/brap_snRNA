# CellChat Analysis

This directory contains scripts and resources for **inferring cell‑cell communication** from single‑cell RNA‑seq data using the [CellChat](https://github.com/sqjin/CellChat) framework. The workflow builds a customised ligand‑receptor database for *Brassica napus* (based on orthology with *Arabidopsis thaliana*), applies CellChat to each sample/timepoint, and generates summary plots to visualise communication patterns across cell types and conditions.

---

## Directory Contents

| File | Description |
|------|-------------|
| `LR_pair_atha_bchi.txt` | Tab‑separated file containing ligand‑receptor pairs from *Arabidopsis* (with *B. rapa* orthologs). Required to build the custom database. |
| `run_cellchat_analysis.R` | R script that builds the custom CellChat database and runs CellChat on a single sample (e.g., CK3D). For other samples, the same code is repeated with appropriate input files. |
| `plot_cellchat_summary.R` | Processes the `net_lr.txt` output files from all samples (CK and TT timepoints) and generates scatter plots showing for each cluster: the number of outgoing (source), incoming (target), and self‑signals. Points are coloured by major cell type (EP/SC/Endosperm) and faceted by timepoint. |
| `plot_selected_LR_pathways.R` | Filters the CellChat output to a curated set of ligand‑receptor interactions (e.g., `ALE1->GSO2`, `RALFL33->AGP16`, `LP1->CAM1`, etc.) and creates a dot plot showing their communication probability (point size) and significance (colour) across all timepoints and cell type pairs. |

---

## Workflow Overview

1. **Prepare the ligand‑receptor database** – `LR_pair_atha_bchi.txt` is read and formatted into a CellChat‑compatible list (`db`).
2. **Run CellChat** – For each sample, the script loads the Seurat object, creates a CellChat object, assigns the custom database, identifies overexpressed genes/interactions, computes communication probabilities, and exports:
   - `*_net_lr.txt` – all ligand‑receptor interactions with probabilities and p‑values.
   - `*_net_pathway.txt` – aggregated communication at the pathway level.
   - A saved CellChat object (`.rds`) for further exploration.
3. **Summarise communication counts** – `plot_cellchat_summary.R` reads all `net_lr.txt` files, calculates per‑cluster incoming/outgoing/self interaction counts, and produces a faceted scatter plot (similar to Figure 5B in the paper).
4. **Visualise selected pathways** – `plot_selected_LR_pathways.R` extracts a predefined set of important ligand‑receptor pairs and generates a dot plot showing their activity across time and cell‑type pairs.

---

## Input Requirements

- **Seurat objects** for each sample (e.g., `CK3D_merge_res1.rds`, `CK5D_merge_res1_subcluster.rds`, etc.) containing the `harmony_clusters_res1` metadata column (or appropriate cluster identities).
- The **`LR_pair_atha_bchi.txt`** file must be present in the specified directory.
- All R packages listed at the top of the scripts (`CellChat`, `Seurat`, `tidyverse`, `NMF`, `ggalluvial`, `patchwork`, `ggplot2`, `svglite`).

> **Note**: The provided `run_cellchat_analysis.R` script only includes code for **CK3D** as an example. To run it for other samples, copy the block and replace the RDS file path and output file names accordingly.

---

## Running the Analysis

### 1. Build database and run CellChat (per sample)

```bash
Rscript run_cellchat_analysis.R
```

This will produce:
- `CK3D_net_lr.txt` – ligand‑receptor table
- `CK3D_net_pathway.txt` – pathway summary
- `cellchat_CK3D.rds` – full CellChat object

Repeat for each sample (CK5D, CK7D, …, TT11D) by modifying the input RDS file and output names.

### 2. Generate summary scatter plot

```bash
Rscript plot_cellchat_summary.R
```

**Output**: `CK_cellchat.svg` and `CK_cellchat.png` (for CK samples) and similar for TT.

The plot shows each cluster as a point, with X = number of outgoing interactions, Y = number of incoming interactions, size = number of self‑interactions, and colour = major cell type (EP, SC, Endosperm). Facets represent different timepoints.

**Example image**:  
<img width="1122" height="461" alt="image" src="https://github.com/user-attachments/assets/4566d1d7-369a-486b-9620-60148b982144" />

*Each point is a cluster; larger circles indicate more self‑signalling; colours show tissue identity.*

### 3. Visualise selected ligand‑receptor pathways

```bash
Rscript plot_selected_LR_pathways.R
```

**Output**: `cellchat_pathway.svg` and `cellchat_pathway.png`.

The dot plot shows, for each selected interaction (facet), the communication probability (point size) and p‑value (colour) for each source‑>target pair across all timepoints.

**Example image**:  
<img width="678" height="927" alt="image" src="https://github.com/user-attachments/assets/26422250-b8fa-4f5b-bb6a-d8f5b79302fc" />

*Points represent a specific ligand‑receptor interaction between a source and target cell type; size = probability, colour = significance.*

---

## Customising the Analysis

- **Select different LR pairs**: Edit the `selected_LR` vector in `plot_selected_LR_pathways.R` (e.g., lines containing `c('ALE1->GSO2', ...)`).
- **Change cell type labels**: The recode mappings in the scripts are manually curated; update them if your cluster annotations change.
- **Use a different database**: By default, the scripts use the custom `db` built from `LR_pair_atha_bchi.txt`. You can also use CellChat's built‑in databases by replacing `cellchat@DB <- db` with `CellChatDB.mouse` or `CellChatDB.human`.
