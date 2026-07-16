# Cluster Correlation Plot

## Purpose

This step visualises the **relationships between cell clusters** in single‑cell data, including:

- The **proportion of cells** in each cluster (represented by node size).
- The **correlation strength** between clusters (represented by edge thickness).
- The **tissue identity** of each cluster (colour‑coded, e.g., endosperm, embryo‑proper, seed coat).

Such network plots help to intuitively understand developmental connections or transcriptional similarities between cell types, providing complementary evidence for cell‑type annotation and trajectory inference.

---

## Input Data

The script requires two Excel files (`.xlsx`) with the following formats:

### 1. Node file (`CK_clusters.xlsx`)

Defines the position, name, group, and size of each cluster (node).

| Column      | Description                                                         |
|-------------|---------------------------------------------------------------------|
| `newX`      | X‑coordinate of the node (the script automatically subtracts 0.5 for fine‑tuning) |
| `newY`      | Y‑coordinate of the node                                            |
| `group`     | Major tissue group (e.g., `Endosperm`, `EP`, `SC`), used for colour |
| `celltype3` | Detailed cell‑type label, displayed as a text annotation            |
| `value3`    | Numeric value used to control node size (usually cell proportion or count) |

**Example** (only a few rows shown):

| newX | newY | group      | celltype3               | value3    |
|------|------|------------|-------------------------|-----------|
| 1    | 1    | Endosperm  | Endosperm-PCD (c18)     | 0.0180    |
| 1    | 2    | Endosperm  | Endosperm-active (c24)  | 0.00486   |
| 1    | 3    | EP         | EP-dividing-cell (c15)  | 0.0310    |
| ...  | ...  | ...        | ...                     | ...       |

### 2. Edge file (`CK_clusters_cor.xlsx`)

Defines the connections (curves) between nodes. Each edge represents a correlation between two clusters.

| Column      | Description                                                         |
|-------------|---------------------------------------------------------------------|
| `startx`    | `newX` of the starting cluster (the script subtracts 0.5)           |
| `starty`    | `newY` of the starting cluster                                      |
| `endx`      | `newX` of the ending cluster                                        |
| `endy`      | `newY` of the ending cluster                                        |
| `curvature` | Curvature of the curve (positive/negative controls direction)       |
| `group`     | Grouping (not actually used in the current code; can be left empty) |
| `size2`     | Numeric value controlling line thickness (usually the correlation coefficient or strength) |

**Example**:

| startx | starty | endx | endy | curvature | group | size2   |
|--------|--------|------|------|-----------|-------|---------|
| 1      | 1      | 2    | 1    | 0.1       | A     | 0.909   |
| 1      | 1      | 2    | 2    | 0.1       | A     | 0.734   |
| ...    | ...    | ...  | ...  | ...       | ...   | ...     |

> Coordinates (`newX`/`newY`) are usually integers or decimals representing the layout. The script subtracts 0.5 from all `x` coordinates to centre the nodes.

---

## R Script Functionality

The script `plot_cluster_links.R` uses a custom geometry (`geom_curve2`, based on an extension of `ggplot2`) to draw curved connections with endpoints. Its main steps are:

1. Read node and edge data.
2. For each edge, keep only the strongest connection per start node (sorted by `size2` descending, take the first).
3. Plot:
   - **Curves**: `geom_curve2` draws arcs from start to end, with line width proportional to `size2`. Colour is determined by `group` (though not actively used in this version).
   - **Nodes**: `geom_point` draws circles, with size proportional to `value3` and colour based on `group` (Endosperm/EP/SC).
   - **Labels**: `annotate(geom="label")` displays the `celltype3` text.
4. Output as an SVG file.

The custom `geom_curve2` function is defined in `geom-curve2.R`; it automatically adds small circles at both ends of each curve to indicate direction.

---

## How to Run

1. **Prepare your data**: Place the node table and edge table as `CK_clusters.xlsx` and `CK_clusters_cor.xlsx` in the same directory (e.g., `../data/CK_batch3_10samples/`).
2. **Ensure required packages are installed**:
   ```r
   install.packages(c("tidyverse", "readxl", "ggplot2", "ggpubr", "grid", "ggsci", "ggnewscale"))
   ```
   Also ensure that `geom-curve2.R` is in the script directory (or sourced).
3. **Run the script**:
   ```r
   source("plot_cluster_links.R")
   ```
4. **Output**: An SVG file `CK_clusters_link.svg` is generated in the specified output directory.

---

## Example Output

The final graphic is a network plot with nodes positioned according to the coordinates, curved edges connecting correlated clusters, node size reflecting cell proportion, colour indicating tissue group, and labels showing detailed cell‑type names.

![Example network plot](CK_clusters_link.svg) *(schematic)*

---

## Parameter Tuning

- **Node size range**: Adjust `scale_size_continuous(range = c(0.5, 15))`.
- **Line width range**: Modify `scale_linewidth_continuous(range = c(0, 0.1))` to control curve thickness.
- **Colours**: Customise `scale_color_manual(values = c("Endosperm"="#e64b35", "EP"="#4dbbd5", "SC"="#00a087"))` as needed.

---

## Notes

- Ensure the node and edge coordinates are consistent (i.e., `startx/starty` match the corresponding `newX/newY` of nodes).
- For multiple time points or conditions, you can build separate network plots (this script only demonstrates the CK dataset).
- Curvature values >0 bend left, <0 bend right; adjust according to your layout preferences.

---

This step serves as a visual complement to cell‑type annotation, helping researchers quickly identify major cell populations and their interrelationships.
