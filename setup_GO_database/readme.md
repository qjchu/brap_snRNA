# GO Database Construction for *Brassica napus*

This directory contains the workflow to build a custom **Gene Ontology (GO) annotation database** (R package) for *Brassica napus* (and related *B. rapa*/*B. oleracea* accessions) using annotations from **eggNOG‑mapper**. The resulting R package (`org.Brapaoleracea.eg.db`) can be installed and used with Bioconductor‑style annotation tools (e.g., `clusterProfiler`).

---

## Workflow Overview

The construction follows three main steps:

1. **Generate annotation files** – use the eggNOG‑mapper web server to produce `*.out.emapper.annotations` for your protein sequences.
2. **Run the R script** – parse the annotation files, extract GO terms and gene names, and build the R package.
3. **Install the package** – install the generated package locally for use in downstream analyses.

---

## Step 1 – Generate eggNOG‑mapper Annotations

- Go to the [eggNOG‑mapper web server](http://eggnog5.embl.de/#/app/home).
- Upload your protein FASTA files (for *B. napus* or related species).
- Select **"GOs"** and **"Preferred_name"** as output options (these are required).
- Run the job and download the resulting `*.out.emapper.annotations` file.
- Repeat for both *B. rapa* and *B. oleracea* accessions (if needed) and name them:
  - `brap.out.emapper.annotations`
  - `bole.out.emapper.annotations`

Place both files in a directory (e.g., `../revised_data/GO_database/`) relative to the R script.

> **Note**: The script expects the files to have the standard column order from eggNOG‑mapper v5. If your output has different columns, adjust the `col_names` vector accordingly.

---

## Step 2 – Build the R Package with the Script

The R script `setup_brap_bole_GO_database.R` performs the following:

- Reads the two annotation files (brap + bole) and combines them.
- Extracts gene identifiers (`query`) and preferred gene names (`Preferred_name`).
- Parses comma‑separated GO terms into a long‑format table with evidence code `IDA`.
- Builds an **AnnotationForge** package using `makeOrgPackage()`.
- Outputs the package source directory `org.Brapaoleracea.eg.db` in the current working directory.

After execution, you will see a folder named `org.Brapaoleracea.eg.db` in the current directory.

---

## Step 3 – Install the R Package

Install the newly built package from the source directory:

```r
install.packages("./org.Brapaoleracea.eg.db/", repos = NULL, type = "source")
```

Once installed, you can load it and use it with `clusterProfiler`, `AnnotationDbi`, etc.

```r
library(org.Brapaoleracea.eg.db)
columns(org.Brapaoleracea.eg.db)   # see available annotation types
```

### Pre‑built Packages (Baidu Cloud)

If you prefer not to build the package yourself, you can download the pre‑built versions directly from Baidu Cloud:

| Package | Description | Baidu Cloud Link | Password |
|---------|-------------|------------------|----------|
| `org.Brapa.eg.db` | GO annotation database for *B. rapa* (Chinese cabbage) genome only | [https://pan.baidu.com/s/1fUUyDEbTNsqNxmxjrrah5w?pwd=brap](https://pan.baidu.com/s/1fUUyDEbTNsqNxmxjrrah5w?pwd=brap) | `brap` |
| `org.BrapoleKEGG.eg.db` | KEGG pathway database for the merged *B. rapa* + *B. oleracea* genome | [https://pan.baidu.com/s/1tuuGfHcn0lJxeE3mnbkyIA?pwd=brap](https://pan.baidu.com/s/1tuuGfHcn0lJxeE3mnbkyIA?pwd=brap) | `brap` |
| `org.Brapaoleracea.eg.db` | GO annotation database for the merged *B. rapa* + *B. oleracea* genome | [https://pan.baidu.com/s/1WfOAke-yrpLmCl4LefPqbQ?pwd=brap](https://pan.baidu.com/s/1WfOAke-yrpLmCl4LefPqbQ?pwd=brap) | `brap` |

**Package descriptions:**

- **`org.Brapa.eg.db`** – GO annotations based solely on the *Brassica rapa* (Chinese cabbage / 白菜) reference genome. Use this if your data align to the *B. rapa* genome only.
- **`org.BrapoleKEGG.eg.db`** – KEGG pathway annotations for the merged *B. rapa* + *B. oleracea* (cabbage / 甘蓝) genome. Suitable for pathway enrichment when working with the combined reference.
- **`org.Brapaoleracea.eg.db`** – GO annotations for the merged *B. rapa* + *B. oleracea* genome. Recommended for GO enrichment when using the merged reference.

After downloading the package source (`.tar.gz` or folder), install it locally:

```r
# For example, install the GO database for the merged genome
install.packages("/path/to/downloaded/org.Brapaoleracea.eg.db/", repos = NULL, type = "source")
```

Choose the package that best matches your reference genome and analysis needs.

---

## Customisation

- To include **KEGG pathway** annotations, uncomment the KEGG parsing block in the script and provide a KEGG hierarchy JSON file (e.g., from the KEGG API for your organism).
- To use only a single annotation file (e.g., only *B. napus*), modify the `files` vector accordingly.
- Adjust `tax_id`, `genus`, and `species` to match your organism’s NCBI taxonomy.

---

## Dependencies

- R (≥ 4.0) with packages: `tidyverse`, `AnnotationForge`.
- eggNOG‑mapper output files must be present.

---

## References

- eggNOG‑mapper: [http://eggnog5.embl.de/](http://eggnog5.embl.de/)
- AnnotationForge: [Bioconductor](https://bioconductor.org/packages/AnnotationForge/)
