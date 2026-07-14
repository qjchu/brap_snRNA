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

You can also download our newly built package from:

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
