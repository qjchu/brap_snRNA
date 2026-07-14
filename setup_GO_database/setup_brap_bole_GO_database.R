###############################################################################
# Script: Build an organism-level annotation package (org.Brapaoleracea.eg.db)
#         from eggNOG-mapper annotation files for Brassica napus (combined 
#         from brap and bole outputs).
# Steps:
#   1. Read two emapper annotation files (brap and bole).
#   2. Extract gene information (GID and Preferred_name).
#   3. Parse GO annotations (split comma-separated GOs).
#   4. Build the annotation package using AnnotationForge::makeOrgPackage.
###############################################################################

# Load required libraries
library(tidyverse)      # includes dplyr, tidyr, etc.
library(AnnotationForge) # to build org packages

# Set working directory if needed (adjust to your path)
# setwd("/path/to/your/project")

# ------------------- 1. Read emapper annotation files -------------------
# Define file paths (modify as needed)
emapper_files <- c(
  "../revised_data/GO_database/brap.out.emapper.annotations",
  "../revised_data/GO_database/bole.out.emapper.annotations"
)

# Define column names (same for both files)
col_names <- c(
  'query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl',
  'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko',
  'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',
  'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'
)

# Read and combine both files
emapper <- map_dfr(emapper_files, ~ read.table(.x, header = FALSE, sep = "\t", 
                                               fill = TRUE, quote = "", 
                                               col.names = col_names))

# ------------------- 2. Prepare gene information -------------------
gene_info <- emapper %>%
  select(GID = query, GENENAME = Preferred_name) %>%
  na.omit()   # keep only rows with a Preferred_name

# ------------------- 3. Parse GO annotations -------------------
# Select only rows with non-empty GOs, then split comma-separated terms
gene2go <- emapper %>%
  select(GID = query, GOs) %>%
  na.omit() %>%
  separate_rows(GOs, sep = ",") %>%          # one GO per row
  filter(GOs != "-") %>%                     # remove placeholder
  mutate(EVIDENCE = "IDA") %>%               # evidence code (Inferred from Direct Assay)
  select(GID, GO = GOs, EVIDENCE)

# ------------------- 4. (Optional) KEGG pathway mapping -------------------
# This block is disabled (if(FALSE)) – kept as reference.
# If needed, uncomment and ensure you have the KEGG JSON file.
# 
# if(FALSE) {
#   gene2ko <- emapper %>%
#     select(GID = query, Ko = KEGG_ko) %>%
#     na.omit() %>%
#     mutate(Ko = gsub("ko:", "", Ko))
#   
#   # Function to parse KEGG hierarchy JSON (requires jsonlite and RCurl)
#   # update_kegg <- function(json = "brp00001.json") { ... }
#   # load("kegg_info.RData")
#   # gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>% ...
# }

# ------------------- 5. Build the annotation package -------------------
tax_id   <- "37113712"          # NCBI taxonomy ID for Brassica rapa (subsp. oleifera?)
genus    <- "Brassica"
species  <- "rapaoleracea"      # adjust if needed

makeOrgPackage(
  gene_info   = gene_info,
  go          = gene2go,
  # ko         = gene2ko,       # uncomment if KEGG data are available
  # pathway    = gene2pathway,
  version     = "0.0.1",
  outputDir   = ".",            # current directory
  tax_id      = tax_id,
  genus       = genus,
  species     = species,
  goTable     = "go"
)

# ------------------- 6. Install the generated package -------------------
install.packages("./org.Brapaoleracea.eg.db/", repos = NULL, type = "source")
