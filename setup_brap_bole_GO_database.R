# Set working directory (adjust path if needed)
# getwd()
# setwd('GO_database')

# Load required libraries
library(tidyverse)
library(stringr)
library(AnnotationForge)
library(dplyr)
# Other packages (KEGGREST, jsonlite, purrr, RCurl) are not needed after removing KEGG

# Read and combine two emapper annotation files
emapper1 <- read.table("GO_database/brap.out.emapper.annotations",
                       header = FALSE, sep = "\t", fill = TRUE, quote = "")
colnames(emapper1) <- c('query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs',
                        'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name',
                        'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module',
                        'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC',
                        'CAZy', 'BiGG_Reaction', 'PFAMs')

emapper2 <- read.table("GO_database/bole.out.emapper.annotations",
                       header = FALSE, sep = "\t", fill = TRUE, quote = "")
colnames(emapper2) <- colnames(emapper1)  # use same column names

emapper <- rbind(emapper1, emapper2)

# Extract gene information (GID and preferred name)
gene_info <- emapper %>%
  dplyr::select(GID = query, GENENAME = Preferred_name) %>%
  na.omit()

# Process GO annotations
gos <- emapper %>%
  dplyr::select(query, GOs) %>%
  na.omit()

gene2go <- data.frame(GID = character(), GO = character(), EVIDENCE = character())

for (row in 1:nrow(gos)) {
  the_gid <- gos[row, "query"][[1]]
  the_gos <- str_split(gos[row, "GOs"], ",", simplify = FALSE)[[1]]
  df_temp <- tibble(GID = rep(the_gid, length(the_gos)),
                    GO = the_gos,
                    EVIDENCE = rep("IDA", length(the_gos)))
  gene2go <- rbind(gene2go, df_temp)
}

# Remove entries with missing GO
gene2go$GO[gene2go$GO == "-"] <- NA
gene2go <- na.omit(gene2go)

# Build the OrgDb package
tax_id <- "37113712"
genus <- "Brassica"
species <- "rapaoleracea"

makeOrgPackage(
  gene_info = gene_info,
  go = gene2go,
  version = "0.0.1",
  maintainer = "Qinjie Chu <qinjiechu@zju.edu.cn>",   # Update with your name and email
  author = "Qinjie Chu <qinjiechu@zju.edu.cn>",       # Update with your name and email
  outputDir = ".",
  tax_id = tax_id,
  genus = genus,
  species = species,
  goTable = "go"
)

# Install the generated package
install.packages("./org.Brapaoleracea.eg.db/", repos = NULL, type = "sources")
