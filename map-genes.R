# Author: Roberto Olvera Hernandez
# Date: 2025-10-13
#
# Description:
# This script maps genes in Biomart and NCBI's
# Entrez' databases from a sumstats file and return
# a tab separated text file.
#
# Usage:
# Rscript map-genes.R sumstats_file output_dir/

suppressPackageStartupMessages(library(rentrez))
suppressPackageStartupMessages(library(biomaRt))

# Process arguments
args <- commandArgs(trailingOnly = TRUE)
sumstats_file <- args[1]
output_dir <- args[2]
model <- if (length(args) < 3) "snipar" else args[3] # Default: snipar (v0.0.22)

# Load personalized functions
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")

# Read the data to be analyzed
df_sumstats <- fancy_process(
  process = read_sumstats_file,
  message = "Reading " %&% sumstats_file,
  # Function parameters
  sumstats_path = args[1],
  chunk_size = 1000000
)

# Change the table format to follow the template
# from https://r-graph-gallery.com/101_Manhattan_plot.html
if (model == "snipar") {
  fgwas_results <- df_sumstats %>% 
    dplyr::filter(!is.na(direct_log10_P)) %>% 
    dplyr::select(
      "CHR" = chromosome,
      "BP" = pos,
      "SNP" = SNP,
      "P" = direct_log10_P
    ) %>% 
    mutate(P = 10^(-P))
} else if (model == "regenie") {
  fgwas_results <- df_sumstats %>% 
    dplyr::filter(!is.na(LOG10P)) %>% 
    dplyr::select(
      "CHR" = CHROM,
      "BP" = GENPOS,
      "SNP" = ID,
      "P" = P
    )
}
fgwas_results$SNP <- gsub("GSA-", "", fgwas_results$SNP)
fgwas_results$CHR <- as.integer(fgwas_results$CHR)
k_snps <- unique(fgwas_results$SNP)
cli_alert_info(scales::comma(length(k_snps)) %&% " SNPs found in `" %&% sumstats_file %&% "`.")

# --- MAPPING GENES ---
source("utils/01_map_genes.R")