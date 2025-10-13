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
cli_alert_info("Reading " %&% sumstats_file %&% " ...")
chunk_size <- 1000000 # Read every 1,000,000 lines
con <- file(sumstats_file, "r")

df_sumstats <- data.table::fread(text = readLines(con, n = chunk_size))
while (TRUE) {
  chunk <- readLines(con, n = chunk_size)

  # When the number of lines left is zero, break
  if (length(chunk) == 0) break 

  c <- data.table::fread(text = chunk)
  if (!identical(names(c), names(df_sumstats))) {
    colnames(c) <- names(df_sumstats)
  }
  df_sumstats <- rbindlist(list(df_sumstats, c))
}
close(con)

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