# Author: Roberto Olvera Hernandez
# Date: 2025-09-11
#
# Description:
# This script consolidates summary statistics files from
# snipar (https://github.com/AlexTISYoung/snipar), generating:
# - Manhattan plots for a single trait
# - QQ-Plots (with lambda inflation factor)
# - Genetic annotations using biomaRt and NCBI databases
# - Effect sizes plots for multiple traits
# - Read log files and return information in a table-like file
#
# Usage:
# Rscript gwas-make-plots.R sumstats output_dir
#
# Note: Required file name should have the next structure:
#
# `<name of file>.<estimator>.txt.gz`
#
# where `estimator` for the type of model that was used in snipar to generate those 
# summary statistics, and should have the `.txt.gz` as the default extension.
#
# Note: For aesthetic reasons, will only show up to the top 50 gene annotations
# among the most significantly-associated loci. But will save a text file with
# all gene annotations and their corresponding GWAS variant
# (i.e. top variant at the locus)

suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(rentrez))
suppressPackageStartupMessages(library(biomaRt))

# Process arguments
args <- commandArgs(trailingOnly = TRUE)
sumstats_file <- args[1]
output_dir <- args[2]

# Check that file and output directory exist
if (!file.exists(sumstats_file)) {
  cli_abort(c(
    "{sumstats_file} does not exist",
    "x" = "You've supplied a file that does not exist."
  ))
}
if (!dir.exists(output_dir)) {
  cli_abort(c(
    "{output_dir} does not exist",
    "x" = "You've supplied a directory that does not exist."
  ))
}

# Load personalized functions
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")

# Plot theme
theme_set(
  theme_bw() +
    theme( 
      panel.border = element_blank(),
      plot.title = element_text(face = "bold", size=16),
      plot.subtitle = element_text(color="#3d3d3d", size=12),
      plot.caption = element_text(color="#3d3d3d", size=10),
      strip.text = element_text(color="#3d3d3d", face="bold", size=12),
      strip.background = element_rect(color="#3d3d3d", fill="white", linewidth=1),
    )
)

# Read the data to be analyzed
df_sumstats <- read.table(gzfile(sumstats_file), header = TRUE)
df_sumstats$SNP <- gsub("GSA-", "", df_sumstats$SNP)
df_sumstats$chromosome <- as.integer(df_sumstats$chromosome)

k_snps <- unique(df_sumstats$SNP)
cli_alert_info(scales::comma(length(k_snps)) %&% " SNPs found in `" %&% sumstats_file %&% '`.')

# Change the table format to follow the template
# from https://r-graph-gallery.com/101_Manhattan_plot.html
fgwas_results <- df_sumstats %>% 
  dplyr::filter(!is.na(direct_log10_P)) %>% 
  dplyr::select(
    "CHR" = chromosome,
    "BP" = pos,
    "SNP" = SNP,
    "P" = direct_log10_P
  ) %>% 
  mutate(P = 10^(-P))

# Make the gene mappings
source("utils/01_map_genes.R")

# Make the QQplot
source("utils/02_qq_plot.R")

# Make the Manhattan plot
source("utils/03_manhattan_plot.R")

# Make effect sizes plot
source("utils/04_effect_sizes.R")