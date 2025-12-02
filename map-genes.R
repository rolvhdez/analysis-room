# Author: Roberto Olvera Hernandez
# Date: 2025-12-02
#
# Based from:
# https://github.com/mcps-analysts/workflows/blob/main/gwas/topmed-imputed/04.2_manhattan-plot.R
#
# Description:
# This script maps genes in Biomart and NCBI's
# Entrez' databases from a sumstats file and return
# a tab separated text file.
#
# Usage:
# Rscript map-genes.R sumstats_file output_dir/

suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

# Load personalized functions
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")
source("utils/01_map_genes.R")

annotate_genes_to_sig_snps <- function(sumstats, gene_range, sig = 5e-8){
  #' Function obtained from MCPS' GitHub organization
  #'
  #' @param sumstats Data frame (formated) with GWAS summary statistics.
  #' @param gene_range List of genetic ranges.
  #' @param sig Significance threshold. Default is genome-wide significance.

  require(dplyr)
  sumstats_sig <- dplyr::filter(sumstats, P <= sig)

  # Create the genetic ranges
  gene_nearest <- c()
  for (i in 1:dim(sumstats_sig)[1]) {
    granges_sig <- GRanges(
      seqnames = sumstats_sig$CHR[i],
      ranges = IRanges(start = sumstats_sig$BP[i], end = sumstats_sig$BP[i])
    )
    # Find the nearest gene
    g <- gene_range[(nearest(granges_sig, gene_range))] %>% names(.)
    gene_nearest <- append(gene_nearest, g)
  }
  sumstats_sig$GENE <- gene_nearest
  out_df <- c()
  for (gene in unique(sumstats_sig$GENE)){
    sub_df <- filter(sumstats_sig, GENE == gene) %>% arrange(P)
    out_df <- rbind(out_df, sub_df[1,])
  }
  return(out_df)
}

# Process arguments
args <- commandArgs(trailingOnly = TRUE)
sumstats_file <- args[1]
output_dir <- args[2]
model <- if (length(args) < 3) "snipar" else args[3] # Default: snipar (v0.0.22)
compute_bonferroni <- if (length(args) < 4) "no" else args[4] # Default: empty

# Check that file and output directory exist
if (!file.exists(sumstats_file)) {
  cli_abort(c(
    "{sumstats_file} does not exist",
    "x" = "You've supplied a file that does not exist."
  ))
}
if (!dir.exists(output_dir)) {
  output_dir <- output_dir %&% "."
}
if (!model %in% c("snipar", "regenie")) {
  cli_abort(c(
    "{model} does not exist",
    "x" = "You've supplied an unsupported type of model. Options are: snipar (default), regenie"
  ))
}

#------------------------------------
# Read the summary statistics ---
raw_sumstats <- fancy_process(
  process = read_sumstats_file,
  message = "Reading " %&% sumstats_file,
  # Function parameters
  sumstats_path = sumstats_file,
  chunk_size = 1000000
)
df_sumstats <- reformat_sumstats(raw_sumstats, model) # Reformat the table
df_sumstats$CHR <- as.integer(df_sumstats$CHR)
k <- length(unique(df_sumstats$SNP)) # Number of lines
if (compute_bonferroni == "yes") {
  bonferroni <- 0.05 / nrow(df_sumstats)
} else {
  bonferroni <- 5e-8 # Bonferroni adjusted P-Value
}
cli_alert_info(scales::comma(k) %&% " SNPs found in `" %&% sumstats_file %&% "`.")
cli_alert_info("Bonferroni adjusted P-value: " %&% scales::scientific(bonferroni))
print(head(df_sumstats))

# Download the Ensembl variant annotations ---
ensembl_url <- "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"
ensembl_path <- "/tmp/Homo_sapiens.GRCh38.115.gtf.gz"
if(!file.exists(ensembl_path)){
  cli_alert_info("Ensembl data base not found, downloading...")
  download.file(
    ensembl_url,
    destfile = ensembl_path,
    method = "wget",
    extra = "-r -p --random-wait"
  )
}
cli_alert_info("Reading " %&% ensembl_path %&% "...")
df_ensembl <- as.data.frame(rtracklayer::import(ensembl_path))
df_genes <- df_ensembl %>%
  filter(
    type == "gene",
    gene_biotype == "protein_coding"
  )
gene_ranges <- GRanges(
  seqnames = df_genes$seqnames,
  IRanges(start = df_genes$start, end = df_genes$end)
)
names(gene_ranges) <- df_genes$gene_name
print(head(gene_ranges))

# Annotate suggested SNPs ---
df_sumstats_annotated <- annotate_genes_to_sig_snps(df_sumstats, gene_ranges, bonferroni)
output_file <- output_dir %&% "annotations.txt"
write.table(df_sumstats_annotate, 
            file = output_file,
            sep = " ",           # White space delimiter
            row.names = FALSE,   # Usually don't want row names
            col.names = TRUE,    # Include column names
            quote = FALSE)       # Don't quote strings