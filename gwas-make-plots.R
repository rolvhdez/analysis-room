# Author: Roberto Olvera Hernandez
# Date: 2025-11-06

suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

# Load personalized functions
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")
source("utils/02_qq_plot.R")
source("utils/03_manhattan_plot.R")

# Process arguments
args <- commandArgs(trailingOnly = TRUE)

sumstats_file <- args[1]
output_dir <- args[2]
model <- if (length(args) < 3) "snipar" else args[3] # Default: snipar (v0.0.22)

compute_bonferroni <- if (length(args) < 4) "no" else args[4] # Default: empty
phenotype <- if (length(args) < 5) "" else args[5] # Default: empty
annotations <- if (length(args) < 6) NULL else args[6] # Default: empty

# Parameter error handling
if (!file.exists(sumstats_file)) {
  # Check if the summary statistics exist
  cli_abort(c(
    "{sumstats_file} does not exist",
    "x" = "You've supplied a file that does not exist."
  ))
}
if (!dir.exists(output_dir)) {
  # Provide a file in the directory for output
  output_dir <- output_dir %&% "."
}
if (!is.null(annotations)) {
  # Read the provided annotations table
  df_annotations <- fancy_process(
    process = data.table::fread,
    message = "Reading " %&% annotations,
    ###
    file = annotations,
    sep = " "
  )
}

#################################################################

# Read the data
df_sumstats <- fancy_process(
  process = read_sumstats_file,
  message = "Reading " %&% sumstats_file,
  # Function parameters
  sumstats_path = sumstats_file,
  chunk_size = 1000000
)

k_snps <- length(unique(df_sumstats$SNP)) # Number of lines
cli_alert_info(scales::comma(k_snps) %&% " SNPs found in `" %&% sumstats_file %&% "`.")

df_sumstats <- reformat_sumstats(df_sumstats, model) # Reformat the table
df_sumstats$CHR <- as.integer(df_sumstats$CHR)
if (compute_bonferroni == "yes") {
  bonferroni <- 0.05 / nrow(df_sumstats)
} else {
  bonferroni <- 5e-8 # Bonferroni adjusted P-Value
}
cli_alert_info("Bonferroni adjusted P-value: " %&% scales::scientific(bonferroni))

### QQ PLOT ###
qq_list <- get_qqvalues(df_sumstats)
pvalues <- qq_list[[1]]
lambda <- qq_list[[2]]
qq_plot <- make_qqplot(pvalues, phenotype, lambda)
cli_alert_info("Lambda genetic inflation factor: " %&% round(lambda, 4))
export_plot(qq_plot, output_dir %&% "qqplot.png")

### MANHATTAN PLOT ###
list_manhattan <- format_for_manhattan(df_sumstats)
df_manhattan <- list_manhattan[[1]]
df_axis <- list_manhattan[[2]]
if (!is.null(annotations)) {
  df_manhattan <- df_manhattan %>%
    left_join(
      df_annotations %>% select(GENE_ID, CHR, BP),
      by = c("CHR", "BP")
    )
}
manhattan_plot <- make_manhattan(df_manhattan, df_axis, phenotype, bonferroni)
if (!is.null(annotations)) {
  manhattan_plot <- manhattan_plot +
    geom_label_repel(
      data = subset(df_manhattan, !is.na(df_manhattan$GENE_ID)),
      aes(label = CHR %&% ":" %&% BP %&% " " %&% GENE_ID),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 16,
      size = 2.25
    )
}
export_plot(manhattan_plot, output_dir %&% "manhattan_plot.png")

### EFFECT SIZES ###
### (only snipar) ###
## Make effect sizes plot (only available for `snipar`)
#if (model == "snipar") source("utils/04_effect_sizes.R")
#