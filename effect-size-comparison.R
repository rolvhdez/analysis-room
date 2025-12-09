# Author:
# Date: 
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")
source("utils/01_map_genes.R")

# Parameter ---
args <- commandArgs(trailingOnly = TRUE)
snipar_sumstats <- args[1]
regenie_sumstats <- args[2]
output_dir <- args[3]
phenotype <- if (length(args) < 4) "" else args[4] # Default: empty

# Check parameters ----
if (!file.exists(snipar_sumstats)) {
  print(snipar_sumstats)
  # Check if the summary statistics exist
  cli_abort(c(
    "{snipar_sumstats} does not exist",
    "x" = "You've supplied a file that does not exist."
  ))
}
if (!file.exists(regenie_sumstats)) {
  print(regenie_sumstats)
  # Check if the summary statistics exist
  cli_abort(c(
    "{regenie_sumstats} does not exist",
    "x" = "You've supplied a file that does not exist."
  ))
}
if (!dir.exists(output_dir)) {
  # Provide a file in the directory for output
  output_dir <- output_dir %&% "."
}

#------------------------------------
# Read the summary statistics -------
df_snipar <- fancy_process(
  process = read_sumstats_file,
  message = "Reading " %&% snipar_sumstats,
  # Function parameters
  sumstats_path = snipar_sumstats,
  chunk_size = 1000000
)
df_snipar <- reformat_sumstats(df_snipar, "snipar")
df_snipar$MODEL <- "snipar"
#
df_regenie <- fancy_process(
  process = read_sumstats_file,
  message = "Reading " %&% regenie_sumstats,
  # Function parameters
  sumstats_path = regenie_sumstats,
  chunk_size = 1000000
)
df_regenie <- reformat_sumstats(df_regenie, "regenie")
df_regenie$SNP <- paste0("chr", df_regenie$SNP)
df_regenie$MODEL <- "regenie"
df_regenie <- df_regenie %>%
  dplyr::filter(SNP %in% unique(df_snipar$SNP))
#
df_sumstats <- dplyr::bind_rows(df_snipar, df_regenie)

bonferroni <- 5e-8
sig_k <- df_sumstats %>% filter(P <= bonferroni) %>% pull(SNP)
if (length(sig_k) > 0) {
  cli::cli_alert_warning(scales::comma(length(sig_k)) %&% " SNPs found at p <= " %&% bonferroni)
  genes <- create_gene_ranges()
  df_annotations <- annotate_genes_to_sig_snps(
    sumstats = df_sumstats,
    gene_range = genes,
    sig = bonferroni
  )
  df_annotations <- df_annotations %>% select(SNP, GENE)
} else {
  cli::cli_alert_warning("No significant SNPs were found at p <= " %&% bonferroni %&% ". Skipping annotation.")
}

# Correlations ---
df_wide <- df_sumstats %>%
  pivot_wider(
    id_cols = SNP,
    names_from = MODEL,
    values_from = BETA
  )
snp_regenie_significant <- df_sumstats %>%
  filter(P < 5e-08 & MODEL == "regenie") %>%
  pull(SNP)
snp_snipar_significant <- df_sumstats %>%
  filter(P < 5e-08 & MODEL == "snipar") %>%
  pull(SNP)
df_wide <- df_wide %>%
  mutate(
    regenie_significant = ifelse(SNP %in% snp_regenie_significant, TRUE, FALSE),
    snipar_significant = ifelse(SNP %in% snp_snipar_significant, TRUE, FALSE)
  )
df_wide <- df_wide %>%
  left_join(df_annotations, by = "SNP")

corr_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "gray", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linewidth = 0.5, linetype = "dashed") +
  geom_point(
    data = df_wide,
    aes(x = regenie, y = snipar),
    color = "gray", alpha = 0.25
  ) +
  geom_point(
    data = filter(df_wide, regenie_significant == TRUE),
    aes(x = regenie, y = snipar, color = "regenie"),
  ) +
  geom_point(
    data = filter(df_wide, snipar_significant == TRUE),
    aes(x = regenie, y = snipar, color = "snipar"),
  ) +
  geom_text_repel(
    data = filter(df_wide, !is.na(GENE) & snipar_significant == TRUE),
    aes(x = regenie, y = snipar, label = GENE),
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 16,
    size = 3
  ) +
  ylab(expression(delta[snipar])) +
  xlab(expression(beta[regenie])) +
  labs(
    title = phenotype
  ) +
  scale_color_manual(
    name = "Genome-wide significance (p < 5e-8):",
    values = c(
      "regenie" = blue,
      "snipar" = red
    )
  ) +
  theme(legend.position = "top")
export_plot(corr_plot, output_dir %&% "correlation_plot.png")