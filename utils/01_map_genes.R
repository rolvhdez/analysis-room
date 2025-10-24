# You can find the functions within the wrapper functions
# in utils/00_utils.R

chr <- fgwas_results %>% filter(P <= bonferroni) %>% pull(CHR)
pos <- fgwas_results %>% filter(P <= bonferroni) %>% pull(BP)
regions <- paste0(chr, ":", pos, ":", pos)

# --- BIOMART QUERY ---
k_query <- nrow(fgwas_results %>% filter(P <= bonferroni))
biomart_message <- paste0(
  "Querying ", scales::comma(k_query),
  " significant SNPs ",
  "in the BiomaRt database."
)
human_variation <- fancy_process(
  process = biomaRt::useMart,
  message = "Connecting to BiomaRt",
  ###
  biomart = "ENSEMBL_MART_SNP",
  dataset = "hsapiens_snp",
  host = "https://www.ensembl.org"
)
query_mart <- fancy_process(
  process = biomaRt::getBM,
  message = biomart_message,
  ###
  attributes = c("refsnp_id", "chr_name", "chrom_start", "associated_gene"),
  filters = "chromosomal_region",
  values = regions,
  mart = human_variation
)
colnames(query_mart)[1] <- "SNP"
ensemble_ids <- query_mart %>%
  filter(if_all(everything(), ~ .x != "")) %>%
  rename(GENE_ID = associated_gene)

# --- NCBI ENTREZ QUERY ---
ncbi_gene_mappings <- fancy_process(
  process = ncbi_query,
  message = "Querying NCBI (Entrez) for extra information",
  # Function arguments
  gene_list = ensemble_ids %>% pull(GENE_ID) %>% unique()
)

# Hard-type in case it returns NULL values
ensemble_ids$SNP <- as.character(ensemble_ids$SNP)
ensemble_ids$GENE_ID <- as.character(ensemble_ids$GENE_ID)
ncbi_gene_mappings$GENE_ID <- as.character(ncbi_gene_mappings$GENE_ID)

# Combine all results
full_gene_map <- fgwas_results %>% 
  inner_join(ensemble_ids, ., by = "SNP") %>% 
  inner_join(ncbi_gene_mappings, ., by = "GENE_ID") %>% 
  arrange(CHR, BP)

# Export to text file (gzipped)
if (nrow(full_gene_map) > 0) {
  out_file <- output_dir %&% "gene_mappings.txt.gz"
  gz <- gzfile(out_file, "w")
  write.table(
    full_gene_map,
    file = gz,
    sep = " ",
    row.names = FALSE
  )
  close(gz)
  cli_alert_success("Exported `" %&% out_file %&% "`")
} else {
  cli_alert_warning('Queries returned no results. No file will be produced.')
}