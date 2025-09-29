ensemble_ids <- fancy_process(
  process = query_biomart,
  message = "Querying BiomaRt for gene mappings",
  # Function arguments
  snps = fgwas_results %>% filter(P <= 1e-4) %>% pull(SNP) %>% unique()
)
ncbi_gene_mappings <- fancy_process(
  process = ncbi_query,
  message = "Querying NCBI (Entrez) for extra information",
  # Function arguments
  gene_list = ensemble_ids %>% pull(GENE_ID) %>% unique()
)

full_gene_map <- fgwas_results %>% 
  inner_join(ensemble_ids, ., by="SNP") %>% 
  inner_join(ncbi_gene_mappings, ., by="GENE_ID") %>% 
  arrange(CHR, BP)
cli_alert_success("Exported `" %&% output_dir %&% "gene_mappings.csv`")
write.csv(
  full_gene_map,
  file = output_dir %&% "gene_mappings.csv",
  row.names = FALSE
)

top5_mapping <- full_gene_map %>% 
  slice_min(order_by = P, n = 5) %>% 
  arrange(P)
cli_alert_success("Exported `" %&% output_dir %&% "top5_mappings.csv`")
write.csv(
  top5_mapping,
  file = output_dir %&% "top5_mappings.csv",
  row.names = FALSE
)