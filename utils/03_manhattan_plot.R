# Input values
bp_window <- 10e4

# Select SNPs by association and distance (will have repeated rows) ----
snps_of_interest <- fgwas_results %>%
  filter(P <= 1e-6) %>%
  dplyr::select(CHR, BP) %>%
  inner_join(fgwas_results, ., by = c("CHR"), relationship = "many-to-many") %>%
  mutate(POS_DIST = abs(BP.x - BP.y)) %>%
  dplyr::select(-BP.y) %>%
  rename(BP = BP.x)

# Define SNP list for Manhattan plot ----
snps_label <- snps_of_interest %>%
  filter(POS_DIST == 0) %>%
  dplyr::select(SNP) %>%
  unique()
snps_significant <- snps_of_interest %>%
  filter(P < 5e-8) %>%
  filter(POS_DIST <= bp_window | POS_DIST == 0) %>%
  dplyr::select(SNP) %>%
  mutate(is_significant = "yes") %>%
  unique()
snps_suggested <- snps_of_interest %>%
  filter(P >= 5e-8 & P < 1e-6) %>%
  filter(POS_DIST <= bp_window | POS_DIST == 0) %>%
  dplyr::select(SNP) %>%
  mutate(is_suggested = "yes") %>%
  unique()

# --- DATA FRAME: MANHATTAN PLOT
df_manhattan <- fgwas_results %>%
  group_by(CHR) %>%
  summarise(CHR_LEN = max(BP)) %>%
  mutate(COORD = cumsum(as.numeric(CHR_LEN))-CHR_LEN) %>%
  dplyr::select(-CHR_LEN) %>%
  right_join(fgwas_results, ., by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BP_CUM = COORD + BP)

df_manhattan <- df_manhattan %>% # add highlighted snps
  right_join(snps_significant, ., by = c("SNP")) %>%
  mutate(is_significant = ifelse(is.na(is_significant), "no", "yes")) %>%
  right_join(snps_suggested, ., by = c("SNP")) %>%
  mutate(is_suggested = ifelse(is.na(is_suggested), "no", "yes"))

# Gene mappings
# df_manhattan <- df_manhattan %>%
#   right_join(ensemble_ids, ., by = "SNP") %>%
#   right_join(ncbi_gene_mappings, ., by = "GENE_ID") %>%
#   mutate(GENE_SNP = paste0(SNP, ": ", GENE_NAME))

df_axis <- df_manhattan %>%
  dplyr::select(CHR, BP_CUM) %>%
  group_by(CHR) %>%
  summarise(CENTER = min(BP_CUM) + (max(BP_CUM) - min(BP_CUM)) / 2)

# Manhattan plot ----
manhattan_plot <- df_manhattan %>%
  ggplot(aes(x = BP_CUM, y = -log10(P)))+
  geom_point(aes(color = as.factor(CHR)), size = 1.3, alpha = 0.8) +
  # Highlight SNPs of interest
  geom_point(data = subset(df_manhattan, is_significant == "yes"),
             color = "green", size = 2.3) +
  geom_point(data = subset(df_manhattan, is_suggested == "yes"),
             color = "orange", size = 1.3) +
  # Significance lines
  geom_hline(yintercept = 6, color = "blue", linetype ="dashed")+
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed")+
  # Mapping labels
#  geom_label_repel(
#    data = subset(df_manhattan, df_manhattan$is_significant=="yes"),
#    aes(label=GENE_SNP),
#    box.padding = 0.5,
#    point.padding = 0.3,
#    max.overlaps = 16,
#    size = 2.25
#  ) +
  # Axis labels
  xlab("Chromosome") +
  ylab(expression(-log[10](italic(p)))) +
  labs(caption = "No. variants: " %&% scales::comma(length(k_snps)))
  # Modify axis
  scale_x_continuous(label = df_axis$CHR, breaks = df_axis$CENTER) +
  scale_color_manual(values = rep(c("gray", "steelblue"), 22)) +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(face = "bold"),
    panel.spacing = unit(1.5, "lines")
  )
png(
  output_dir %&% "manhattan_plot.png",
  width = 1920 * 2, height = 1920,
  res = 300, units = "px"
)
suppressMessages(print(manhattan_plot))
dev.off()
cli_alert_success("Exported `" %&% output_dir %&% "manhattan_plot.png`")