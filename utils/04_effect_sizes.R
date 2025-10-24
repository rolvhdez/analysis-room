df_direct_plot <- df_sumstats %>%
  mutate(label = ifelse(direct_log10_P > -log10(5e-8), "yes", "no")) %>%
  dplyr::select(
    SNP, direct_Beta, direct_N, direct_log10_P, freq
  )

if (!is.null(annotations)) {
  df_direct_plot <- df_direct_plot %>%
    right_join(ensemble_ids, ., by = "SNP") %>%
    right_join(ncbi_gene_mappings, ., by = "GENE_ID") %>%
    dplyr::select(
      SNP, direct_Beta, direct_N, direct_log10_P, freq, GENE_NAME, label
    ) %>%
    mutate(GENE_SNP = paste0(SNP, ": ", GENE_NAME))
}

dge_caption <- paste0(
  "No. variants: ", scales::comma(length(k_snps)), "\n",
  "Bonferroni adjusted p-value < " %&% scales::scientific(bonferroni,
                                                          digits = 4)
)
direct_effects_plot <- df_direct_plot %>%
  ggplot(aes(x = freq, y = direct_log10_P)) +
  # All points
  geom_point(aes(color = direct_Beta), size = 2.33, alpha = 0.25) +
  # Significant points
  geom_point(
    data = subset(df_direct_plot, df_direct_plot$label == "yes"),
    aes(color = direct_Beta),
    size = 2.33, alpha = 1
  ) +
  geom_hline(
    yintercept = 6, color = "blue", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(bonferroni), color = "red", linetype = "dashed"
  ) +
  ylab(expression(-log[10](italic(p)))) +
  xlab("Allele frequency (%)") +
  labs(
    color = "DGE (δ)",
    caption = dge_caption
  ) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(limits = c(0, 12)) +
  scale_color_gradientn(
    colors = hcl.colors(20, "Spectral", rev = TRUE)
  )

if (!is.null(annotations)) {
  direct_effects_plot <- direct_effects_plot +
    geom_label_repel(
      data = subset(df_direct_plot, df_direct_plot$label == "yes"),
      aes(label = GENE_SNP),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 32,
      size = 2.25
    )
}

png(
  output_dir %&% "direct_effects_plot.png",
  width = 1080, height = 1080 * 0.75,
  res = 150, units = "px"
)
suppressMessages(print(direct_effects_plot))
dev.off()
cli_alert_success("Exported `" %&% output_dir %&% "direct_effects_plot.png`")