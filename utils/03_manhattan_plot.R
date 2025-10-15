df_manhattan <- fgwas_results %>%
  group_by(CHR) %>%
  summarise(CHR_LEN = max(BP)) %>%
  mutate(COORD = cumsum(as.numeric(CHR_LEN)) - CHR_LEN) %>%
  dplyr::select(-CHR_LEN) %>%
  right_join(fgwas_results, ., by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BP_CUM = COORD + BP)

df_axis <- df_manhattan %>%
  dplyr::select(CHR, BP_CUM) %>%
  group_by(CHR) %>%
  summarise(CENTER = min(BP_CUM) + (max(BP_CUM) - min(BP_CUM)) / 2)

if (!is.null(annotations)) {
  df_manhattan <- df_manhattan %>%
    left_join(
      df_annotations %>% select(GENE_ID, CHR, BP),
      by = c("CHR", "BP")
    )
}

# Manhattan plot ----
manhattan_caption <- paste0(
  "No. variants: ", scales::comma(length(k_snps)), "\n",
  "Bonferroni adjusted p-value < " %&% scales::scientific(bonferroni,
                                                          digits = 4)
)
manhattan_plot <- df_manhattan %>%
  ggplot(aes(x = BP_CUM, y = -log10(P))) +
  geom_point(
    aes(color = as.factor(CHR)),
    size = 1.3, shape = 1
  ) +
  # Significance line
  geom_hline( # High
    yintercept = -log10(bonferroni),
    color = "red",
    linetype = "dashed"
  ) +
  geom_hline( # Moderate
    yintercept = -log10(1e-6),
    color = "blue",
    linetype = "dashed"
  ) +
  # Axis labels
  xlab("Chromosome") +
  ylab(expression(-log[10](italic(p)))) +
  labs(
    title = phenotype,
    caption = manhattan_caption
  ) +
  # Modify axis
  scale_x_continuous(
    label = df_axis$CHR,
    breaks = df_axis$CENTER
  ) +
  scale_color_manual(
    values = rep(c("gray", "steelblue"), 22)
  ) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(face = "bold"),
    panel.spacing = unit(1.5, "lines"),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5, hjust = 1
    )
  )

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

png(
  output_dir %&% "manhattan_plot.png",
  width = 1080, height = 1080 * 0.75,
  res = 150, units = "px"
)
suppressMessages(print(manhattan_plot))
dev.off()
cli_alert_success("Exported `" %&% output_dir %&% "manhattan_plot.png`")