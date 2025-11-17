make_effectsizes_plot <- function(df, title, bonferroni) {
  k <- length(unique(df$SNP))
  caption <- paste0(
    "No. variants: ", scales::comma(k), "\n",
    "Bonferroni adjusted p-value < " %&% scales::scientific(bonferroni, digits = 4)
  )
  p <- ggplot() +
    # Significant variants
    geom_hline(aes(yintercept = -log10(bonferroni)), color = "red", linetype = "dashed") +
    geom_point(
      data = subset(df, df$P < bonferroni),
      aes(y = -log10(P), x = BETA, color = BETA),
      shape = 1
    ) +
    # Non-Significant variants
    geom_hline(aes(yintercept = -log10(bonferroni)), color = "red", linetype = "dashed") +
    geom_point(
      data = subset(df, df$P >= bonferroni),
      aes(y = -log10(P), x = BETA),
      shape = 1, color = "gray", alpha = 0.3
    ) +
    xlab(expression(Effect~size~(beta))) +
    ylab(expression(-log[10](italic(p)))) +
    labs(
      title = title,
      caption = caption
    ) +
    scale_color_gradientn(
      colors = hcl.colors(20, "Spectral", rev = TRUE)
    ) +
    theme(
      legend.position = "none"
    )
}
make_effectmaf_plot <- function(df, title, bonferroni) {
  k <- length(unique(df$SNP))
  caption <- paste0(
    "No. variants: ", scales::comma(k), "\n",
    "Bonferroni adjusted p-value < " %&% scales::scientific(bonferroni, digits = 4)
  )
  effects_maf_plot <- ggplot() +
  # Significant variants
    geom_hline(aes(yintercept = -log10(bonferroni)), color = "red", linetype = "dashed") +
    geom_point(
      data = subset(df, df$P < bonferroni),
      aes(x = MAF, y = -log10(P), color = BETA)
    ) +
  # Non-Significant variants
    geom_point(
      data = subset(df, df$P >= bonferroni),
      aes(x = MAF, y = -log10(P), color = BETA),
      alpha = 0.3, shape = 1
    ) +
  ylab(expression(-log[10](italic(p)))) +
  xlab("Allele frequency") +
  labs(
    title = title,
    caption = caption,
    color = expression(Effect~size~(beta))
  ) +
  scale_color_gradientn(
    colors = hcl.colors(20, "Spectral", rev = TRUE)
  )
}