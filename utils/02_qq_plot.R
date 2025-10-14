# QQ plot
get_qqvalues <- function(data, df = 1) {
  # Compute expected values
  p_values <- data$P
  chi <- qchisq(1 - p_values, df = df)

  # Lambda inflation correction
  lambda <- median(chi, na.rm = TRUE) / qchisq(0.5, df)
  p_values <- chi / lambda
  n <- length(p_values)
  probs <- (1:n - 0.5) / n
  theoretical <- qchisq(probs, df)
  df <- data.frame(
    observed = sort(p_values),
    theoretical = sort(theoretical)
  )
  return(list(qqvalues = df, lambda = lambda))
}
qq_results <- get_qqvalues(fgwas_results)
lambda_value <- qq_results$lambda
cli_alert_info("Lambda genetic inflation factor: " %&% round(lambda_value, 4))

# Make the plot
df_qqvalues <- qq_results$qqvalues
qq_plot <- df_qqvalues %>%
  ggplot(aes(x = theoretical, y = observed)) +
  geom_point(size = 2.33, alpha = 0.85, shape = 1) +
  geom_point(
    data = subset(df_qqvalues, observed >= -log10(bonferroni)),
    size = 2.33, alpha = 0.85, shape = 1, color = "red"
  ) +
  geom_smooth(
    method = "lm", linetype = "dashed", color = "blue",
    alpha = 0.85, linewidth = 1, se = FALSE
  ) +
  labs(
    x = "Theoretical -log10(P)",
    y = "Observeded -log10(P)",
    caption = "No. variants: " %&% scales::comma(length(k_snps)) %&% "\nLambda inflation: " %&% round(lambda_value, 4)
  ) +
  theme(
    plot.caption = element_text(hjust = 0),
    plot.caption.position = "plot",
  )

# Export
png(
  output_dir %&% "qqplot.png",
  width = 1920 * 2, height = 1920,
  res = 300, units = "px"
)
suppressMessages(print(qq_plot))
dev.off()
cli_alert_success("Exported `" %&% output_dir %&% "qqplot.png`")