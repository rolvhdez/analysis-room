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

qq_caption <- paste0(
  "No. variants: ", scales::comma(length(k_snps)), "\n",
  "Lambda inflation factor: ", round(lambda_value, 4)
)

qq_plot <- df_qqvalues %>%
  filter(observed != Inf) %>%
  ggplot(aes(x = theoretical, y = observed)) +
  geom_point(size = 1.33, shape = 1) +
  geom_smooth(method = "lm",
    linetype = "dashed", color = "red"
  ) +
  xlab(expression(Theoretical ~ -log[10](italic(p)))) +
  ylab(expression(Observed ~ -log[10](italic(p)))) +
  labs(
    title = phenotype,
    caption = qq_caption
  )

# Export
png(
  output_dir %&% "qqplot.png",
  width = 1080, height = 1080 * 0.75,
  res = 150, units = "px"
)
suppressMessages(print(qq_plot))
dev.off()
cli_alert_success("Exported `" %&% output_dir %&% "qqplot.png`")