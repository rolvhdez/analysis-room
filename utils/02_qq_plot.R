# QQ plot
get_qqvalues <- function(data, df = 1) {
  # Compute expected values
  p_values <- data$P
  chi <- qchisq(1 - p_values, df = 1)

  # Lambda inflation correction
  lambda <- median(chi, na.rm = TRUE) / qchisq(0.5, 1)
  p_values <- chi/lambda
  n <- length(p_values)
  probs <- (1:n - 0.5) / n
  theoretical <- qchisq(probs, df = df)
  df <- data.frame(
    observed = sort(p_values),
    theoretical = theoretical
  )
  return(list(qqvalues = df, lambda = lambda))
}
qq_results <- get_qqvalues(fgwas_results)
lambda <- qq_results$lambda
cli_alert_info("Lambda genetic inflation factor: " %&% round(lambda, 4))

# Make the plot
df_qqvalues <- qq_results$qqvalues
qq_plot <- df_qqvalues %>%
  ggplot(aes(x = theoretical, y = observed)) +
  geom_point(size = 2.33, alpha = 0.85, shape = 1) +
  geom_smooth(
    method = "lm", linetype = "dashed", color = "red",
    alpha = 0.85, linewidth = 1, se = FALSE
  ) +
  labs(
    x = "Theoretical",
    y = "Observeded",
    caption = "No. variants: " %&% scales::comma(length(k_snps))
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