siblings_ibd_plots <- function(x) {
  require(ggplot2)
  # Take two random siblings from the IBD file
  sib_sample <- x[x$id1 == sample(x$id1, 1), ]
  sib_sample <- sib_sample[sib_sample$id2 == sample(sib_sample$id2, 1)]
  # Get the IID's of the siblings
  sib1 <- unique(sib_sample$id1)
  sib2 <- unique(sib_sample$id2)
  # Make the plot
  ibd_colors <- c("IBD0" = "white", "IBD1" = "dodgerblue2", "IBD2" = "firebrick2")
  g <- ggplot(data = sib_sample) +
      geom_rect(aes(xmin = start, xmax = stop, ymin = 0, max = 0.9, fill = ibd), color = "black", linewidth = 0.80) +
      facet_grid(chr ~ .) +
      xlab("Position (Mbp)") +
      labs(title = paste0(sib1, ":", sib2), fill = "IBD Type") +
      scale_fill_manual(values = ibd_colors, drop = FALSE) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 8)
      )
  return(g)
}
chr_density_length_plot <- function(x) {
  g <- ggplot(data = x) +
    geom_density(aes(x = length, color = chr), linewidth = 1) +
    facet_wrap(. ~ chr) +
    guides(color = "none") +
    xlab("Segment length (cM)") +
    theme_minimal()
  return(g)
}
histogram_length_plot <- function(x) {
  median_length = median(x$length)
  g <- ggplot(data = x) +
    geom_histogram(
      aes(x = length),
      fill = "steelblue", color = "black",
      boundary = 0, binwidth = 50 / 10
    ) +
    geom_vline(
      xintercept = median_length,
      color = "red", linetype = "dashed"
    ) +
    xlab("Segment length (cM)") +
    ylab("No. segments") +
    labs(
      caption = paste0("Median segment length: ", round(median_length, 2), "cM")
    ) +
    scale_y_continuous(
      labels = scales::scientific
    ) +
    scale_x_continuous(
      breaks = seq(0, 250, 25)
    ) +
    theme_minimal()
}
create_plots <- function(data) {
  require(dplyr)
  require(ggplot2)
  plots <- list()
  data <- data %>%
    select( # simpler column names
      "chr" = "Chr", "ibd" = "IBDType", "length" = "length",
      "id1" = "ID1", "id2" = "ID2",
      "start" = "start_coordinate", "stop" = "stop_coordinate",
      "length" = "length"
    ) %>%
    mutate(
      # Make chr a discrete factor
      chr = as.factor(chr),
      # Change the level names to `IBD0`, `IBD1` and `IBD2`
      ibd = as.factor(paste0("IBD", ibd)),
      # Change to Mbp for more readable plots
      start = start / 1e6,
      stop = stop / 1e6
    )
  # Create the plots
  plots[[1]] <- chr_density_length_plot(data)
  names(plots)[1] <- "shared_length_chr_density.png"

  plots[[2]] <- histogram_length_plot(data)
  names(plots)[2] <- "shared_length_histogram.png"

  # Generate plots for 10 random sibling pairs
  num_fs_plots <- 10
  start_index <- length(plots) + 1
  end_index <- start_index + num_fs_plots - 1
  for(i in start_index:end_index) {
      j <- i - start_index + 1
      fs_num <- sprintf("%02d", j)
      plots[[i]] <- siblings_ibd_plots(data)
      names(plots)[i] <- paste0("ibd-mosaic_FS", fs_num, ".png")
  }
  return(plots)
}