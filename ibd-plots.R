# Author: Roberto Olvera Hernandez
# Date: 2025-10-21
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
ibd_path <- args[1]
output_dir <- args[2]

ibd <- data.table::fread(ibd_path, header = TRUE, sep = "\t")

# DECLARING FUNCTIONS
"%&%" <- function(a, b) paste0(a, b)

check_out <- function(x) if (!dir.exists(x)) dir.create(x)

siblings_ibd_plots <- function(x) {
  require(ggplot2)
  # Take two random siblings from the IBD file
  sib_sample <- x[x$id1 == sample(x$id1, 1), ]
  sib_sample <- sib_sample[sib_sample$id2 == sample(sib_sample$id2, 1)]

  # Get the IID's of the siblings
  sib1 <- unique(x$id1)
  sib2 <- unique(x$id2)

  # Make the plot
  ibd_colors <- c("IBD0" = "white", "IBD1" = "dodgerblue2", "IBD2" = "firebrick2")
  g <- ggplot(data = sib_sample) +
      geom_rect(aes(xmin = start, xmax = stop, ymin = 0, max = 0.9, fill = ibd), color = "black", linewidth = 0.80) +
      facet_grid(chr ~ .) +
      xlab("Position (Mbp)") +
      labs(title = paste0(sib1, "_", sib2), fill = "IBD Type") +
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

density_length_plot <- function(x) {
  g <- ggplot(data = x) +
    geom_density(aes(x = length, color = chr), linewidth = 1) +
    facet_wrap(. ~ chr) +
    guides(color = "none") +
    xlab("Segment length (cM)") +
    theme_minimal()
  return(g)
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
  plots[[1]] <- density_length_plot(data)
  for (i in 2:6) plots[[i]] <- siblings_ibd_plots(data)

  # Name the plots
  names(plots)[1] <- "chr_density.png"
  for (i in 2:6) {
    j <- i - 1
    names(plots)[i] <- paste0("ibd-map_FS", sprintf("%02d", j), ".png")
  }
  return(plots)
}

export_plot <- function(plot, file_name) {
  png(file_name, height = 1080, width = 1920, res = 150, units = "px")
  print(plot)
  dev.off()
}

# --- MAIN EXECUTION ---
check_out(output_dir)

# Same result as Map but with slightly different syntax
plots <- create_plots(ibd)
print(names(plots))
mapply(export_plot, plots, output_dir %&% names(plots))