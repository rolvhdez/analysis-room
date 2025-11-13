# Author: Roberto Olvera Hernandez
# Date: 2025-10-21
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# DECLARING FUNCTIONS
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")

check_out <- function(x) if (!dir.exists(x)) dir.create(x)
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
histogram_length_plot <- function(x, chrom = 1:22) {
  x <- x[x$chr %in% chrom]
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

  # Histogram plot but for each chromosome
  start_index <- length(plots) + 1
  end_index <- start_index + 22 - 1
  for(i in start_index:end_index) {
    j <- i - start_index + 1
    chr_num <- sprintf("%02d", j)
    plots[[i]] <- histogram_length_plot(data, chrom = j)
    names(plots)[i] <- paste0("shared_length_chr", chr_num, "_histogram.png")
  }

  # Generate plots for 5 random sibling pairs
  num_fs_plots <- 5
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
export_plot <- function(plot, file_name) {
  tryCatch(
    {
      #' Try exporting the 
      png(file_name, height = 1080, width = 1920, res = 150, units = "px")
      print(plot)
      suppressMessages(invisible(dev.off()))
    },
    error = function(cond) {
      cli::cli_alert_danger("Something happened during exportation.")
      message(conditionMessage(cond))
      NA
    },
    finally = {
      cli::cli_alert_success("Exported " %&% file_name)
    }
  )
}

# --- MAIN EXECUTION ---
args <- commandArgs(trailingOnly = TRUE)
ibd_path <- args[1]
output_dir <- args[2]
check_out(output_dir)
ibd <- fancy_process(
  process = data.table::fread,
  message = "Reading " %&% ibd_path,
  # Function parameters
  file = ibd_path,
  header = TRUE,
  sep = "\t"
)

ibd_types <- ibd %>%
  mutate(IBDType = as.factor(IBDType)) %>%
  group_by(Chr, IBDType) %>%
  summarise(Count = n()) %>%
  group_by(Chr) %>%
  mutate(IBDProp = Count / sum(Count)) %>%
  ungroup()

print(ibd_types %>% filter(Chr == 22))

segment_number <- ibd %>%
  group_by(Chr) %>%
  summarise(Count = n()) %>%
  ggplot() +
  geom_bar(aes(x = Chr, y = Count), stat = "identity",
  color = "black", fill = "white") +
  ylab("No. segments") +
  scale_x_continuous(breaks = seq(1,22,1))
export_plot(segment_number, output_dir %&% "number_segments_chr.png")

ibd_colors <- c("IBD0" = "white", "IBD1" = "dodgerblue2", "IBD2" = "firebrick2")

ibd_prop_plot <- ibd_types %>%
  mutate(IBDType = as.factor(paste0("IBD", IBDType))) %>%
  ggplot(aes(x = Chr, fill = IBDType, y = IBDProp)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = ibd_colors, drop = FALSE) +
  labs(y = "Proportion within Chromosome", x = "Chromosome") +
  scale_x_continuous(breaks = seq(1,22,1))

export_plot(ibd_prop_plot, output_dir %&% "ibd_proportions_chr.png")
q()

list_plots <- fancy_process(
  process = create_plots,
  message = "Creating plots",
  # Function parameters
  data = ibd
)
# Export the plots
mapply(export_plot, list_plots, output_dir %&% names(list_plots))