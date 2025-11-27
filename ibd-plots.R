# Author: Roberto Olvera Hernandez
# Date: 2025-10-21
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# DECLARING FUNCTIONS
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")
source("utils/05_ibd_inference.R")

# --- MAIN EXECUTION ---
args <- commandArgs(trailingOnly = TRUE)
ibd_path <- args[1]
output_dir <- args[2]; check_out(output_dir)

ibd <- fancy_process(
  process = data.table::fread,
  message = "Reading " %&% ibd_path,
  # Function parameters
  file = ibd_path,
  header = TRUE,
  sep = "\t"
)

# IBD mosaic ----
ibd_colors <- c("IBD0" = "white", "IBD1" = "dodgerblue2", "IBD2" = "firebrick2")
list_plots <- fancy_process(
  process = create_plots,
  message = "Creating plots",
  # Function parameters
  data = ibd
)
mapply(export_plot, list_plots, output_dir %&% names(list_plots))