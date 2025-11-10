# Author: Roberto Olvera Hernandez
# Date: 2025-11-10

suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# Load personalized functions
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")
table_h2_logs <- function(log_pattern, output_file = NULL) {
  #' Read heritability log files and create a data frame
  #'
  #' @param log_pattern A file pattern to match log files (e.g., "*.log")
  #' @param output_file Optional: path to output file to save the table
  #' @return A data frame with heritability results
  #'
  #' @examples
  #' table_h2_logs("~/results/heritability/*.log")
  #' table_h2_logs("~/results/heritability/*.log", "~/results/heritability/table_h2.txt")
  #' 
  # Get list of log files
  log_files <- Sys.glob(log_pattern)
  
  if (length(log_files) == 0) {
    stop("No log files found matching the pattern: ", log_pattern)
  }
  
  # Initialize empty vectors to store results
  traits <- character()
  h2 <- numeric()
  h2_se <- numeric()
  lambda <- numeric()
  chi <- numeric()
  intercept <- numeric()
  intercept_se <- numeric()
  
  # Process each log file
  for (file in log_files) {
    # Read file content
    content <- readLines(file)
    
    # 1. Extract trait name from filename
    trait_name <- tools::file_path_sans_ext(basename(file))
    traits <- c(traits, trait_name)
    
    # 2. Extract values using regular expressions
    # h2 and h2_se
    h2_line <- grep('Total Observed scale h2:', content, value = TRUE)
    h2_values <- as.numeric(stringr::str_extract_all(h2_line, '[0-9]+\\.[0-9]+')[[1]])
    current_h2 <- ifelse(length(h2_values) >= 1, h2_values[1], NA)
    current_h2_se <- ifelse(length(h2_values) >= 2, h2_values[2], NA)
    
    # lambda
    lambda_line <- grep('Lambda GC:', content, value = TRUE)
    current_lambda <- ifelse(length(lambda_line) > 0, 
                            as.numeric(stringr::str_extract(lambda_line, '[0-9]+\\.[0-9]+')), 
                            NA)
    
    # chi
    chi_line <- grep('Mean Chi\\^2:', content, value = TRUE)
    current_chi <- ifelse(length(chi_line) > 0, 
                         as.numeric(stringr::str_extract(chi_line, '[0-9]+\\.[0-9]+')), 
                         NA)
    
    # intercept and intercept_se
    intercept_line <- grep('Intercept:', content, value = TRUE)
    intercept_values <- ifelse(length(intercept_line) > 0,
                              as.numeric(stringr::str_extract_all(intercept_line, '[0-9]+\\.[0-9]+')[[1]]),
                              NA)
    current_intercept <- ifelse(length(intercept_values) >= 1, intercept_values[1], NA)
    current_intercept_se <- ifelse(length(intercept_values) >= 2, intercept_values[2], NA)
    
    # Append to vectors
    h2 <- c(h2, current_h2)
    h2_se <- c(h2_se, current_h2_se)
    lambda <- c(lambda, current_lambda)
    chi <- c(chi, current_chi)
    intercept <- c(intercept, current_intercept)
    intercept_se <- c(intercept_se, current_intercept_se)
  }
  
  # Create data frame
  result_df <- data.frame(
    trait = traits,
    h2 = h2,
    h2_se = h2_se,
    lambda = lambda,
    chi = chi,
    intercept = intercept,
    intercept_se = intercept_se,
    stringsAsFactors = FALSE
  )
  
  # Save to file if output_file is specified
  if (!is.null(output_file)) {
    write.table(result_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    message("Done writing to ", output_file, " ...\n")
  }
  
  return(result_df)
}

# Process arguments ---
args <- commandArgs(trailingOnly = TRUE)
h2_files <- args[1]
output_dir <- args[2]

# Read the files ---
h2 <- table_h2_logs(h2_files)

# Create the plot

h2_plot <- ggplot(data = h2, aes(color = trait)) +
  geom_point(aes(x = h2, y = reorder(trait, h2)), size = 2, shape = 1) +
  geom_errorbarh(  # Use geom_errorbarh for horizontal error bars
    aes(xmin = h2 - h2_se, xmax = h2 + h2_se, y = reorder(trait, h2)),
    height = 0.2
  ) +
  xlab(expression(Mean~SNP~heritability~(h[SNP]^2))) +
  ylab("") +
  scale_x_continuous(labels = scales::percent) +
  theme(legend.position = "none")

png(
  output_dir %&% "h2_plot.png",
  width = 1080, height = 1080 * 0.75,
  res = 150, units = "px"
)
suppressMessages(print(h2_plot))
dev.off()
cli_alert_success("Exported `" %&% output_dir %&% "h2_plot.png`")