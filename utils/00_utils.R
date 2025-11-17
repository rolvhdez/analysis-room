# Helper functions -----
check_out <- function(x) if (!dir.exists(x)) dir.create(x)
fancy_process <- function(
  process,
  spinner_type = "simpleDotsScrolling",
  message = "Processing",
  ...
) {
  #' Creates an environment to be executed with a 
  #' spinner function to show progress for a process
  #'
  #' @param process
  #' @param message Message to be displayed when process is executed
  # Define the wrapper function
  wrapper <- function() {
    tryCatch({
      # Start message
      cli_process_start(message)
      # Execute the process
      result <- do.call(process, list(...))
      cli_process_done()
      return(result)
    }, error = function(e) {
      # Finish with error message
      cli_alert_danger(paste("Error:", e$message))
      stop(e)
    })
  }
  wrapper()
}
read_sumstats_file <- function(sumstats_path, chunk_size = 1000000) {
  #' Read the summary statistics file
  #' 
  #' By default read in batches of 1 million lines
  #' to reduce memory usage.
  #'
  #' @param sumstats_path
  #' @param chunk_size
  #'
  con <- file(sumstats_path, "r")
  df <- data.table::fread(text = readLines(con, n = chunk_size))
  while (TRUE) {
    chunk <- readLines(con, n = chunk_size)

    # When the number of lines left is zero, break
    if (length(chunk) == 0) break

    c <- data.table::fread(text = chunk)
    if (!identical(names(c), names(df))) {
      colnames(c) <- names(df)
    }
    df <- data.table::rbindlist(list(df, c))
  }
  close(con)
  df # Return
}

reformat_sumstats <- function(sumstats, model) {
  #' Change the table format to follow the template
  #' from https://r-graph-gallery.com/101_Manhattan_plot.html
  # 
  #' @param sumstats Raw GWAS summary statistics file
  #' @param model Software where it comes from
  
  require(dplyr)
  require(stats)
   
  # Column mappings for all available models
  models <- list(
    snipar = list(
      filter_col = "direct_log10_P",
      select_cols = c("chromosome", "pos", "SNP", "direct_log10_P", "direct_Beta", "freq", "direct_N"),
      new_names = c("CHR", "BP", "SNP", "P", "BETA", "MAF", "N")
    ),
    regenie = list(
      filter_col = "LOG10P",
      select_cols = c("CHROM", "GENPOS", "ID", "P", "BETA", "N"),
      new_names = c("CHR", "BP", "SNP", "P", "BETA", "N")
    )
  )
  if (!model %in% names(models)) {
    stop("Supported models are: ", paste(names(models), collapse = ", "))
  }
  spec <- models[[model]]
  x <- sumstats %>%
    dplyr::filter(!is.na(.data[[spec$filter_col]])) %>%
    dplyr::select(dplyr::all_of(spec$select_cols)) %>%
    stats::setNames(spec$new_names)
  # Apply model-specific transformations
  if (model == "snipar") {
    x <- x %>% dplyr::mutate(P = 10^(-P))
  }
  
  return(x)
}

ncbi_query <- function(gene_list){
  #' Get an NCBI query from Ensembl ID's using Entrez.
  #' Returns: Gene name, description, and summary
  #'
  #' @param gene_list List with Ensembl ID's

  # Create an empty data frame
  if (length(gene_list) == 0) {
    ncbi_annotations <- data.frame(
      GENE_ID = NA,
      GENE_NAME = NA,
      GENE_DESCRIPTION = NA,
      GENE_SUMMARY = NA
    )
    return(ncbi_annotations)
  }

  ncbi_annotations <- lapply(gene_list, function(gene_id){
    search_res <- entrez_search(
      db = "gene", term = paste0(gene_id, "[Ensembl ID]")
    )
    if(search_res$count == 0) return(NULL)
    summary <- rentrez::entrez_summary(db = "gene", id = search_res$ids)

    gene_name <- tryCatch({
      if(!is.null(summary$name)) {
        summary$name
      } else {
        NA_character_
      }
    }, error = function(e) NA_character_)

    gene_descriptor <- tryCatch({
      if (!is.null(summary$description)) {
        summary$description
      } else {
        NA_character_
      }
    }, error = function(e) NA_character_)

    gene_summary <- tryCatch({
      if(!is.null(summary$summary)) {
        summary$summary
      } else {
        NA_character_
      }
    }, error = function(e) NA_character_)

    data.frame(
      GENE_ID = gene_id,
      GENE_NAME = gene_name,
      GENE_DESCRIPTION = gene_descriptor,
      GENE_SUMMARY = gene_summary
    )
  })
  ncbi_annotations <- dplyr::bind_rows(ncbi_annotations)
  return(ncbi_annotations)
}

# Plot theme
theme_set(
  theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black",
                                 linewidth = 0.5),
      axis.line.y = element_line(color = "black",
                                 linewidth = 0.5),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(color = "#3d3d3d", size = 8),
      plot.caption = element_text(color = "#3d3d3d", size = 8),
      strip.text = element_text(color = "#3d3d3d", face = "bold", size = 12),
      strip.background = element_rect(
        color = "#3d3d3d", fill = "white", linewidth = 1
      )
    )
)

export_plot <- function(plot_obj, file_path, 
                       width = 1080, height = 1080 * 0.75, 
                       res = 150, units = "px") {
  #' Export a plot to file
  #'
  #' @param plot_obj The plot object to export
  #' @param plot_name Name of the plot (without extension)
  #' @param output_dir Output directory path
  #' @param width Plot width
  #' @param height Plot height
  #' @param res Plot resolution
  #' @param units Plot units
  
  # Export the plot
  png(
    filename = file_path,
    width = width, height = height,
    res = res, units = units
  )
  suppressMessages(print(plot_obj))
  dev.off()
  
  # Success message
  cli::cli_alert_success("Exported `{file_path}`")
  
  # Return the file path invisibly
  invisible(file_path)
}