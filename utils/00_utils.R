# Helper functions -----
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
query_biomart <- function(snps) {
  #' Obtain Ensemble ID's through BiomaRt
  #'
  #' @param snps List with rsid's to query
  human_variation <- biomaRt::useMart(biomart = "ENSEMBL_MART_SNP",
                                      dataset = "hsapiens_snp")
  snp_attributes <- c("refsnp_id",
                     "ensembl_gene_stable_id")
  query_mart <- biomaRt::getBM(attributes = snp_attributes,
                              filters = "snp_filter",
                              values = snps,
                              mart = human_variation)
  colnames(query_mart)[1] <- "SNP"

  ensemble_ids <- query_mart %>%
    filter(if_all(everything(), ~ .x != "")) %>%
    rename(GENE_ID = ensembl_gene_stable_id)

  ensemble_ids
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
  ncbi_annotations
}