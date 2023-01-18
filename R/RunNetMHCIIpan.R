#' @export
# THE MAJOR STEPS
# 1. Keeping valid peptides
# 2. Keeping valid alleles
# 3. Keeping valid output value types
# 4. Keeping valid number of threads
# 5. Doing paired analysis, if asked
# 5.1 Forcing long format
# 5.2 Checking for allele and peptide vector lengths
# 5.3 Keeping only valid allele-peptide pairs
# 5.4 Collecting peptides per alleles in a list format
# 5.5 Running RunNetMHCpan function on the peptide subsets per alleles
# 5.6 Assessing the outputs

# 6. Writing temporary pepfile to its location
# 7. Generating commands to be run
# 8. Modifying commands before running
# 9. Running the predictions
# 9.a Getting back results into R
# 9.b Just running NetMHCpan 4.1 and leaving the result files
# 10. Notifying the user about the successful run

RunNetMHCIIpan <- function(alleles,
                           peptides,
                           paired_input = F,
                           value_type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"),
                           output_format = "long",
                           threads = 1,
                           result_files_location = NULL,
                           keep_pep = FALSE,
                           software_path = NULL,
                           tmppep_loc = NULL) {
  
  # setting tmppep_loc
  software_path_splt <- strsplit(software_path, "\\/")[[1]]
  if(is.null(tmppep_loc)) {tmppep_loc <- paste0(paste0(software_path_splt[-length(software_path_splt)], collapse = "/"), "/tmp.pep")}
  
  # filtering for valid peptides (only if it is not a paired analysis)
  if(!paired_input) {
    is_valid <- IsValidPeptide(peptides, hla_type = 2)
    invalid_peptides <- peptides[!is_valid]
    
    peptides <- peptides[is_valid]
    if(length(invalid_peptides) > 0) {
      msg <- paste("Skipped the following invalid peptides:", paste0(invalid_peptides, collapse = ", "))
      message(msg)
    }
  }
  
  # filtering for valid alleles (only if it is not a paired analysis)
  alleles_supported <- unlist(strsplit(readLines(paste0(c(software_path_splt[-length(software_path_splt)], "data", "allelelist.txt"), collapse = "/")), "\t"))
  if(!paired_input) {
    is_valid_allele <- alleles %in% alleles_supported
    invalid_alleles <- alleles[!is_valid_allele]
    alleles <- unique(alleles[is_valid_allele])
    if(length(invalid_alleles) > 0) {
      msg <- "Skipped the following unsupported alleles:"
      message(paste(msg, paste0(invalid_alleles, collapse = ", ")))
    }
  }
  
  # filtering for supported output value types
  bvaltypes_unsupported <- setdiff(value_type, c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"))
  if(!is.null(value_type) & all(value_type %in% bvaltypes_unsupported)) {
    msg <- "None of the binding value types are supported by netMHCIIpan 4.0."
    stop(msg)
  } else if (length(bvaltypes_unsupported) > 0) {
    msg <- paste0("The following binding value types are not supported by netMHCIIpan 4.0, they will be skipped: ",
                  paste0(bvaltypes_unsupported, collapse = ", "))
    warning(msg)
  }
  
  # if number of threads is higher than CPUs available, reset number of threads
  n_cores <- parallel::detectCores()
  if(threads > n_cores | threads < 1) {
    msg <- paste0("Resetting number of threads to be used (", threads,") to the number of available threads (", n_cores, ")")
    message(msg)
    threads <- n_cores
  }
  
  # performing analyses for paired input
  if(paired_input) {
    if(output_format == "wide") {
      msg <- "'wide' output format is not supported in paired mode, using 'long' output."
      warning(msg)
      output_format <- "long"
    }
    if(length(alleles) != length(peptides)) {
      msg <- "Pairwise prediction cannot be performed because length of allele vector is not equal to the length of peptide vector"
      stop(msg)
    }
    
    # storing the original input data before filtering for valid pairs
    # this is needed to return the corresponding binding values based on the original input of the function
    alleles_input <- alleles
    peptides_input <- peptides
    
    # the following snippet till the "checking if..." should be checked. what is the function of this?
    is_valid_allele <- alleles %in% alleles_supported
    invalid_alleles <- alleles[!is_valid_allele]
    alleles <- alleles[is_valid_allele]
    if(length(invalid_alleles) > 0) {
      msg <- paste("Skipped the following unsupported alleles:", paste0(invalid_alleles, collapse = ", "))
      message(msg)
    }
    
    # checking if all peptides are valid, keeping only valid pairs
    is_valid_peptide <- IsValidPeptide(peptides, hla_type = 2)
    is_valid_allele <- alleles %in% alleles_supported
    is_valid_pair <- is_valid_peptide & is_valid_allele
    
    invalid_alleles <- alleles[!is_valid_pair]
    invalid_peptides <- peptides[!is_valid_pair]
    alleles <- alleles[is_valid_pair]
    peptides <- peptides[is_valid_pair]
    if(any(!is_valid_pair)) {
      msg <- paste("Skipped the following invalid peptide-allele pairs:", paste0(invalid_peptides, "-", invalid_alleles, collapse = ", "))
      message(msg)
    }
    
    # collecting peptides per alleles
    peptides_per_alleles <- aggregate(peptides ~ alleles, FUN =  function(x) {x}, simplify = FALSE)
    peptides_per_alleles <- magrittr::set_names(peptides_per_alleles$peptides, peptides_per_alleles$alleles)
    
    if(threads == 1) {
      suppressWarnings(outlist <- purrr::imap(peptides_per_alleles, ~RunNetMHCIIpan(alleles = .y, peptides = .x,
                                                                                    value_type = value_type, output_format = output_format,
                                                                                    result_files_location = result_files_location,
                                                                                    threads = 1, keep_pep = keep_pep, software_path = software_path,
                                                                                    tmppep_loc = tmppep_loc)))
    } else {
      future::plan(multisession, workers = threads)
      suppressWarnings(outlist <- furrr::future_imap(peptides_per_alleles, ~RunNetMHCIIpan(alleles = .y, peptides = .x,
                                                                                         value_type = value_type, output_format = output_format,
                                                                                         result_files_location = result_files_location,
                                                                                         threads = 1, keep_pep = keep_pep, software_path = software_path,
                                                                                         tmppep_loc = paste0(tmppep_loc, ".", gsub("\\:", "-", .y)))))
      future:::ClusterRegistry("stop")
    }
    
    # assembling output
    outdf <- magrittr::set_rownames(do.call(rbind.data.frame, outlist), NULL)
    outdf <- outdf[fastmatch::fmatch(paste0(gsub("\\*", "", alleles_input), ".", peptides_input),
                                     paste0(gsub("\\*", "", outdf$allele), ".", outdf$peptide)), ]
    return(outdf)
  }
  
  # writing temporary pepfile to its location
  if(!paired_input) {write(peptides, file = tmppep_loc)}
  
  # generating commands to run
  cmds <- sapply(alleles, function(allele) {
    paste0(software_path, " -inptype 1 -f ", tmppep_loc  ," -BA -a ", allele, sep = "")
  }, USE.NAMES = F)
  
  # if a result file directory is specified
  if(!is.null(result_files_location)) {
    suppressWarnings(dir.create(result_files_location, recursive = TRUE))  # creating result files directory
    alleles_underscore <- gsub("\\:", "_", alleles) # windows cannot handle colons in filenames
    alleles_underscore <- gsub("\\-", "_", alleles_underscore)
    if(is.null(value_type)) {
      cmds <- paste0(cmds, " > ", result_files_location, "/", alleles_underscore, ".txt")
    } else {
      cmds <- paste0(cmds, " | tee ", result_files_location, "/", alleles_underscore, ".txt")
    }
  }
  
  # RUNNING PREDICTIONS - and optionally collecting binding values
  if(!is.null(value_type)) {
    results <- RunCommand(cmds, threads = threads, intern = T)
    if(!keep_pep) {file.remove(tmppep_loc)} # removing temporary pepfile
    outobj <- CollectBindingResults(results, value_type = value_type, output_format = output_format, hla_type = 2)
    
    # stroring version number
    predictor.version <- unlist(strsplit(results[[1]][grepl("# NetMHCIIpan version", results[[1]])], " "))
    attr(outobj, "predictor_version") <- predictor.version[length(predictor.version)]
    attr(outobj, "HLA_type") <- "II"
    
    return(outobj)
  } else {
    RunCommand(cmds, threads = threads, intern = F)
    if(!keep_pep) {file.remove(tmppep_loc)} # removing temporary pepfile
  }
  
  # printing final message
  if(is.null(value_type)) {
    msg <- paste0("Predictions are done, see result files in ", result_files_location)
    message(msg)
  }
}

