### RunNetMHCIIpan function
### starts one or multiple threads of RunNetMHCIIpan processes to start epitope predictions

##################### SETTING INPUTS
configs <- list()
configs$software_loc_MHC_II <- "~/Programok/netMHCIIpan-4.0/netMHCIIpan"
configs$tmppep_loc_MHC_II <- "~/Programok/netMHCIIpan-4.0/tmp.pep"
peptides <- sapply(1:30, function(x) {sample(c("A", "E", "K", "G", "D", "S", "C"), 15, replace = T) %>% paste0(collapse = "")})
alleles <- sample(c(paste0("DRB1_010", 1:9)), size = 30, replace = T)
# # alleles <- unique(alleles)
# substr(peptides[1], 1, 1) <- "U"
# substr(peptides[5], 1, 1) <- "U"
# substr(peptides[10], 1, 1) <- "U"

### quick-setting the defaults for testing
# pepfile_loc = NULL
# paired_input = T
# type = c("score_el", "score_ba", "rank_el", "rank_ba", "aff")
# output_format = "wide"
# result_files_location = NULL
# threads = 1
# keep_pep = FALSE
# software_loc = NULL
# tmppep_loc = NULL
# logfile_loc = "log.txt"

##################### FUNCTIONS
RunNetMHCIIpan <- function(alleles,
                           peptides = NULL,
                           pepfile_loc = NULL,
                           paired_input = F,
                           type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"),
                           output_format = "long",
                           result_files_location = NULL,
                           threads = 1,
                           keep_pep = FALSE,
                           software_loc = NULL,
                           tmppep_loc = NULL,
                           logfile_loc = "log.txt") {
  
  # getting netMHCIIpan location from config file if not specified by the user
  if(is.null(software_loc)) {software_loc <- configs$software_loc_MHC_II}
  
  # getting temporary pepfile location from config file if not specified by the user
  if(is.null(tmppep_loc) & is.null(pepfile_loc)) {
    tmppep_loc <- configs$tmppep_loc_MHC_II
    tmppep_loc_for_paired <- configs$tmppep_loc_MHC_II
  }
  
  # setting temporary pepfile location for paired analysis
  tmppep_loc_for_paired <- configs$tmppep_loc_MHC_II
  
  
  # do not allow paired mode in case of pepfile input
  if(!is.null(pepfile_loc) & paired_input) {
    msg <- "Paired mode is not allowed when epitopes are loaded from pep file using argument pepfile_loc"
    LogIt(msg, "RunNetMHCIIpan", logfile_loc, "stop")
    stop(msg)
  }
  
  # if pepfile_loc is defined, then use it as a source of input peptides instead of creating a temporary pepfile
  if(!is.null(pepfile_loc)) {
    tmppep_loc <- pepfile_loc
    peptides <- readLines(pepfile_loc)
    keep_pep <- TRUE
    tmppep_loc_for_paired <- configs$tmppep_loc
    
    # checking if all peptides are OK in the pepfile, break if not
    invalid_peptides <- peptides[!IsValidPeptide(peptides, hla_type = 2)]
    if(length(invalid_peptides) > 0) {
      msg <- paste0("input contains too short peptides and/or sequences containing invalid characters. \nPlease remove the following peptides before starting the prediction:\n",
                    paste0(invalid_peptides, collapse = ", "))
      LogIt(msg, "RunNetMHCIIpan", logfile_loc, "stop")
      stop(msg)
    }
  }
  
  # filtering for valid peptides (only if peptides are loaded directly from R, not from pepfile, and if it is not a paired analysis)
  if(is.null(pepfile_loc) & !paired_input) {
    is_valid <- IsValidPeptide(peptides, hla_type = 2)
    invalid_peptides <- peptides[!is_valid]

    peptides <- peptides[is_valid]
    if(length(invalid_peptides) > 0) {
      msg <- paste("Skipped the following invalid peptides:", paste0(invalid_peptides, collapse = ", "))
      LogIt(msg, funcname = "RunNetMHCIIpan", logfile_loc, "message")
      message(msg)
    }
  }
  
  # filtering for valid alleles (only if it is not a paired analysis)
  software_loc_splt <- strsplit(software_loc, "\\/")[[1]]
  alleles_supported <- readLines(paste0(c(software_loc_splt[-length(software_loc_splt)], "data", "allele.list"), collapse = "/")) %>% strsplit(" ") %>% unlist()
  if(!paired_input) {
    is_valid_allele <- alleles %in% alleles_supported
    invalid_alleles <- alleles[!is_valid_allele]
    alleles <- alleles[is_valid_allele]
    if(length(invalid_alleles) > 0) {
      msg <- "Skipped the following unsupported alleles:"
      LogIt(paste(msg, paste0(invalid_alleles, collapse = ", ")), funcname = "RunNetMHCpan", logfile_loc, "message")
      message(msg)
    }
  }
  
  # if any of the bval_types is not supported, return a warning message
  bvaltypes_unsupported <- setdiff(type, c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"))
  if(all(type %in% bvaltypes_unsupported)) {
    msg <- "None of the binding value types are supported by netMHCIIpan 4.0."
    LogIt(msg, "RunNetMHCIIpan", logfile_loc, "stop")
    stop(msg)
  } else if (length(bvaltypes_unsupported) > 0) {
    msg <- paste0("The following binding value types are not supported by netMHCIIpan 4.0, so they have been skipped: ",
                  paste0(bvaltypes_unsupported, collapse = ", "))
    LogIt(msg, "RunNetMHCIIpan", logfile_loc, "warning")
    warning(msg)
  }
  
  # if number of threads is higher than CPUs available, reset number of threads
  n_cores <- detectCores()
  if(threads > n_cores | threads < 1) {
    msg <- paste0("Resetting number of threads to be used (", threads,") to the number of available CPUs (", n_cores, ")")
    LogIt(msg, "RunNetMHCIIpan", logfile_loc, "message")
    message(msg)
    threads <- n_cores
  }
  
  # if paired_input is TRUE treat alleles and peptides as pairs
  # forces the usage of "long" output format
  if(paired_input) {
    if(output_format == "wide") {
      msg <- "'wide' output format is not supported in paired mode, using 'long' output."
      LogIt(msg, "RunNetMHCIIpan", logfile_loc, "warning")
      warning(msg)
      output_format <- "long"
    }
    
    if(length(alleles) != length(peptides)) {
      msg <- "Pairwise prediction cannot be performed because length of allele vector is not equal to the length of peptide vector"
      LogIt(msg, "RunNetMHCIIpan", logfile_loc, "stop")
      stop(msg)
    }
    
    is_valid_allele <- alleles %in% alleles_supported
    invalid_alleles <- alleles[!is_valid_allele]
    alleles <- alleles[is_valid_allele]
    if(length(invalid_alleles) > 0) {
      msg <- paste("Skipped the following unsupported alleles:", paste0(invalid_alleles, collapse = ", "))
      LogIt(msg, funcname = "RunNetMHCpan", logfile_loc, "message")
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
      LogIt(msg, "RunNetMHCIIpan", logfile_loc, "message")
      message(msg)
    }
    
    # collecting peptides per alleles
    peptides_per_alleles <- aggregate(peptides ~ alleles, FUN =  function(x) {x})
    peptides_per_alleles <- set_names(peptides_per_alleles$peptides, peptides_per_alleles$alleles)
    
    
    if(threads == 1) {
      outlist <- imap(peptides_per_alleles, ~RunNetMHCIIpan(alleles = .y, peptides = .x,
                                                            type = type, output_format = output_format,
                                                            result_files_location = result_files_location,
                                                            threads = 1, keep_pep = keep_pep, software_loc = software_loc,
                                                            tmppep_loc = tmppep_loc_for_paired))
    } else {
      plan(multisession, workers = threads)
      outlist <- future_imap(peptides_per_alleles, ~RunNetMHCIIpan(alleles = .y, peptides = .x,
                                                                   type = type, output_format = output_format,
                                                                   result_files_location = result_files_location,
                                                                   threads = 1, keep_pep = keep_pep, software_loc = software_loc,
                                                                   tmppep_loc = paste0(tmppep_loc_for_paired, ".", gsub("\\:", "-", .y))))
    }
    
    return(set_rownames(do.call(rbind.data.frame, outlist), NULL))
  }
  
  # writing temporary pepfile to its location
  if(is.null(pepfile_loc)) {write(peptides, file = tmppep_loc)}
  
  # generating commands to run
  cmds <- sapply(alleles, function(allele) {
    paste0(software_loc, " -inptype 1 -f ", tmppep_loc  ," -BA -a ", allele, sep = "")
  }, USE.NAMES = F)
  
  # if a result file directory is specified
  if(!is.null(result_files_location)) {
    suppressWarnings(dir.create(result_files_location, recursive = TRUE))  # creating result files directory
    alleles_underscore <- gsub("\\:", "_", alleles) # windows cannot handle colons in filenames...
    alleles_underscore <- gsub("\\-", "_", alleles_underscore)
    if(is.null(type)) {
      cmds <- paste0(cmds, " > ", result_files_location, "/", alleles_underscore, ".txt")
    } else {
      cmds <- paste0(cmds, " | tee ", result_files_location, "/", alleles_underscore, ".txt")
    }
  }
  
  
  # RUNNING PREDICTIONS - and optionally collecting binding values
  if(!is.null(type)) {
    results <- RunCommand(cmds, threads = threads, intern = T)
    if(!keep_pep) {file.remove(tmppep_loc)} # removing temporary pepfile
    outobj <- InstantRecognMatrixII(results, type = type, output_format = output_format)
    
    # stroring version number
    predictor.version <- results[[1]][grepl("# NetMHCIIpan version", results[[1]])] %>% strsplit(" ") %>% unlist()
    attr(outobj, "predictor_version") <- predictor.version[length(predictor.version)]
    attr(outobj, "HLA_type") <- "II"
    
    return(outobj)
  } else {
    RunCommand(cmds, threads = threads, intern = F)
    if(!keep_pep) {file.remove(tmppep_loc)} # removing temporary pepfile
  }
  
  # printing final message
  if(is.null(type)) {
    msg <- paste0("Predictions are done, see result files in ", result_files_location)
    LogIt(msg, "RunNetMHCIIpan", logfile_loc, "message")
    message(msg)
  }
  
  
  # printing final message
  if(is.null(type)) {
    msg <- paste0("Predictions are done, see result files in ", result_files_location)
    LogIt(msg, "RunNetMHCIIpan", logfile_loc, "message")
    message(msg)
  }
}


##################### TESTING
# RunNetMHCIIpan(alleles, peptides, paired_input = T)
