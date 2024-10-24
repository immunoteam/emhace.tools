#' @export
RunNetMHCpan_ng <- function(alleles,
                                 peptides,
                                 paired_input = F,
                                 value_type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"),
                                 output_format = "long",
                                 threads = 1,
                                 version_number = "4.1",
                                 result_files_location = NULL,
                                 software_path = NULL,
                                 tmppep_loc = NULL) {
  
  # checking whether version number is correct
  if (!version_number %in% c("4.0", "4.1")) {
    stop("the version number of NetMHCpan is either not inputed as a string or is not supported. Make sure to use either NetMHCpan 4.0 or 4.1.")
  }
  
  # removing asterisks from allele names
  # asterisk-containing allele names are being misinterpreted by bash
  alleles <- gsub("\\*", "", alleles)
  
  # setting tmppep_loc if there was no pre-set
  software_path_splt <- strsplit(software_path, "\\/")[[1]]
  if(is.null(tmppep_loc)) {tmppep_loc <- paste0(paste0(software_path_splt[-length(software_path_splt)], collapse = "/"), "/tmp.pep")}
  
  # loading supported alleles
  alleles_supported <- unlist(strsplit(readLines(paste0(c(software_path_splt[-length(software_path_splt)],
                                                          "data", "allelenames"), collapse = "/")), " "))
  
  # performing input pre-filtering
  if(!paired_input) {
    # checking for peptide validity
    is_valid <- IsValidPeptide(peptides)
    invalid_peptides <- peptides[!is_valid]
    
    peptides <- peptides[is_valid]
    peptides <- unique(peptides)
    if(length(invalid_peptides) > 0) {
      msg <- paste("Skipped the following invalid peptides:", paste0(invalid_peptides, collapse = ", "))
      message(msg)
    }
    
    # checking for allele validity
    is_valid_allele <- alleles %in% alleles_supported
    invalid_alleles <- alleles[!is_valid_allele]
    alleles <- unique(alleles[is_valid_allele])
    if(length(invalid_alleles) > 0) {
      msg <- "Skipped the following unsupported alleles:"
      message(paste(msg, paste0(invalid_alleles, collapse = ", ")))
    }
  } else {
    # check if wide format was used instead of long
    if(output_format == "wide") {
      msg <- "'wide' output format is not supported in paired mode, using 'long' output."
      warning(msg)
      output_format <- "long"
    }
    
    # check whether the equal number of alleles and peptides are inputed
    if(length(alleles) != length(peptides)) {
      msg <- "Pairwise prediction cannot be performed because length of allele vector is not equal to the length of peptide vector"
      stop(msg)
    }
    
    # checking for allele-peptide pair validity together
    is_valid_peptide <- IsValidPeptide(peptides)
    is_valid_allele <- alleles %in% alleles_supported
    is_valid_pair <- is_valid_peptide & is_valid_allele
    invalid_pairs <- paste0(alleles, "-", peptides)[!is_valid_pair]
    
    if(length(invalid_pairs) > 0) {
      msg <- "Skipped the following invalid allele-peptide pairs:"
      message(paste(msg, paste0(invalid_pairs, collapse = ", ")))
    }
  }

  # filtering for supported output value types
  if(version_number == "4.1") {
    bvaltypes_unsupported <- setdiff(value_type, c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"))
  } else if (version_number == "4.0") {
    bvaltypes_unsupported <- setdiff(value_type, c("Score_BA", "Rank_BA", "Aff_nm", "Exp"))
  }
  
  if(!is.null(value_type) & all(value_type %in% bvaltypes_unsupported)) {
    msg <- "None of the binding value types are supported by the selected version of netMHCpan."
    stop(msg)
  } else if (length(bvaltypes_unsupported) > 0) {
    msg <- paste0("The following binding value types are not supported by the selected version of netMHCpan, they will be skipped: ",
                  paste0(bvaltypes_unsupported, collapse = ", "))
    value_type <- setdiff(value_type, bvaltypes_unsupported)
    warning(msg)
  }
  
  # if number of threads is higher than CPUs available, reset number of threads
  n_cores <- parallel::detectCores()
  if(threads > n_cores | threads < 1) {
    msg <- paste0("Resetting number of threads to be used (", threads,") to the number of available threads (", n_cores, ")")
    message(msg)
    threads <- n_cores
  }
  
  # creating the temporary input pepfiles
  if(paired_input) {
    # creating pepfiles
    peptides_per_threads <- data.frame(peptides, alleles)
    peptides_per_threads <- peptides_per_threads[order(peptides_per_threads$alleles), ]
    if(threads > 1) {
      peptides_per_threads$thread <- cut(seq_along(peptides), threads, labels = FALSE)
    } else {
      peptides_per_threads$thread <- 1
    }
    
    peptides_per_threads$pepfile_id <- paste0(peptides_per_threads$alleles, ".", peptides_per_threads$thread)
    pepfile_ids_unq <- unique(peptides_per_threads$pepfile_id)
    
    tmppep_paths <- purrr::map_chr(pepfile_ids_unq, function(.x) {
      crnt_path <- paste0(tmppep_loc, ".", .x)
      write(peptides_per_threads$peptides[peptides_per_threads$pepfile_id == .x], file = crnt_path)
      return(crnt_path)
    })
    
    # creating inputdf
    inputdf <- data.frame(tmppep_paths, alleles = purrr::map_chr(strsplit(tmppep_paths, "\\."), ~tail(.x, 2)[1]))
  } else if(!paired_input) {
    # creating pepfiles
    if(threads > 1) {thread_column <- cut(seq_along(peptides), threads, labels = FALSE)} else {thread_column <- 1}
    peptides_per_threads <- data.frame(peptides, thread = thread_column)
    tmppep_paths <- purrr::map_chr(1:threads, function(.x) {
      crnt_path <- paste0(tmppep_loc, ".", .x)
      write(peptides_per_threads$peptides[peptides_per_threads$thread == .x], file = crnt_path)
      return(crnt_path)
    })
    
    # creating inputdf
    inputdf <- expand.grid(tmppep_paths, alleles, stringsAsFactors = FALSE)
    colnames(inputdf) <- c("tmppep_paths", "alleles")
  }
  
  # adding commands and threads to inputdf
  inputdf$cmd <- paste0(software_path, " -inptype 1 -f ", inputdf[, 1]," -BA -a ", inputdf[, 2], sep = "")
  inputdf$thread <- purrr::map_chr(strsplit(inputdf[, 1], "\\."), ~tail(.x, 1))
  
  # modifying commands to include result files output location
  if(!is.null(result_files_location)) {
    suppressWarnings(dir.create(result_files_location, recursive = TRUE))  # creating result files directory
    inputdf$alleles_underscore <- gsub("\\:", "_", inputdf[, 2]) # windows cannot handle colons in filenames
    inputdf$alleles_underscore <- gsub("\\-", "_", inputdf$alleles_underscore)
    if(is.null(value_type)) {
      inputdf
      inputdf$cmd <- paste0(inputdf$cmd, " > ", result_files_location, "/", inputdf$alleles_underscore, ".", inputdf$thread, ".txt")
    } else {
      inputdf$cmd <- paste0(inputdf$cmd, " | tee ", result_files_location, "/", inputdf$alleles_underscore, ".", inputdf$thread, ".txt")
    }
  }
  
  # running predictions and collecting outputs (if asked)
  if(!is.null(value_type)) {
    results <- RunCommand(inputdf$cmd, threads = threads, intern = T)
    file.remove(unique(inputdf$tmppep_paths)) # removing temporary pepfiles
    outobj <- CollectBindingResults(results, value_type = value_type, output_format = output_format, version_number = version_number)
    outobj$allele <- gsub("\\*", "", outobj$allele)
    
    if(paired_input) outobj <- outobj[match(paste0(alleles, ".", peptides), paste0(outobj$allele, ".", outobj$peptide)), ]
    
    # stroring version number
    predictor.version <- unlist(strsplit(results[[1]][grepl("# NetMHCpan version", results[[1]])], " "))
    attr(outobj, "predictor_version") <- predictor.version[length(predictor.version)]
    attr(outobj, "HLA_type") <- "I"
    
    return(outobj)
  } else {
    RunCommand(cmds, threads = threads, intern = F)
    file.remove(inputdf$tmppep_paths) # removing temporary pepfiles
  }
}