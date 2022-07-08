#### THE FOLLOWING FILE CONTAINS FUNCTIONS FROM THE OLD RecognMatrixTools PACKAGE
#### DEVELOPED FOR THE WS
#### USE ONLY AS A SOURCE OF IDEAS

require(magrittr)
require(reshape2)


RunNetMHCpan <- function(alleles_loc, pep_loc, results_out_loc,
                         parallel_threads = FALSE, send_notification_to = FALSE, netmhcipan_loc = "~/Programok/netMHCpan-4.0/netMHCpan") {
  alleles = readLines(alleles_loc)
  strings = NULL
  
  dir.create(results_out_loc, recursive = TRUE)
  
  for (j in 1:(length(alleles))) {
    strings = c(strings, paste0(netmhcipan_loc, " -inptype 1 -f ", pep_loc  ," -BA -a ", alleles[j], " > ", results_out_loc, alleles[j], ".txt", sep = ""))
  }
  
  if (parallel_threads == FALSE) {
    for (i in 1:length(strings)) {
      system(strings[i])
    }
  } else {
    library(doParallel)
    library(foreach)
    c3 <- makeCluster(20)
    registerDoParallel(c3)
    
    foreach (i = 1:length(strings)) %dopar% {
      system(strings[i])
    }
    stopCluster(c3)
  }
  
  if (send_notification_to != FALSE) {NotifyByMail(s = "netMHCpan process finished", b = "Yeah, it really finished.", r = send_notification_to)}
}

RunNetMHCIIpan <- function(alleles_loc, pep_loc, results_out_loc,
                           parallel_threads = FALSE, send_notification_to = FALSE, tmploc = "tmp/",
                           netmhciipan_loc = "~/Programok/netMHCIIpan-3.2/netMHCIIpan") {
  alleles = readLines(alleles_loc)
  strings = NULL
  
  dir.create(results_out_loc, recursive = TRUE)
  dir.create(tmploc, recursive = TRUE)
  
  for (j in 1:(length(alleles))) {
    strings = c(strings, paste0(netmhciipan_loc, " -inptype 1 -f ", pep_loc  ," -a ", alleles[j], " -t ", paste0(tmploc, "/", alleles[j]),
                                " >", results_out_loc, alleles[j], ".txt", sep = ""))
  }
  
  if (parallel_threads == FALSE) {
    for (i in 1:length(strings)) {
      system(strings[i])
    }
  } else {
    library(doParallel)
    library(foreach)
    c3 <- makeCluster(parallel_threads)
    registerDoParallel(c3)
    
    foreach (i = 1:length(strings)) %dopar% {
      system(strings[i])
    }
    stopCluster(c3)
  }
  
  if (send_notification_to != FALSE) {NotifyByMail(s = "netMHCIIpan process finished", b = "Yeah, it really finished.", r = send_notification_to)}
}

CreateRecognMatricesXmers <- function(result_files_maindir, pep_loc, rp_or_aff = "both") {
  peps <- list.files(pep_loc)
  mers <- gsub("\\.pep", "", peps)
  dirs <- list.dirs(result_files_maindir, recursive = F)
  recogn_matrices_per_mers <- list()
  for (i in 1:length(mers)) {
    epitopes_loc <- paste0(pep_loc, "/", peps[i])
    result_files_dir <- dirs[i]
    CreateRecognMatrix(result_files_dir = result_files_dir,
                       epitopes_loc = epitopes_loc,
                       rp_or_aff = rp_or_aff) -> recogn_matrices_per_mers[i]
  }
  return(AssembleRecognMatrices(recogn_matrices_per_mers)) # MÉG HIBÁS!!!!!!!!!!!!!!! AZ AFFOKAT ÉS RPKET NEM KEZELI KÜLÖN!!!
}


CreateRecognMatrix <- function(result_files_dir, rp_or_aff = "both") {
  alleles <- list.files(result_files_dir) %>% gsub("\\.txt", "", .)
  resfiles <- list.files(result_files_dir, full.names = T)
  results <- pblapply(resfiles, function(x) {LoadnetMHCpanResfile(x)}) %>% set_names(alleles) %>% do.call(rbind.data.frame, .)
  results$HLA <- results$HLA %>% gsub("\\*", "", .)
  
  if(rp_or_aff == "aff") {
    recogn_matrices <- list(aff = acast(results, HLA ~ peptide, value.var = "aff"))
  } else if(rp_or_aff == "rp") {
    recogn_matrices <- list(rp = acast(results, HLA ~ peptide, value.var = "rp"))
  } else {
    if(rp_or_aff != "both") {print("Parameter value is not defined correctly: creating both matrices")}
    recogn_matrices <- list(aff = acast(results, HLA ~ peptide, value.var = "aff"),
                            rp = acast(results, HLA ~ peptide, value.var = "rp"))
  }
  return(recogn_matrices)
}

CreateRecognMatrixOld <- function(result_files_dir, epitopes_loc, rp_or_aff, fastmode = TRUE) {
  recogn_matrix <- CreateEmptyRecognMtx(result_files_dir, epitopes_loc) # creates a template to be filled in the next step
  if(fastmode == FALSE) {
    FillRecognMtx(empty_recogn_mtx = recogn_matrix,
                  result_files_dir = result_files_dir,
                  rp_or_aff = rp_or_aff) -> recogn_matrices
  } else {
    FillRecognMtxFast(empty_recogn_mtx = recogn_matrix,
                      result_files_dir = result_files_dir,
                      rp_or_aff = rp_or_aff) -> recogn_matrices
  }
  return(recogn_matrices)
}

LoadnetMHCpanResfile <- function(resfile_path, simple = TRUE) {
  resfile <- readLines(resfile_path)
  resfile <- resfile[51:(length(resfile)-5)]
  strsplit(resfile, " ") %>% lapply(function(x) {x[x != ""]}) %>% lapply(function(x) {x[1:14]}) %>% do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>% set_colnames(c("pos", "HLA", "peptide", "Core", "Of", "Gp", "Gl", "Ip", "Il", "Icore", "Identity", "score", "aff", "rp")) -> resdf
  resdf$aff <- as.numeric(resdf$aff)
  resdf$rp <- as.numeric(resdf$rp)
  if (simple) {
    return(resdf[, c("HLA", "peptide", "aff", "rp")])
  } else {
    resdf$pos %<>% as.numeric()
    resdf$Of %<>% as.numeric()
    resdf$Gp %<>% as.numeric()
    resdf$Gl %<>% as.numeric()
    resdf$Ip %<>% as.numeric()
    resdf$Il %<>% as.numeric()
    resdf$score %<>% as.numeric()
    return(resdf)
  }
}

CreateRecognMatrixNetMHCIIpan <- function(result_files_dir, epitopes_loc, rp_or_aff) {
  recogn_matrix <- CreateEmptyRecognMtx(result_files_dir, epitopes_loc) # creates a template to be filled in the next step
  FillRecognMtxNetMHCIIpan(empty_recogn_mtx = recogn_matrix,
                           result_files_dir = result_files_dir,
                           rp_or_aff = rp_or_aff) -> recogn_matrices
  return(recogn_matrices)
}

CreateEmptyRecognMtx <- function(result_files_dir, epitopes_loc) {
  alleles <- list.files(result_files_dir)
  alleles <- gsub("\\.txt", "", alleles)
  epitopes <- read.delim(epitopes_loc, stringsAsFactors = FALSE, header = FALSE)[, 1]
  recogn_matrix <- matrix(nrow = length(alleles), ncol = length(epitopes),
                          dimnames = list(alleles, epitopes))
  return(recogn_matrix)
}

FillRecognMtx <- function(empty_recogn_mtx, result_files_dir, rp_or_aff) {
  alleles <- rownames(empty_recogn_mtx)
  epitopes <- colnames(empty_recogn_mtx)
  recogn_matrix_aff <- empty_recogn_mtx
  recogn_matrix_rp <- empty_recogn_mtx
  files <- list.files(result_files_dir, full.names = TRUE)
  for (i in 1:length(alleles)) { # borrowed from Mate
    print(i)
    data = readLines(paste(files[i], sep = ""))
    data = subset(data, grepl("HLA-", data) == TRUE & grepl("PEPLIST", data) == TRUE)
    if (length(data) == 0) {next} # if result file doesn't contain any results, skip it
    data = data[1:(length(data)-1)]
    data_processed = matrix(NA, ncol = 2, nrow = length(data))
    colnames(data_processed) = c("aff", "rank")
    rownames(data_processed) = epitopes
    for(j in 1:length(data)) {
      a = unlist(strsplit(data[j], " "))
      a = a[nchar(a) > 0]
      data_processed[j,] = c(aff = a[13], rank = a[14])
    }
    recogn_matrix_aff[i,] = data_processed[,1]
    recogn_matrix_rp[i, ] = data_processed[,2]
  }
  if (rp_or_aff == "aff") {class(recogn_matrix_aff) <- "numeric"; return(list(aff = recogn_matrix_aff))}
  if (rp_or_aff == "rp") {class(recogn_matrix_rp) <- "numeric"; return(list(aff = recogn_matrix_rp))}
  if (rp_or_aff == "both") {
    class(recogn_matrix_rp) <- "numeric"
    class(recogn_matrix_aff) <- "numeric"
    return(list(aff = recogn_matrix_aff, rp = recogn_matrix_rp))
  }
}

FillRecognMtxFast <- function(empty_recogn_mtx, result_files_dir, rp_or_aff) { # same as FillRecognMtx, but instead of splitting rows, it uses substr indecis to find the aff and rp values. Probably faster.
  alleles <- rownames(empty_recogn_mtx)
  epitopes <- colnames(empty_recogn_mtx)
  recogn_matrix_aff <- empty_recogn_mtx
  recogn_matrix_rp <- empty_recogn_mtx
  files <- list.files(result_files_dir, full.names = TRUE)
  for (i in 1:length(alleles)) { # borrowed from Mate
    data = unlist(data.table::fread(files[i], sep = "?", skip = 50, header = FALSE))
    data = data[1:(length(data)-5)]
    if (length(data) == 0) {next} # if result file doesn't contain any results, skip it
    data_processed = matrix(NA, ncol = 2, nrow = length(data))
    colnames(data_processed) = c("aff", "rank")
    rownames(data_processed) = epitopes
    begin_3 <- Sys.time()
    for (j in 1:length(data)) {
      a = data[j]
      data_processed[j, ] = c(aff = substr(a, 97, 103), rp = substr(a, 105, 111))
    }
    end_3 <- Sys.time()
    
    recogn_matrix_aff[i,] = data_processed[,1]
    recogn_matrix_rp[i, ] = data_processed[,2]
  }
  if (rp_or_aff == "aff") {class(recogn_matrix_aff) <- "numeric"; return(list(aff = recogn_matrix_aff))}
  if (rp_or_aff == "rp") {class(recogn_matrix_rp) <- "numeric"; return(list(aff = recogn_matrix_rp))}
  if (rp_or_aff == "both") {
    class(recogn_matrix_rp) <- "numeric"
    class(recogn_matrix_aff) <- "numeric"
    return(list(aff = recogn_matrix_aff, rp = recogn_matrix_rp))
  }
}

FillRecognMtxNetMHCIIpan <- function(empty_recogn_mtx, result_files_dir, rp_or_aff) {
  alleles <- rownames(empty_recogn_mtx)
  epitopes <- colnames(empty_recogn_mtx)
  recogn_matrix_aff <- empty_recogn_mtx
  recogn_matrix_rp <- empty_recogn_mtx
  files <- list.files(result_files_dir, full.names = TRUE)
  for (i in 1:length(alleles)) { # borrowed from Mate
    data = readLines(files[i])
    if (length(data) <= 1) {next} # if result file doesn't contain any results, skip it
    data = data[12:(length(data)-3)]
    
    strsplit(data, split = " ") %>%
      lapply(function(x) {x[x != ""]}) %>%
      sapply(function(x) {c(aff = x[9], rp = x[10]) %>% as.numeric()}) %>%
      t() -> data_processed
    
    recogn_matrix_aff[i, ] = data_processed[,1]
    recogn_matrix_rp[i, ] = data_processed[,2]
  }
  if (rp_or_aff == "aff") {class(recogn_matrix_aff) <- "numeric"; return(list(aff = recogn_matrix_aff))}
  if (rp_or_aff == "rp") {class(recogn_matrix_rp) <- "numeric"; return(list(aff = recogn_matrix_rp))}
  if (rp_or_aff == "both") {
    class(recogn_matrix_rp) <- "numeric"
    class(recogn_matrix_aff) <- "numeric"
    return(list(aff = recogn_matrix_aff, rp = recogn_matrix_rp))
  }
}

AssembleRecognMatrices <- function(recogn_matrices_per_mers) {
  return(do.call("cbind", recogn_matrices_per_mers))
}

ReduceEpitopeSimilarity <- function(epitopes, similarity) {
  library(seqinr)
  epitope_names <- paste0(" epitope_", 1:length(epitopes))
  epitopes_list <- lapply(epitopes, FUN = s2c)
  write.fasta(epitopes_list, names = epitope_names, file.out = "tmp.fasta")
  print("Running ClustalO, creating the distmat")
  system("clustalo -i tmp.fasta -o tmp_out.fasta -v --distmat-out=dist.txt --threads=2 --full --force")
  
  dist_matr <- data.table::fread(input = "dist.txt", skip = 1,  header = F)
  dist_matr <- subset(dist_matr, select = colnames(dist_matr)[2:ncol(dist_matr)])
  dist_matr <- as.matrix(dist_matr)
  
  dist_matr[dist_matr == 0] = NA
  
  minimum = min(dist_matr, na.rm = TRUE)
  if(length(epitopes) != nrow(dist_matr)) {print("Something went wrong during clustalo run!")}
  rownames(dist_matr) = epitopes
  colnames(dist_matr) = epitopes
  gc()
  begin <- Sys.time()
  print(begin)
  while (min(dist_matr, na.rm = TRUE) < similarity) {
    print(paste(min(dist_matr, na.rm = TRUE), length(epitopes)))
    wm = which(dist_matr == min(dist_matr, na.rm = TRUE), arr.ind = TRUE)
    means = colMeans(dist_matr, na.rm = TRUE)
    means = means[wm[,1]]
    index = wm[which.min(means), 1]
    dist_matr = dist_matr[-index, -index]
    epitopes = epitopes[-index]
    gc()
  }
  end <- Sys.time()
  print(end-begin)
  return(epitopes)
}

MHCFlurrySupportedAlleles <- function() { # returns the names of alleles, which are supported by MHCFlurry
  return(system("mhcflurry-predict --list-supported-alleles", intern = TRUE))
}

RunMHCFlurry <- function(alleles = NULL, peptides = NULL, outfile = "out.csv",
                         parallel_threads = 1, tmpfile_name = "tmp_mhcflurry_input.csv") { # Run MHCFlurry prediction. Input alleles in netmhcpan format!
  
  input_df <- expand.grid(alleles, peptides, stringsAsFactors = FALSE)
  colnames(input_df) <- c("allele", "peptide")
  
  if (parallel_threads == 1) {
    write.table(input_df, file = tmpfile_name, sep = ",", quote = FALSE, row.names = FALSE)
    cmd <- paste0("mhcflurry-predict ",  tmpfile_name, " --out ", outfile)
    system(cmd)
    file.remove(tmpfile_name)
  } else { # running on single thread
    nr <- nrow(input_df)
    input_df_list <- split(input_df, f = rep(1:4, each = ceiling(nr/4)))
    
    tmpfile_names <- NULL
    outfile_names <- NULL
    for (i in 1:length(input_df_list)) {
      tmpfile_names[i] <- paste(tmpfile_name, i, sep = "_")
      outfile_names[i] <- paste(outfile, i, sep = "_")
      write.table(input_df_list[[i]], file = tmpfile_names[i], sep = ",",
                  quote = FALSE, row.names = FALSE)
    }
    
    # running the predictions parallelly
    cmd <- paste0("mhcflurry-predict ",  tmpfile_names, " --out ", outfile_names)
    library(doParallel)
    library(foreach)
    c3 <- makeCluster(parallel_threads)
    registerDoParallel(c3)
    
    foreach (i = 1:length(cmd)) %dopar% {
      system(cmd[i])
    }
    stopCluster(c3)
    file.remove(tmpfile_names)
    
    # collecting the results into one file
    results_list <- list()
    for (i in 1:length(outfile_names)) {
      results_list[[i]] <- read.delim(outfile_names[i], stringsAsFactors = FALSE, sep = ",")
    }
    results <- do.call(rbind, results_list)
    write.table(results, file = outfile, sep = ",", quote = FALSE, row.names = FALSE)
    file.remove(outfile_names)
  }
}

CreateRecognMatrixMHCFlurry <- function(result_file, rp_or_aff) {
  resdf <- read.delim(result_file, stringsAsFactors = FALSE, sep = ",")
  alleles <- unique(resdf$allele)
  peptides <- unique(resdf$peptide)
  
  if (grepl("rp", rp_or_aff)) {
    recogn_matrix_rp <- acast(resdf[, c("allele", "peptide", "mhcflurry_prediction_percentile")],
                              allele ~ peptide, value.var = "mhcflurry_prediction_percentile")
    return(recogn_matrix_rp)
  }
  
  if (grepl("aff", rp_or_aff)) {
    recogn_matrix_aff <- acast(resdf[, c("allele", "peptide", "mhcflurry_prediction")],
                               allele ~ peptide, value.var = "mhcflurry_prediction")
    return(recogn_matrix_aff)
  }
  
  if (grepl("both", rp_or_aff)) {
    recogn_matrix_rp <- acast(resdf[, c("allele", "peptide", "mhcflurry_prediction_percentile")],
                              allele ~ peptide, value.var = "mhcflurry_prediction_percentile")
    recogn_matrix_aff <- acast(resdf[, c("allele", "peptide", "mhcflurry_prediction")],
                               allele ~ peptide, value.var = "mhcflurry_prediction")
    return(list(aff = recogn_matrix_aff, rp = recogn_matrix_rp))
  }
}

RunNetMHCstabpan<- function(alleles_loc, pep_loc, result_files_dir, parallel_threads = FALSE, send_notification_to = FALSE, netmhcstabpan_loc = "~/Programok/netMHCstabpan-1.0/netMHCstabpan") {
  dir.create(result_files_dir, recursive = TRUE)
  
  alleles = readLines(alleles_loc)
  
  peps = readLines(pep_loc)
  peps8 = peps[nchar(peps) == 8]
  write.table(peps8, file = paste0(result_files_dir, "/temppeps8"), quote = F, row.names = F, col.names = F)
  peps9 = peps[nchar(peps) == 9]
  write.table(peps9, file = paste0(result_files_dir, "/temppeps9"), quote = F, row.names = F, col.names = F)
  peps10 = peps[nchar(peps) == 10]
  write.table(peps10, file = paste0(result_files_dir, "/temppeps10"), quote = F, row.names = F, col.names = F)
  peps11 = peps[nchar(peps) == 11]
  write.table(peps11, file = paste0(result_files_dir, "/temppeps11"), quote = F, row.names = F, col.names = F)
  peps12 = peps[nchar(peps) == 12]
  write.table(peps12, file = paste0(result_files_dir, "/temppeps12"), quote = F, row.names = F, col.names = F)
  peps13 = peps[nchar(peps) == 13]
  write.table(peps13, file = paste0(result_files_dir, "/temppeps13"), quote = F, row.names = F, col.names = F)
  peps14 = peps[nchar(peps) == 14]
  write.table(peps14, file = paste0(result_files_dir, "/temppeps14"), quote = F, row.names = F, col.names = F)
  
  strings = NULL
  for (j in 1:(length(alleles))) {
    strings = c(strings,
                paste0(netmhcstabpan_loc, " -p ", result_files_dir, "/temppeps8 -a ", alleles[j], " > ", result_files_dir, "/", alleles[j], "_8.txt", sep = ""),
                paste0(netmhcstabpan_loc, " -p ", result_files_dir, "/temppeps9 -a ", alleles[j], " > ", result_files_dir, "/", alleles[j], "_9.txt", sep = ""),
                paste0(netmhcstabpan_loc, " -p ", result_files_dir, "/temppeps10 -a ", alleles[j], " > ", result_files_dir, "/", alleles[j], "_10.txt", sep = ""),
                paste0(netmhcstabpan_loc, " -p ", result_files_dir, "/temppeps11 -a ", alleles[j], " > ", result_files_dir, "/", alleles[j], "_11.txt", sep = ""),
                paste0(netmhcstabpan_loc, " -p ", result_files_dir, "/temppeps12 -a ", alleles[j], " > ", result_files_dir, "/", alleles[j], "_12.txt", sep = ""),
                paste0(netmhcstabpan_loc, " -p ", result_files_dir, "/temppeps13 -a ", alleles[j], " > ", result_files_dir, "/", alleles[j], "_13.txt", sep = ""),
                paste0(netmhcstabpan_loc, " -p ", result_files_dir, "/temppeps14 -a ", alleles[j], " > ", result_files_dir, "/", alleles[j], "_14.txt", sep = ""))
  }
  
  if (parallel_threads == FALSE) {
    for (i in 1:length(strings)) {
      system(strings[i])
    }
  } else {
    library(doParallel)
    library(foreach)
    c3 <- makeCluster(parallel_threads)
    registerDoParallel(c3)
    
    foreach (i = 1:length(strings)) %dopar% {
      system(strings[i])
    }
    stopCluster(c3)
  }
  
  #if (send_notification_to != FALSE) {NotifyByMail(s = "netMHCpan process finished", b = "Yeah, it really finished.", r = send_notification_to)}
  invisible(file.remove(paste0(result_files_dir, "/temppeps8"),
                        paste0(result_files_dir, "/temppeps9"),
                        paste0(result_files_dir, "/temppeps10"),
                        paste0(result_files_dir, "/temppeps11"),
                        paste0(result_files_dir, "/temppeps12"),
                        paste0(result_files_dir, "/temppeps13"),
                        paste0(result_files_dir, "/temppeps14")))
}

LoadStabpanResfile <- function(resfile_path) {
  resfile <- readLines(resfile_path)
  resfile <- resfile[50:(length(resfile)-5)]
  strsplit(resfile, " ") %>% lapply(function(x) {x[x != ""]}) %>% lapply(function(x) {x[1:7]}) %>% do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>% set_colnames(c("pos", "HLA", "peptide", "Identity", "Pred", "Thalf", "Rank_Stab")) -> resdf
  resdf$pos %<>% as.numeric()
  resdf$Pred %<>% as.numeric()
  resdf$Thalf %<>% as.numeric()
  resdf$Rank_Stab %<>% as.numeric()
  return(resdf)
}


CreateStabRecognMatrix <- function(pep_loc, result_files_dir, thalf_or_rank) {
  alleles = sort(unique(gsub("_.*", "", list.files(result_files_dir))))
  peps = readLines(pep_loc)
  recogn_matrix <- matrix(nrow = length(alleles), ncol = length(peps),
                          dimnames = list(alleles, peps))
  recogn_matrix_thalf <- recogn_matrix
  recogn_matrix_rank <- recogn_matrix
  files <- list.files(result_files_dir, full.names = TRUE)
  for (i in 1:length(alleles)) { # borrowed from Mate
    #print(i)
    dataofallele = matrix(ncol = 2, nrow = 0)
    colnames(dataofallele) = c("thalf", "rank")
    for(j in 1:length(grep(alleles[i], files, value = T))) {
      data = readLines(grep(alleles[i], files, value = T)[j])
      data = subset(data, grepl("HLA-", data) == TRUE & grepl("PEPLIST", data) == TRUE)
      if (length(data) == 0) {next} # if result file doesn't contain any results, skip it
      data = data[1:(length(data)-1)]
      data_processed = matrix(NA, ncol = 2, nrow = length(data))
      colnames(data_processed) = c("thalf", "rank")
      rn = NULL
      for(k in 1:length(data)) {
        a = unlist(strsplit(data[k], " "), use.names = F)
        a = a[nchar(a) > 0]
        rn = c(rn, a[3])
        data_processed[k,] = c(thalf = a[6], rank = a[7])
      }
      rownames(data_processed) = rn
      recogn_matrix_thalf[i,rn] = data_processed[,1]
      recogn_matrix_rank[i,rn] = data_processed[,2]
    }
  }
  class(recogn_matrix_thalf) = "numeric"
  class(recogn_matrix_rank) = "numeric"
  if (thalf_or_rank == "thalf") {
    return(list(thalf = recogn_matrix_thalf))}
  if (thalf_or_rank == "rp") {
    return(list(rank = recogn_matrix_rank))}
  if (thalf_or_rank == "both") {
    return(list(thalf = recogn_matrix_thalf, rank = recogn_matrix_rank))
  }
}
