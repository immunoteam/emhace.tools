#' @export

ReduceEpitopeSimilarity <- function(epitopes, similarity = 0.5, threads = 7,
                                   clustalo_path, using_windows = TRUE, keep_temp = F,
                                   tempdir = getwd(), save_histograms = T, 
                                   verbose = F) {
  strings <- as.vector(rbind(paste0(">", epitopes), epitopes))
  writeLines(strings, "epitopes.fasta")
  
  message("1 - Running ClustalO, creating the distrance matrix...")
  if(using_windows) clustalo_path <- paste0('"', clustalo_path, '"')
  invisible(system(paste0(clustalo_path, " -i epitopes.fasta -o tmp_out.fasta -v --distmat-out=dist.txt --threads=", 
                threads,  " --full --force"), show.output.on.console = F, wait = T))
  
  message("2 - ClustalO finished.")
  dist_matr <- data.table::fread(input = "dist.txt", skip = 1,  header = F)
  dist_matr <- subset(dist_matr, select = colnames(dist_matr)[2:ncol(dist_matr)])
  dist_matr <- as.matrix(dist_matr)
  dist_matr[dist_matr == 0] <- NA
  if(length(epitopes) != nrow(dist_matr)) {print("Something went wrong during ClustalO run!")}
  dimnames(dist_matr) <- list(epitopes, epitopes)
  gc()
  values_before <- as.vector(dist_matr)
  values_before <- values_before[!is.na(values_before)]
  num <- length(epitopes)
  
  message("3 - Reducing similarity... (this may take a while)")
  while (min(dist_matr, na.rm = TRUE) < similarity) {
    txtProgressBar(0, num, initial = -length(epitopes) + num, width = 2)
    wm <- which(dist_matr == min(dist_matr, na.rm = TRUE), arr.ind = TRUE)
    means <- colMeans(dist_matr[,wm[,2]], na.rm = T)
    index <- wm[which.min(means), 2]
    dist_matr <- dist_matr[-index, -index]
    epitopes <- epitopes[-index]
  }
  
  message("4 - Creating plots...")
  values_after <- as.vector(dist_matr)
  values_after <- values_after[!is.na(values_after)]
  df <- data.frame(values = c(values_before, values_after),
                  group = factor(c(rep("Before", length(values_before)), rep("After", length(values_after))), 
                                 levels = c("Before", "After")))
  g <- ggplot2::ggplot(df, ggplot2::aes(x = values)) +
    ggplot2::geom_histogram(alpha=0.5, position = "identity") +
    ggplot2::facet_wrap(~group) +
    ggplot2::labs(x="K-tuple distance", y = "Count", fill= "") + 
    ggplot2::scale_fill_manual(values=c("Blue", "Red")) +
    ggplot2::scale_y_continuous(trans = "sqrt") +
    ggplot2::geom_vline(xintercept = similarity, linetype = "dashed")
  print(g)
  return(epitopes)
}

# TEST
# out <- sapply(1:1000, function(x) sample(rownames(protr::AABLOSUM62), 9, replace = FALSE) %>% paste0(collapse = ""))
# ReduceEpitopeSimilarity(out, clustalo_path = "clustalo")
