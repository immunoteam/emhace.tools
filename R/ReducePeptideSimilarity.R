# TO BE DONE!


# setwd("C:/Munka/netmhcpan_tools/netMHCpan-tools/ReduceEpitopeSimilarity_playground/")
# load("../../../tumor_promiscuity/_objects/mhc_ligands_for_kl")
# epitopes = unlist(mhc_ligand[1:500,1])
# epitopes = unlist(mhc_ligand$Description[mhc_ligand$`Allele Name` == "HLA-A*02:01"][1:1000])

#B verzió: 2 csoport denzitását két külön ábrán

ReduceEpitopeSimilarity = function(epitopes, input_type, similarity = 0.5, threads = 7, clustalo_dir, keep_temp = F, tempdir = getwd(), save_histograms = T, 
                                   verbose = F) {
  strings = as.vector(rbind(paste0(">", epitopes), epitopes))
  writeLines(strings, "epitopes.fasta")
  
  cat("1 - Running ClustalO, creating the distmat")
  invisible(system(paste0("\"C:\\Program Files\\Clustal Omega\\clustalo.exe\" -i epitopes.fasta -o tmp_out.fasta -v --distmat-out=dist.txt --threads=", 
                threads,  " --full --force"), show.output.on.console = F, wait = T))
  cat("\n2 - Clustalo finished")
  dist_matr = fread(input = "dist.txt", skip = 1,  header = F)
  dist_matr = subset(dist_matr, select = colnames(dist_matr)[2:ncol(dist_matr)])
  dist_matr = as.matrix(dist_matr)
  dist_matr[dist_matr == 0] = NA
  
  if(length(epitopes) != nrow(dist_matr)) {print("Something went wrong during clustalo run!")}
  rownames(dist_matr) = epitopes
  colnames(dist_matr) = epitopes
  gc()
  values_before = as.vector(dist_matr)
  values_before = values_before[!is.na(values_before)]
  num = length(epitopes)
  cat("\n3 - Reducing similarity (this may take a while)\n")
  while (min(dist_matr, na.rm = TRUE) < similarity) {
    txtProgressBar(0, num, initial = -length(epitopes) + num, width = 2)
    wm = which(dist_matr == min(dist_matr, na.rm = TRUE), arr.ind = TRUE)
    means = colMeans(dist_matr[,wm[,2]], na.rm = T)
    index = wm[which.min(means), 2]
    dist_matr = dist_matr[-index, -index]
    epitopes = epitopes[-index]
  }
  cat("\n4 - Creating plots")
  values_after = as.vector(dist_matr)
  values_after = values_after[!is.na(values_after)]
  df = data.frame(values = c(values_before, values_after), group = factor(c(rep("Before", length(values_before)), rep("After", length(values_after))), 
                                                                          levels = c("Before", "After")))
  g = ggplot(df, aes(x = values)) + geom_histogram(alpha=0.5, position = "identity") + facet_wrap(~group) + labs(x="K-tuple distance", y = "Count", fill= "")  +  scale_fill_manual(values=c("Blue", "Red")) +
    scale_y_continuous(trans = "sqrt") + geom_vline(xintercept = similarity, linetype = "dashed")
  print(g)
  return(epitopes)
}

# epitopes_nonsim = ReduceEpitopeSimilarity(epitopes, similarity = 0.5)
