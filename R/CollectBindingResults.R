CollectBindingResults <- function(results, value_type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"), output_format = "long"){
  readin <- results[[1]][grep("PEPLIST ",results[[1]])] 
  peplist <- sapply(readin, function(x) {unlist(strsplit(x, "\\s+"))[4]}, USE.NAMES = FALSE)
  
  # create long format
  resultsdf <- as.data.frame(do.call(rbind, lapply(results, function(x){
    x <- x[grep("PEPLIST ", x)]
    x <- gsub(" <= WB", "", x)
    x <- gsub(" <= SB", "", x)
    mat <- do.call(rbind, strsplit(x, "\\s+"))
    mat[, c(3,4,13,14,15,16,17)]
  })), stringsAsFactors = F)
  for (i in 3:7) {resultsdf[, i] <- as.numeric(resultsdf[, i])}
  
  colnames(resultsdf) <- c("allele", "peptide", "Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm")
  resultsdf <- cbind.data.frame(resultsdf[, 1:2], resultsdf[, colnames(resultsdf) %in% value_type, drop = FALSE])
  
  # create wide format
  if(output_format == "wide"){
    resultsdf <- lapply(value_type, function(x) reshape2::acast(resultsdf, allele ~ peptide, value.var = x, fun.aggregate = sum))
    names(resultsdf) <- value_type
  }
  
  return(resultsdf)
}
