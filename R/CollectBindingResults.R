#' @export
CollectBindingResults <- function(results, value_type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"), output_format = "long", hla_type = 1, version_number = "4.1"){
  if(hla_type == 1) {
    datarow_identifier <- "PEPLIST "
    to_be_trimmed_wb <- " <= WB"
    to_be_trimmed_sb <- " <= SB"
    if(version_number == "4.1") {
      coords <- c(2, 3, 4, 12, 13, 14, 15, 16)
      clnms_output <- c("allele", "peptide", "core", "Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm")
    } else if(version_number == "4.0") {
      coords <- c(2, 3, 4, 12, 13, 14)
      clnms_output <- c("allele", "peptide", "core", "Score_BA", "Aff_nm", "Rank_BA")
    } else {
      stop("Version number is not supported")
    }
    
  } else if (hla_type == 2) {
    datarow_identifier <- "Sequence "
    to_be_trimmed_wb <- " <=WB"
    to_be_trimmed_sb <- " <=SB"
    coords <- c(2, 3, 5, 8, 9, 11, 13, 12)
  } else {
    stop("The value of hla_type can be either 1 or 2!")
  }
  readin <- results[[1]][grep(datarow_identifier, results[[1]])]
  peplist <- sapply(readin, function(x) {unlist(strsplit(x, "\\s+"))[4]}, USE.NAMES = FALSE)
  
  # create long format
  resultsdf <- as.data.frame(do.call(rbind, lapply(results, function(x){
    x <- stringr::str_trim(x, "left")
    x <- x[grep(datarow_identifier, x)]
    x <- gsub(to_be_trimmed_wb, "", x)
    x <- gsub(to_be_trimmed_sb, "", x)
    mat <- do.call(rbind, strsplit(x, "\\s+"))
    mat[, coords]
  })), stringsAsFactors = F)
  for (i in 4:ncol(resultsdf)) {resultsdf[, i] <- as.numeric(resultsdf[, i])}
  
  colnames(resultsdf) <- clnms_output
  resultsdf <- cbind.data.frame(resultsdf[, 1:3], resultsdf[, colnames(resultsdf) %in% value_type, drop = FALSE])
  
  # create wide format
  if(output_format == "wide"){
    resultsdf <- lapply(value_type, function(x) reshape2::acast(resultsdf, allele ~ peptide, value.var = x, fun.aggregate = sum))
    names(resultsdf) <- value_type
  }
  
  return(resultsdf)
}
