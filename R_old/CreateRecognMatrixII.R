#CreateRecognMatrix

library(magrittr)
library(reshape2)

#### TESTING ----------
# dir = "RunNetMHCpan_playground/TESLA_results/"
# # type = aff, rp, both
# type = c("Score_EL", "Rank_EL")
# write_to_csv = "createRecognMatrix_playground/results.csv"
# output_format = "wide"
# CreateRecognMatrix(dir = "RunNetMHCpan_playground/TESLA_results/", type = c("Score_EL", "Score_BA"), output_format = "wide", write_to_csv = "createRecognMatrix_playground/valami.txt")

#### FUNCTIONS ----------
#elutionScore-t beletenni
CreateRecognMatrixII = function(dir, type = c("Score_EL", "Rank_EL"), output_format = "long", write_to_csv = NULL) {
  result_files = list.files(dir)
  if (type[1] == "all") {type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm")}
  readin =  readLines(paste0(dir[1], result_files[1]))
  readin = readin[grepl("Sequence", readin)]
  peplist = unlist(lapply(readin, function(x) {unlist(strsplit(x, "\\s+"))[4]}))
  
  if(length(result_files) > 0) { 
    
    # create long format
    mtx = as.data.frame(do.call(rbind, lapply(result_files, function(x){
      temp = readLines(paste0(dir, x))
      temp = temp[grepl("Sequence ", temp)]
      temp = gsub(" <= SB| <= WB", "", temp)
      #this gives a warning massage 
      mat =  do.call(rbind, strsplit(temp, "\\s+"))
      mat[, c(3,4,9,10,12,13,14)]
    })), stringsAsFactors = F)
    
    for(i in 3:7) {mtx[,i] = as.numeric(mtx[,i])}
    colnames(mtx) = c("allele", "peptide", "Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm")
    
    mtx = cbind.data.frame(mtx[,1:2], mtx[, colnames(mtx) %in% type, drop = FALSE])
    
    # create wide format
    if (output_format == "wide"){
      mtx = lapply(type, function(x) acast(mtx, allele ~ peptide, value.var = x)) 
      names(mtx) = type
    }
    
    # return output
    if (!is.null(write_to_csv)) {
      if(output_format == "long") {
        write.csv(mtx, file = write_to_csv, row.names = F, quote = F)
      } else {
        lapply(names(mtx), function(x){
          write.csv(mtx[[x]], file = paste0(gsub(".csv", "", write_to_csv), "_", x, ".csv"), row.names = T, quote = FALSE)
        })
      }
    } else {mtx}
  } else "No file in directory"
}


#test
#letesztelni gsub v invisible() paranccsal /supress.warning() gyorsabb-e

#runFunction ouput:
#load("RunNetMHCpan_playground/results_list.RData")

InstantRecognMatrixII = function(results, type = c("Score_EL", "Rank_EL"), output_format = "long"){
  if (type[1] == "all") {type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm")}
  readin = results[[1]][grep("Sequence ",results[[1]])]
  peplist = sapply(readin, function(x) {unlist(strsplit(str_trim(x, side = "left"), "\\s+"))[3]}, USE.NAMES = FALSE) #Modified by Balazs
  
  # create long format
  mtx = as.data.frame(do.call(rbind, lapply(results, function(x){
    x = x[grep("Sequence ", x)]
    x = gsub("<=WB|<=SB", "", x) #Inserted by Balazs
    mat = do.call(rbind, strsplit(str_trim(x, side = "left"), "\\s+")) #Modified by Balazs
    mat[, c(2,3,8,9,11,12,13)] #Modified by Balazs
  })), stringsAsFactors = F)
  for (i in 3:7) {mtx[, i] = as.numeric(mtx[, i])}
  
  colnames(mtx) = c("allele", "peptide", "Score_EL", "Rank_EL", "Score_BA", "Aff_nm", "Rank_BA") #Modified by Balazs
  mtx = cbind.data.frame(mtx[, 1:2], mtx[, colnames(mtx) %in% type, drop = FALSE])
  
  # create wide format
  if(output_format == "wide"){
    mtx = lapply(type, function(x) acast(mtx, allele ~ peptide, value.var = x))
    names(mtx) = type
  }
  
  return(mtx)
}

