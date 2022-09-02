#epitope-allele pair predictor
setwd("AlleleEpitopePairPredictor_playground")

#necessary functions

library(fastmatch)
row.fmatch = function (x, table, nomatch = NA) 
{
  if (class(table) == "matrix") 
    table <- as.data.frame(table)
  if (is.null(dim(x))) 
    x <- as.data.frame(matrix(x, nrow = 1))
  cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
  fmatch(cx, ct, nomatch = nomatch)
}


#Input: chr. vectors:
  #alleles
  #epitopes
  #result_dir
#First, create prediction files


#how to create temporary files:
writeLines(missing_peptides, con = )
writeLines(missing_alleles, con = )


#runNetMHCPan

#Recogn matrices in result dir






