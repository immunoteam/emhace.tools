### IsValidPeptide
### checks if the given peptide has a correct length and all the characters in it are valid ones
amino_acids <- rownames(protr::AABLOSUM45)
IsValidPeptide <- function(peptides, hla_type = 1) {
  if(hla_type == 1) {
    sapply(strsplit(peptides, ""), function(x) {all(x %in% amino_acids)}) & nchar(peptides) >= 8 & nchar(peptides) <= 14
  } else if(hla_type == 2) {
    sapply(strsplit(peptides, ""), function(x) {all(x %in% amino_acids)}) & nchar(peptides) >= 9
  } else {stop("incorrect HLA type")}
}