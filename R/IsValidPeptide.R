#' @export

IsValidPeptide <- function(peptides, hla_type = 1) {
  if(hla_type == 1) {
    sapply(strsplit(peptides, ""), function(x) {all(x %in% rownames(protr::AABLOSUM45))}) & nchar(peptides) >= 8 & nchar(peptides) <= 14
  } else if(hla_type == 2) {
    sapply(strsplit(peptides, ""), function(x) {all(x %in% rownames(protr::AABLOSUM45))}) & nchar(peptides) >= 9
  } else {stop("incorrect HLA type")}
}
