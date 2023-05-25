#' @export
# PARAMETERS:
# x: list, containing the binding matrices (wide binding result objects) to be merged

# test codes ----------
# GenerateRandomPeptide <- function() {
#   paste0(sample(rownames(protr::AABLOSUM45), 9), collapse = "")
# }

# pepset_1 <- sapply(1:5, function(x) GenerateRandomPeptide())
# pepset_2 <- sapply(1:5, function(x) GenerateRandomPeptide())
# alleles_1 <- c("HLA-A02:01", "HLA-C07:02")
# alleles_2 <- c("HLA-A11:01", "HLA-B07:02")
# 
# outmtx_1 <- RunNetMHCpan(alleles = alleles_1, peptides = pepset_1, output_format = "wide", software_path = "/home/lhgergo/Programok/netMHCpan-4.1/netMHCpan")$Rank_EL
# outmtx_2 <- RunNetMHCpan(alleles = alleles_1, peptides = pepset_2, output_format = "wide", software_path = "/home/lhgergo/Programok/netMHCpan-4.1/netMHCpan")$Rank_EL
# outmtx_3 <- RunNetMHCpan(alleles = alleles_2, peptides = pepset_1, output_format = "wide", software_path = "/home/lhgergo/Programok/netMHCpan-4.1/netMHCpan")$Rank_EL
# outmtx_4 <- RunNetMHCpan(alleles = alleles_2, peptides = pepset_2, output_format = "wide", software_path = "/home/lhgergo/Programok/netMHCpan-4.1/netMHCpan")$Rank_EL
# 
# list(outmtx_1, outmtx_2, outmtx_3, outmtx_4, outmtx_4) %>% MergeBindingResults()
# list(outmtx_1, outmtx_2, outmtx_3, outmtx_4) %>% MergeBindingResults()

# the function itself ----------
MergeBindingResults <- function(x) {
  reshape2::acast(data = reshape2::melt(x), Var1 ~ Var2, value.var = "value", fun.aggregate = mean)
}
