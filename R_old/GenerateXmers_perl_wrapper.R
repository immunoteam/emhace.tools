# GenerateXmers.pl function
# an R function to use GenerateXmers.pl script from R using system calls

# ##### measuring epitope splitting runtime compared to the R function
# protseqs <- sapply(1:100, function(x) {sample(rownames(protr::AABLOSUM100), 1000, replace = TRUE) %>% paste0(collapse = "")})
# source("../SARS-CoV-2-immunoadaptation/functions.R")
# GenerateXmersStringVctr <- Vectorize(GenerateXmersString, vectorize.args = c("szeki", "merslngth"))
# 
# begin_1 <- Sys.time()
# peps_1 <- GenerateXmersStringVctr(protseqs, 9)
# end_1 <- Sys.time()
# end_1 - begin_1
# 
# sapply(1:40, function(x) {
#   begin_2 <- Sys.time()
#   peps_2 <- sapply(protseqs, function(x){GenerateXmersString(x, 9)})
#   end_2 <- Sys.time()
#   
#   
#   begin_3 <- Sys.time()
#   peps <- system("./GenerateXmers.pl 9", input = protseqs, intern = TRUE)
#   peps <- SplitAt(peps, which(peps == "--")) %>% lapply(function(x) {x[-1]})
#   end_3 <- Sys.time()
#   c(runtime_2 = end_2 - begin_2, runtime_3 = end_3-begin_3)
# }) %>% t() -> runtimes
# 
# boxplot(runtimes[, 1], runtimes[, 2], log = "y", ylab = "second", names = c("R", "perl"))
# 
# mean(runtimes[, 1]) / mean(runtimes[, 2])

##### 
GenerateXmers <- function(x, merlngth) {
  peps <- system(paste0("./GenerateXmers.pl ", merlngth), input = x, intern = TRUE)
  map(SplitAt(peps, which(peps == "--")), ~.x[-1])[-(sum(peps == "--")+1)]
}
GenerateXmers <- Vectorize(GenerateXmers, vectorize.args = "merlngth", SIMPLIFY = FALSE)