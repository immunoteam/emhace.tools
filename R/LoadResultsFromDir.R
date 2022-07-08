LoadResultsFromDir <- function(x) {
  file_paths <- list.files(x, full.names = TRUE)
  lapply(file_paths, function(x) {readLines(x)})
}
