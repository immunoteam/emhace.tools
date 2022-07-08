### RunCommand
### starts commands (stored in a character vector) on one or multiple threads

RunCommand <- function(cmds, threads = 1, intern = F) {
  if(threads == 1) {
    lapply(cmds, function(cmd) {system(cmd, intern = intern)})
  } else if (threads > 1) {
    future::plan("multisession", workers = threads)
    out <- furrr::future_map(cmds, ~system(.x, intern = intern))
    future:::ClusterRegistry("stop")
    return(out)
  } else {
    message("Number of threads cannot be smaller than 1.")
  }
}
