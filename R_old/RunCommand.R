### RunCommand
### starts commands (stored in a character vector) on one or multiple threads

RunCommand <- function(cmds, threads = 1, intern = F) {
  if(threads == 1) {
    lapply(cmds, function(cmd) {system(cmd, intern = intern)})
  } else if (threads > 1) {
    cl <- makeCluster(threads)
    registerDoParallel(cl)
    out <- parLapply(cl, cmds, function(cmd) {system(cmd, intern = intern)})
    stopCluster(cl)
    return(out)
  } else {
    message("Number of threads cannot be smaller than 1.")
  }
}
