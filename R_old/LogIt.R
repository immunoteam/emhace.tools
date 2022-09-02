# LogIt
# pastes a new entry in a logfile, marked with timestamp and function name causing the issue
LogIt <- function(msg, funcname, file, action = "silent") {
  write(paste0("# ", Sys.time(), " ", funcname, " ", action, "\n", msg), file = file, append = TRUE)
}

