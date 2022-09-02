#' @export

### GetNetMHCIIpan function
### helps in fethching, installation and configuration of netMHCIIpan

GetNetMHCIIpan <- function() {
  # checking if the user performed the function call from Linux
  if(Sys.info()['sysname'] != "Linux") {
    stop("Installation is only available on Linux systems!")
  }
  
  # checking if tcsh is installed on the Linux system (dependency of netMHCIIpan 4.0)
  if(class(try(system("command -v tcsh", intern = T), silent = T)) == "try-error") {
    stop("Installation failed: tcsh is not installed on this Linux system.\nPlease use your system package manager (apt, yum..) to download tcsh!")
  }
  
  # setting download path
  message("Welcome! This function helps you download and install netMHCIIpan.")
  message("During the installation process please always use full paths!")
  download_path <- readline("Set the directory where netMHCIIpan installation files should be downloaded: ")
  suppressWarnings(dir.create(download_path, recursive = T))
  
  # opening browser to fill in the form
  message("On the webpage appearing soon please choose netMHCIIpan 4.0 for download, and fill in the form!")
  message("Please use your institutional e-mail adress for registration! You will get a link to an FTP repository in an e-mail.")
  message("When you receive your e-mail, paste the link from the mail in this console below!")
  readline("Press <ENTER> to open the download website")
  browseURL("https://services.healthtech.dtu.dk/software.php")

  # downloading main archive
  ftplink <- readline("Link from the _e-mail_: ")
  download.file(paste0(ftplink, "/netMHCIIpan-4.0.Linux.tar.gz"), destfile = paste0(download_path, "/netmhcIIpan.tar.gz"))
  
  # extracting main archive
  install_path <- readline("Set the installation directory, where you would like your 'netMHCIIpan' directory to reside:")
  suppressWarnings(dir.create(install_path, recursive = T))
  system(paste0("tar -xf ",  download_path, "/netmhcIIpan.tar.gz -C ", install_path))
  message(paste0("Succesfully uncompressed archive to ", install_path))
  
  # setting netMHCIIpan script up
  netmhcpan_file <- readLines(paste0(install_path, "/netMHCIIpan-4.0/netMHCIIpan"))
  netmhcpan_file[14] <- paste0("setenv NMHOME ", install_path, "/netMHCIIpan-4.0")
  tmpdir_path_default <- paste0(install_path, "/netMHCIIpan-4.0/tmp/")
  tmpdir_path <- readline(paste0("Set a directory for temporary files, used by netMHCpan (for default press <ENTER>: ", install_path, "/netMHCIIpan-4.0/tmp/)"))
  tmpdir_path <- ifelse(tmpdir_path == "", tmpdir_path_default, tmpdir_path)
  suppressWarnings(dir.create(tmpdir_path, recursive = T))
  netmhcpan_file[19] <- paste0("	setenv  TMPDIR  ", tmpdir_path)
  write(netmhcpan_file, file = paste0(install_path, "/netMHCIIpan-4.0/netMHCIIpan"))
  
  # downloading and extracting data archive
  message("Downloading and uncompressing data needed by netMHCIIpan-4.0")
  download.file("https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/data.tar.gz", destfile = paste0(install_path, "/netMHCIIpan-4.0/data.tar.gz"))
  system(paste0("tar -xf ", install_path, "/netMHCIIpan-4.0/data.tar.gz -C ", install_path, "/netMHCIIpan-4.0/"))
  file.remove(paste0(install_path, "/netMHCIIpan-4.0/data.tar.gz"))
  
  # running test
  message("Performing prediction on a test file")
  resfile_crnt <- system(paste0(install_path, "/netMHCIIpan-4.0/netMHCIIpan -inptype 1 -f ", install_path, "/netMHCIIpan-4.0/test/example.pep"), intern = T)
  resfile_test <- readLines(paste0(install_path, "/netMHCIIpan-4.0/test/example.pep.out"))
  
  same_output_lgl <- all(resfile_crnt[14:(length(resfile_crnt)-4)] == resfile_test[14:(length(resfile_test)-4)])
  if(same_output_lgl) {
    message("Test passed. Successfully installed netMHCIIpan 4.0!")
  } else {
    message("Test and installation of netMHCIIpan 4.0 failed!")
  }
}
