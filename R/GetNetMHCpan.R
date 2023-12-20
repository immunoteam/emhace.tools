#' @export

### GetNetMHCpan function
### helps in fethching, installation and configuration of netMHCpan

GetNetMHCpan <- function(version_number = "4.1") {
  # checking if the user performed the function call from Linux
  if(Sys.info()['sysname'] != "Linux") {
    stop("Installation is only available on Linux systems!")
  }
  
  # checking if tcsh is installed on the Linux system (dependency of netMHCpan 4.1)
  if(class(try(system("command -v tcsh", intern = T), silent = T)) == "try-error") {
    stop("Installation failed: tcsh is not installed on this Linux system.\nPlease use your system package manager (apt, yum..) to download tcsh!")
  }
  
  # setting download path
  message("Welcome! This function helps you download and install netMHCpan.")
  message("During the installation process please always use full paths!")
  download_path <- readline("Set the directory where netMHCpan installation files should be downloaded: ")
  suppressWarnings(dir.create(download_path, recursive = T))
  
  # opening browser to fill in the form
  message(paste0("On the webpage appearing soon please choose netMHCpan ", version_number, " for download, and fill in the form!"))
  message("Please use your institutional e-mail adress for registration! You will get a link to an FTP repository in an e-mail.")
  message("When you receive your e-mail, paste the link from the mail in this console below!")
  readline("Press <ENTER> to open the download website")
  browseURL("https://services.healthtech.dtu.dk/software.php")
  
  # downloading main archive
  ftplink <- readline("Link from the _e-mail_: ")
  download_page <- rvest::read_html(ftplink)
  pkg_file <- download_page %>% rvest::html_elements("a") %>% rvest::html_attr("href") %>% extract(grepl("\\.tar.gz", .))
  
  download.file(paste0(ftplink, "/", pkg_file), destfile = paste0(download_path, "/netmhcpan.tar.gz"))
  
  # extracting main archive
  install_path <- readline("Set the installation directory, where you would like your 'netMHCpan' directory to reside:")
  suppressWarnings(dir.create(install_path, recursive = T))
  system(paste0("tar -xf ",  download_path, "netmhcpan.tar.gz -C ", install_path))
  message(paste0("Succesfully uncompressed archive to ", install_path))
  
  # setting netMHCpan script up
  sw_name <- paste0("netMHCpan-", version_number)
  netmhcpan_file <- readLines(paste0(install_path, "/", sw_name, "/netMHCpan"))
  netmhcpan_file[14] <- paste0("setenv NMHOME ", install_path, "/", sw_name)
  tmpdir_path_default <- paste0(install_path, "/", sw_name, "/tmp/")
  tmpdir_path <- readline(paste0("Set a directory for temporary files, used by netMHCpan (for default press <ENTER>: ", install_path, "/", sw_name, "/tmp/)"))
  tmpdir_path <- ifelse(tmpdir_path == "", tmpdir_path_default, tmpdir_path)
  suppressWarnings(dir.create(tmpdir_path, recursive = T))
  netmhcpan_file[19] <- paste0("	setenv  TMPDIR  ", tmpdir_path)
  write(netmhcpan_file, file = paste0(install_path, "/", sw_name, "/netMHCpan"))
  
  # downloading and extracting data archive
  message(paste0("Downloading and uncompressing data needed by netMHCpan ", version_number))
  if(version_number == "4.0") {
    download.file(paste0("https://services.healthtech.dtu.dk/services/NetMHCpan-", version_number, "/data.Linux.tar.gz"),
                  destfile = paste0(install_path, "/netMHCpan-", version_number, "/data.tar.gz"))
    system(paste0("tar -xf ", install_path, "/netMHCpan-", version_number, "/data.tar.gz -C ",
                  install_path, "/", sw_name))
    file.remove(paste0(install_path, "/", sw_name, "/data.tar.gz"))
  } else if(version_number == "4.1") {
    download.file(paste0("https://services.healthtech.dtu.dk/services/NetMHCpan-", version_number, "/data.tar.gz"),
                  destfile = paste0(install_path, "/netMHCpan-", version_number, "/data.tar.gz"))
    system(paste0("tar -xf ", install_path, "/netMHCpan-", version_number, "/data.tar.gz -C ",
                  install_path, "/", sw_name))
    file.remove(paste0(install_path, "/", sw_name, "/data.tar.gz"))
  }
  
  # running test
  message("Performing prediction on a test file")
  resfile_crnt <- system(paste0(install_path, "/", sw_name, "/netMHCpan -p ", install_path, "/", sw_name, "/test/test.pep"), intern = T)
  resfile_test <- readLines(paste0(install_path, "/", sw_name, "/test/test.pep.out"))
  
  same_output_lgl <- all(resfile_crnt[48:(length(resfile_crnt)-5)] == resfile_test[48:(length(resfile_test)-5)])
  if(same_output_lgl) {
    message("Test passed. Successfully installed netMHCpan ", version_number, "!")
  } else {
    message("Test and installation of netMHCpan ", version_number, " failed!")
  }
}
