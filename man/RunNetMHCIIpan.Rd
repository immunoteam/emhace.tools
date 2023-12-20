\name{RunNetMHCIIpan}
\alias{RunNetMHCIIpan}
\title{Runs MHC-II binding predictions using NetMHCIIpan and collects the outputs of the software.}
\usage{
RunNetMHCIIpan(alleles, peptides, paired_input = F,
             value_type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"),
             output_format = "long", threads = 1, version_number = "4.1",
             result_files_location = NULL, keep_pep = FALSE,
             software_path = NULL, tmppep_loc = NULL)
}
\arguments{
  \item{alleles}{a character vector of MHC-II alleles (i.e. DRB1_0101, HLA-DQA10101-DQB10101)}
  \item{value_type}{the types of binding values to be kept in the output of the function. If \code{NULL}, no output will be returned inside R, only the result files if specified in \code{result_files_location}.}
  \item{paired_input}{logical; if \code{TRUE}, predictions will be performed on corresponding allele-peptide pairs. In this case, vectors \code{alleles} and \code{peptides} must be of the same length.}
  \item{output_format}{determines the format of the output object, whether it should be a long data frame ("long") or a matrix ("wide").}
  \item{threads}{the number of threads used for the predictions.}
  \item{version_number}{a character indicating the version number of NetMHCIIpan.}
  \item{result_files_location}{path to a directory where the NetMHCIIpan result files (txt) will be stored.}
  \item{keep_pep}{logical, determines whether to keep the temporary pepfile (default: \code{FALSE}).}
  \item{software_path}{the full path to the executable file of NetMHCIIpan.}
  \item{tmppep_loc}{the path where the temporary pepfile should be stored. The default is the directory derived from \code{software_path}.}
}
\description{The function performs the binding predictions between MHC molecules and peptides.}

\value{
The output is either a list of matrices (\code{output_format = "wide"}), containing data for separate binding parameters, or a long data frame showing information for allele-peptide pairs (\code{output_format = "long"}). If \code{output_format} is \code{NULL}, then only raw result files are stored in a directory specified in \code{result_files_location}
}

\examples{
# returns a long data frame containing all combinations of alleles and peptides (output not shown)
RunNetMHCIIpan(alleles = c("DRB1_0101", "HLA-DQA10401-DQB10307"),
             peptides = c("KLKLKLKLKL", "ALALKLAKLKA", "PQPQPQPQPQ"), version_number = "4.3",
             software_path = "/home/lhgergo/Programok/netMHCIIpan-4.3/netMHCIIpan")

# returns a wide matrix containing all combinations of alleles and peptides
RunNetMHCIIpan(alleles = c("DRB1_0101", "HLA-DQA10401-DQB10307"),
             peptides = c("KLKLKLKLKL", "ALALKLAKLKA", "PQPQPQPQPQ"), version_number = "4.3",
             output_format = "wide", value_type = "Rank_EL",
             software_path = "/home/lhgergo/Programok/netMHCIIpan-4.3/netMHCIIpan")
    
# output:
# $Rank_EL
#                       ALALKLAKLKA KLKLKLKLKL PQPQPQPQPQ
# DRB1_0101                   98.00        100     100.00
# HLA-DQA10401-DQB10307       99.48        100      99.92
}
