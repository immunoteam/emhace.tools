\name{RunNetMHCpan}
\alias{RunNetMHCpan}
\title{Runs MHC-I binding predictions using NetMHCpan and collects the outputs of the software.}
\usage{
RunNetMHCpan(alleles, peptides, paired_input = F,
             value_type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"),
             output_format = "long", threads = 1, version_number = "4.1",
             result_files_location = NULL, keep_pep = FALSE,
             software_path = NULL, tmppep_loc = NULL)
}
\arguments{
  \item{alleles}{a character vector of MHC-I alleles (i.e. HLA-A02:01)}
  \item{value_type}{the types of binding values to be kept in the output of the function. If \code{NULL}, no output will be returned inside R, only the result files if specified in \code{result_files_location}.}
  \item{paired_input}{logical; if \code{TRUE}, predictions will be performed on corresponding allele-peptide pairs. In this case, vectors \code{alleles} and \code{peptides} must be of the same length.}
  \item{output_format}{determines the format of the output object, whether it should be a long data frame ("long") or a matrix ("wide").}
  \item{threads}{the number of threads used for the predictions.}
  \item{version_number}{a character indicating the version number of NetMHCpan.}
  \item{result_files_location}{path to a directory where the NetMHCpan result files (txt) will be stored.}
  \item{keep_pep}{logical, determines whether to keep the temporary pepfile (default: \code{FALSE}).}
  \item{software_path}{the full path to the executable file of NetMHCpan.}
  \item{tmppep_loc}{the path where the temporary pepfile should be stored. The default is the directory derived from \code{software_path}.}
}
\description{The function performs the binding predictions between MHC molecules and peptides.}

\value{
The output is either a list of matrices (\code{output_format = "wide"}), containing data for separate binding parameters, or a long data frame showing information for allele-peptide pairs (\code{output_format = "long"}). If \code{output_format} is \code{NULL}, then only raw result files are stored in a directory specified in \code{result_files_location}
}

\examples{
# returns a long data frame containing all combinations of alleles and peptides (output not shown)
RunNetMHCpan(alleles = c("HLA-A02:01", "HLA-B07:02", "HLA-C03:01"),
             peptides = c("KLKLKLKLKL", "ALALKLAKLKA", "PQPQPQPQPQ"), version_number = "4.1",
            software_path = "/home/user/Programs/netMHCpan-4.1/netMHCpan")

# returns a wide matrix containing all combinations of alleles and peptides
RunNetMHCpan(alleles = c("HLA-A02:01", "HLA-B07:02", "HLA-C03:01"),
             peptides = c("KLKLKLKLKL", "ALALKLAKLKA", "PQPQPQPQPQ"), version_number = "4.1", 
             output_format = "wide", value_type = "Rank_EL",
             software_path = "/home/user/Programs/netMHCpan-4.1/netMHCpan")
             
# output:
# $Rank_EL
#             ALALKLAKLKA KLKLKLKLKL PQPQPQPQPQ
# HLA-A*02:01       4.531      2.175     73.750
# HLA-B*07:02      53.077     11.000     27.619
# HLA-C*03:01      75.000     33.455     65.000
}
