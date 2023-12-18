\name{CollectBindingResults}
\alias{CollectBindingResults}
\title{Collects the outputs of NetMHCpan/NetMHCIIpan in a long data frame or a wide matrix format.}
\usage{
CollectBindingResults(results, value_type = c("Score_EL", "Rank_EL", "Score_BA", "Rank_BA", "Aff_nm"), output_format = "long", hla_type = 1, version_number = "4.1")
}
\arguments{
  \item{results}{a list of character vectors containing the results of MHC-binding predictions.}
  \item{value_type}{the types of binding values to be kept in the output of the function.}
  \item{output_format}{determines the format of the output object, whether it should be a long data frame ("long") or a matrix ("wide").}
  \item{hla_type}{an integer indicating whether the input of the function is from the results of MHC-I (1) or MHC-II (2) binding predictions.}
  \item{version_number}{a character indicating the version number of NetMHCpan/NetMHCIIpan being used to produce the input of this function.}
}
\description{The function transforms the outputs of NetMHCpan/NetMHCIIpan stored in lists of character vectors into a data frame- ("long") or matrix-based ("wide") format.}

\value{
The output is either a list of matrices (\code{output_format = "wide"}), containing data for separate binding parameters, or a long data frame showing information for allele-peptide pairs (\code{output_format = "long"}).
}

\examples{
# collecting outputs of NetMHCpan from a directory called "results"
results_raw <- LoadResultsFromDir("results/")

# transforming raw prediction outputs into a wide matrix format
results_wide <- CollectBindingResults(results_raw, value_type = c("Rank_EL", "Rank_BA"), output_format = "wide", hla_type = 1, version_number = "4.1")

# transforming raw prediction outputs into a long data frame format
resultsdf <- CollectBindingResults(results_raw, value_type = c("Rank_EL", "Rank_BA"), output_format = "long", hla_type = 1, version_number = "4.1")
}
