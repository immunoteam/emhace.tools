\name{IsValidPeptide}
\alias{IsValidPeptide}
\title{Checks if the input peptide sequence is suitable for the chosen version of NetMHCpan/NetMHCIIpan}
\usage{
IsValidPeptide(peptides, hla_type = 1)
}
\arguments{
  \item{peptides}{a vector of peptide sequences.}
  \item{hla_type}{an integer indicating whether the peptides should be valid candidates for either MHC-I (1) or MHC-II (2) predictions.}
}
\description{The function checks if each peptide in a character-type vector is suitable for MHC-binding predictions, containing only supported proteogenic amino acids and peptide lengths.}

\value{
The output is a logical vector with the same length of the input vector \code{peptides}, where \code{TRUE} values indicate that the peptide is a valid input for prediction.
}

\examples{
peptides <- c("KLKLKLKLK", "KLKLKLKLU", "KLKLKLKLKKLKLKLKLK")
IsValidPeptide(peptides, hla_type = 1)
# output:
# [1]  TRUE FALSE FALSE
}