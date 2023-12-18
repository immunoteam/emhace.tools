\name{LoadResultsFromDir}
\alias{LoadResultsFromDir}
\title{Loads all NetMHCpan/NetMHCIIpan result text files from a directory into a string vector}
\usage{
LoadResultsFromDir(x)
}
\arguments{
  \item{x}{the path to the directory containing prediction result files.}
}
\description{The function loads text-file outputs of either NetMHCpan or NetMHCIIpan into a list of character vectors, which can be further processed using \code{CollectBindingResults} function.}

\value{
A list of character vectors, where each list element belongs to a loaded file, and each element of the string vectors correspond to a line in the output text file.
}


\examples{
result_files <- LoadResultsFromDir("/home/user/results/")
}