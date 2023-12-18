\name{ReduceEpitopeSimilarity}
\alias{ReduceEpitopeSimilarity}
\title{Iteratively removes similar peptides to create a set of diverse amino acid sequences}
\usage{
ReduceEpitopeSimilarity(epitopes, similarity = 0.5, threads = 7, clustalo_path,
                        using_windows = TRUE, keep_temp = F,
                        tempdir = getwd(), save_histograms = T, verbose = F)
}
\arguments{
  \item{epitopes}{a character vector containing the initial set of peptide sequences.}
  \item{similarity}{the similarity value to be reached during the iterative elimination of similar sequences}
  \item{threads}{the number of threads to be used by the algorithm}
  \item{clustalo_path}{the path to the Clustal Omega software}
  \item{using_windows}{logical, should be set to \code{TRUE} if the function is used in a Windows-base environment.}
  \item{keep_temp}{logical, if \code{TRUE}, the temporary Clustal Omega output file will be kept}
  \item{tempdir}{the directory where the temporary file should be stored. If not defined, the current working directory will be used.}
  \item{save_histograms}{logical, determining whether to save the histogram visualizing the changes in similarity values in the peptide set after the process.}
}
\description{The function iteratively removes similar peptides from an input peptide set until it reaches a pre-defined minimum similarity threshold to create a diverse set of amino acid sequences.}
\details{
In the process, highly similar peptide sequences are identified using Clustal Omega and excluded as follows. A protein distance matrix is being created and peptides are discarded iteratively. In each iteration, the peptide pairs with the lowest k-tuple distance are identified. Then, the peptide with the highest average similarity to all other sequences is excluded. Iterations are repeated until distance values less than 0.5 (corresponding to greater than approximately 50\% sequence identity) were eliminated from the matrix.
}
\value{
The function returns a reduced peptide set (as a character vector) and two \code{ggplot2} histograms indicating the distribution of similarity values between peptides before and after the elimination process.
}

\references{
Manczinger, Máté, et al. "Pathogen diversity drives the evolution of generalist MHC-II alleles in human populations." PLoS biology 17.1 (2019): e3000131.

Yang, Kuan, and Liqing Zhang. "Performance comparison between k-tuple distance and four model-based distances in phylogenetic tree reconstruction." Nucleic acids research 36.5 (2008): e33.
}


\examples{
\emph{To be written...}
}