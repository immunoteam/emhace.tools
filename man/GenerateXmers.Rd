\name{GenerateXmers}
\alias{GenerateXmers}
\title{Generating overlapping subsequences of pre-defined lengths from a protein sequence}
\usage{
GenerateXmers(x, l, string_input = TRUE)
}
\arguments{
  \item{x}{the input protein sequence to be splitted into subsequences; the input type depends on the value of \code{string_input})}
  \item{l}{a single numeric value, indicating the length of subsequences to be generated.}
  \item{string_input}{if \code{TRUE}, GenerateXmers accepts input protein sequence (\code{x}) as a single character variable; otherwise it requires the input sequence in a character vector format.}
}
\description{The function generates all possible peptide fragments of pre-defined length from a protein sequence.}

\value{
The output is a character vector of short peptide sequences with lengths determined in parameter \code{l}.
}

\examples{
# generating all possible nonamers based on a character vector of amino acids
print(protseq)
# output:
# [1] "K" "L" "T" "D" "V" "I" "T" "H" "I" "Y" "Q" "Y" "G" "Y" "T" "G" "A" "Q" "F" "G"
GenerateXmers(protseq, 9, string_input = FALSE)
# output:
# [1] "KLTDVITHI" "LTDVITHIY" "TDVITHIYQ" "DVITHIYQY" "VITHIYQYG" "ITHIYQYGY" "THIYQYGYT" "HIYQYGYTG" "IYQYGYTGA" "YQYGYTGAQ" "QYGYTGAQF" "YGYTGAQFG"

# generating all possible nonamers based on a sequence stored in a single character string
print(protseq)
# output:
# [1] "YNFMMINHHYQHRTPWVGLH"
GenerateXmers(protseq, 9, string_input = TRUE)
# output:
# [1] "YNFMMINHH" "NFMMINHHY" "FMMINHHYQ" "MMINHHYQH" "MINHHYQHR" "INHHYQHRT" "NHHYQHRTP" "HHYQHRTPW" "HYQHRTPWV" "YQHRTPWVG" "QHRTPWVGL" "HRTPWVGLH"

# splitting up sequences into subsequences of different lengths in one turn
GenerateXmersMultilength <- Vectorize(GenerateXmers, vectorize.args = "l")
GenerateXmersMultilength(protseq, 8:12)
# output: a 5-element list containing 8, 9... 12 amino acid long subsequences of the input protein
}
