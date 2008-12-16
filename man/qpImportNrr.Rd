\name{qpImportNrr}
\alias{qpImportNrr}

\title{
Import non-rejection rates
}
\description{
Imports non-rejection rates from an external flat file.
}
\usage{
qpImportNrr(filename, nTests)
}
\arguments{
  \item{filename}{name of the flat file with the data on the non-rejection rates.}
  \item{nTests}{number of tests performed in the estimation of these
       non-rejection rates.}
}
\details{
This function expects a flat file with three tab-separated columns corresponding
to, respectively, 0-based index of one of the variables, 0-based index of the
other variable, number of non-rejected tests for the pair of variables of that
row in the text file. An example of a few lines of that file would be:

\preformatted{
6       3       95
6       4       98
6       5       23
7       0       94
7       1       94
}

After reading the file the function builds a matrix of non-rejection rates by
dividing the number of non-rejected tests by \code{nTests}. Note that if the flat
file to be imported would eventually have directly the rates instead of the
number of tests, these can be also imported by setting \code{nTests=1}.

This function is thought to be used to read files obtained from the standalone
parallel version of \code{\link{qpNrr}} which can be downloaded from
\url{http://functionalgenomics.upf.edu/qp}.
}
\value{
A symmetric matrix of non-rejection rates with the diagonal set to the \code{NA}
value.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
}
\keyword{models}
\keyword{multivariate}
