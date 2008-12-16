\name{qpK2R}
\alias{qpK2R}

\title{
Partial correlation coefficients
}
\description{
Obtains partial correlation coefficients from a given concentration matrix.
}
\usage{
qpK2R(K)
}
\arguments{
  \item{K}{concentration matrix.}
}
\details{
Partial correlation coefficients are obtained from a concentration matrix by
scaling it by its diagonal and swapping the sign of the off-diagonal elements.
}
\value{
A matrix with partial correlation coefficients.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpI2K}}
}
\examples{
nVar <- 50 # number of variables
maxCon <- 5  # maximum connectivity per variable

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

R <- qpK2R(K)

# true partial correlation coefficients of the present edges
summary(abs(R[upper.tri(R) & I]))

# true partial correlation coefficients of the missing edges
summary(abs(R[upper.tri(R) & !I]))

}
\keyword{models}
\keyword{multivariate}
