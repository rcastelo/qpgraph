\name{qpI2K}
\alias{qpI2K}

\title{
Random concentration matrix
}
\description{
Builds a random concentration matrix containing zeroes on those entries
associated to pairs of variables that are disconnected on a given undirected
graph.
}
\usage{
qpI2K(I, verbose=FALSE, R.code.only=FALSE)
}
\arguments{
  \item{I}{incidence matrix of an undirected graph.}
  \item{verbose}{show progress on the calculations.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
The random concentration matrix is built by first generating a matrix of
random correlations using the method from Marsaglia and Oltkin (1984). Second,
this matrix is inverted to obtain an initial random covariance matrix. Third,
this covariance matrix is adjusted to the independence constraints of the input
undirected graph by using the function \code{\link{qpIPF}} and finally is
inverted to obtain the final random concentration matrix.
}
\value{
A random concentration matrix with zeroes at the empty adjacencies of the
undirected graph defined by the input incidence matrix.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpSampleMvnorm}}
  \code{\link{qpK2R}}
}
\examples{
nVar <- 50 # number of variables
maxCon <- 5  # maximum connectivity per variable

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

realI <- K != 0
diag(realI) <- FALSE

sum(realI) / 2

sum(I) / 2

# all present edges (dependencies) in realI must be in I
identical(I & realI, realI)

# all missing edges (independencies) in I must be in realI
identical(!I & !realI, !I)

}
\keyword{models}
\keyword{multivariate}
