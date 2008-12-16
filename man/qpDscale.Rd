\name{qpDscale}
\alias{qpDscale}

\title{
Scale a matrix by its diagonal
}
\description{
Scales the values of a matrix by the values on the diagonal.
}
\usage{
qpDscale(V)
}
\arguments{
  \item{V}{matrix with numerical values.}
}
\details{
The scaling of the values of the matrix by its diagonal is made
using the outer product of the inverted square root of the diagonal.
}
\value{
A scaled matrix such that the diagonal values become the unit.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpK2R}}
}
\examples{
m <- matrix(1:9, nrow=3)

m

qpDscale(m)
}
\keyword{models}
\keyword{multivariate}
