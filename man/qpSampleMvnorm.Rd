\name{qpSampleMvnorm}
\alias{qpSampleMvnorm}

\title{
Sample multivariate normal observations
}
\description{
Samples independent observations from a multivariate normal distribution with
a given mean vector and a given concentration matrix.
}
\usage{
qpSampleMvnorm(K, N, mean = rep(0, nrow(K)))
}
\arguments{
  \item{K}{concentration matrix of the multivariate normal distribution.}
  \item{N}{number of observations to sample.}
  \item{mean}{mean vector of the multivariate normal distribution.}
}
\details{
This function requires the \code{mvtnorm} package. This function is designed
to be used to sample multivariate normal observations from a randomly
generated concentration matrix that has a zero-structure reflecting the
conditional independencies from a, possibly also randomly generated, undirected
graph. The function \code{\link{qpI2K}} can be employed to obtain such a
concentration matrix.
}
\value{
A matrix where rows correspond to observations and columns to random variables.
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
nObs <- 30 # number of observations to simulate

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

X <- qpSampleMvnorm(K, nObs)

dim(X)
}
\keyword{models}
\keyword{multivariate}
