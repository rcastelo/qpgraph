\name{qpRndHMGM}
\alias{qpRndHMGM}

\title{
Random homogeneous mixed graphical Markov model
}
\description{
Builds a random homogeneous mixed graphical Markov model (experimental feature).
}
\usage{
qpRndHMGM(nDiscrete=1, nContinuous=3, d=2, mixedIntStrength=5, rho=0.5, G=NULL)
}
\arguments{
  \item{nDiscrete}{number of discrete variables.}
  \item{nContinuous}{number of continuous variables.}
  \item{d}{degree of every vertex.}
  \item{mixedIntStrength}{strength of the mixed interactions.}
  \item{rho}{marginal correlation of the quadratic interactions.}
  \item{G}{input graph, if we don't want the function to simulate one.}
}
\details{
This function builds a random homogeneous mixed graphical model. It uses
\code{\link{qpRndGraph}} to simulate a random d-regular graph and then
builds a set of parameters that encode the conditional independencies
encoded by the graph and the given number of discrete and continuous
vertices. This is still an experimental feature and by now it generates
only models where the discrete variables are marginally independent.

}
\value{
A list with the graph and the parameters of the homogeneous mixed graphical
model, ready to be used with the function \code{\link{qpSampleFromHMGM}} for
sampling synthetic data using this model.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpRndGraph}}
  \code{\link{qpSampleFromHMGM}}
}
\examples{
## qpRndHMGM()
}
\keyword{models}
\keyword{multivariate}
