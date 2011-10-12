\name{qpBoundary}
\alias{qpBoundary}

\title{
Maximum boundary size of the resulting qp-graphs
}
\description{
Calculates and plots the size of the largest vertex boundary
as function of the non-rejection rate.
}
\usage{
qpBoundary(nrrMatrix, N=NA, threshold.lim=c(0,1), breaks=5, plot=TRUE,
           qpBoundaryOutput=NULL, density.digits=0,
           logscale.bdsize=FALSE,
           titlebd="Maximum boundary size as function of threshold",
           verbose=FALSE)
}
\arguments{
  \item{nrrMatrix}{matrix of non-rejection rates.}
  \item{N}{number of observations from where the non-rejection rates were
       estimated.}
  \item{threshold.lim}{range of threshold values on the non-rejection rate.}
  \item{breaks}{either a number of threshold bins or a vector of threshold
       breakpoints.}
  \item{plot}{logical; if TRUE makes a plot of the result; if FALSE it does not.}
  \item{qpBoundaryOutput}{output from a previous call to \code{\link{qpBoundary}}.
       This allows one to plot the result changing some of the plotting
       parameters without having to do the calculation again.}
  \item{density.digits}{number of digits in the reported graph densities.}
  \item{logscale.bdsize}{logical; if TRUE then the scale for the maximum boundary
       size is logarithmic which is useful when working with more than 1000
       variables; FALSE otherwise (default).}
  \item{titlebd}{main title to be shown in the plot.}
  \item{verbose}{show progress on calculations.}
}
\details{
The maximum boundary is calculated as the largest degree among all vertices of
a given qp-graph.
}
\value{
A list with the maximum boundary size and graph density as function of threshold,
the threshold on the non-rejection rate that provides a maximum boundary size
strictly smaller than the sample size N and the resulting maximum boundary size.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpHTF}}
  \code{\link{qpGraphDensity}}
}
\examples{
require(mvtnorm)

nVar <- 50  ## number of variables
maxCon <- 5 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

A <- qpRndGraph(p=nVar, d=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

## the higher the q the less complex the qp-graph

nrr.estimates <- qpNrr(X, q=1, verbose=FALSE)

qpBoundary(nrr.estimates, plot=FALSE)

nrr.estimates <- qpNrr(X, q=5, verbose=FALSE)

qpBoundary(nrr.estimates, plot=FALSE)
}
\keyword{models}
\keyword{multivariate}
