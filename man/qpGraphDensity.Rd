\name{qpGraphDensity}
\alias{qpGraphDensity}

\title{
Densities of resulting qp-graphs
}
\description{
Calculates and plots the graph density as function of the non-rejection rate.
}
\usage{
qpGraphDensity(nrrMatrix, threshold.lim=c(0,1), breaks=5,
               plot=TRUE, qpGraphDensityOutput=NULL,
               density.digits=0,
               titlegd="graph density as function of threshold")

}
\arguments{
  \item{nrrMatrix}{matrix of non-rejection rates.}
  \item{threshold.lim}{range of threshold values on the non-rejection rate.}
  \item{breaks}{either a number of threshold bins or a vector of threshold
       breakpoints.}
  \item{plot}{logical; if TRUE makes a plot of the result; if FALSE it does not.}
  \item{qpGraphDensityOutput}{output from a previous call to
       \code{\link{qpGraphDensity}}. This allows one to plot the result changing
       some of the plotting parameters without having to do the calculation
       again.}
  \item{density.digits}{number of digits in the reported graph densities.}
  \item{titlegd}{main title to be shown in the plot.}
}
\details{
The estimate of the sparseness of the resulting qp-graphs is calculated as
one minus the area enclosed under the curve of graph densities.
}
\value{
A list with the graph density as function of threshold and an estimate of the
sparseness of the resulting qp-graphs across the thresholds.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
  \code{\link{qpAvgNrr}}
  \code{\link{qpEdgeNrr}}
  \code{\link{qpClique}}
}
\examples{
nVar <- 50 # number of variables
maxCon <- 5  # maximum connectivity per variable
nObs <- 30 # number of observations to simulate

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

X <- qpSampleMvnorm(K, nObs)

# the higher the q the sparser the qp-graphs

nrr.estimates <- qpNrr(X, q=1, verbose=FALSE)

qpGraphDensity(nrr.estimates, plot=FALSE)$sparseness

nrr.estimates <- qpNrr(X, q=5, verbose=FALSE)

qpGraphDensity(nrr.estimates, plot=FALSE)$sparseness

}
\keyword{models}
\keyword{multivariate}
