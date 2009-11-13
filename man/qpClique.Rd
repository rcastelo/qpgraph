\name{qpClique}
\alias{qpClique}

\title{
Complexity of the resulting qp-graphs
}
\description{
Calculates and plots the size of the largest maximal clique (the so-called
clique number or maximum clique size) as function of the non-rejection
rate.
}
\usage{
qpClique(nrrMatrix, N=NA, threshold.lim=c(0,1), breaks=5, plot=TRUE,
         exact.calculation=TRUE, approx.iter=100,
         qpCliqueOutput=NULL, density.digits=0,
         logscale.clqsize=FALSE,
         titleclq="maximum clique size as function of threshold",
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
  \item{exact.calculation}{logical; if TRUE then the exact clique number is
       calculated; if FALSE then a lower bound is given instead.}
  \item{approx.iter}{number of iterations to be employed in the calculation of
       the lower bound (i.e., only applies when \code{exact.calculation=FALSE}).}
  \item{qpCliqueOutput}{output from a previous call to \code{\link{qpClique}}.
       This allows one to plot the result changing some of the plotting
       parameters without having to do the calculation again.}
  \item{density.digits}{number of digits in the reported graph densities.}
  \item{logscale.clqsize}{logical; if TRUE then the scale for the maximum clique
       size is logarithmic which is useful when working with more than 1000
       variables; FALSE otherwise (default).}
  \item{titleclq}{main title to be shown in the plot.}
  \item{verbose}{show progress on calculations.}
}
\details{
The estimate of the complexity of the resulting qp-graphs is calculated as
the area enclosed under the curve of maximum clique sizes.

The maximum clique size, or clique number, is obtained by calling the function
\code{\link{qpCliqueNumber}}
The calculation of the clique number of an undirected graph is an NP-complete
problem which means that its computational cost is bounded by an exponential
running time (Pardalos and Xue, 1994). Therefore, giving breakpoints between
0.95 and 1.0 may result into very dense graphs which can lead to extremely
long execution times. If it is necessary to look at that range of breakpoints
it is recommended either to use the lower bound on the clique number
(\code{exact.calculation=FALSE}) or to look at \code{\link{qpGraphDensity}}.
}
\value{
A list with the maximum clique size and graph density as function of threshold,
an estimate of the complexity of the resulting qp-graphs across the thresholds,
the threshold on the non-rejection rate that provides a maximum clique size
strictly smaller than the sample size N and the resulting maximum clique size.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Pardalos, P.M. and Xue, J. The maximum clique problem.
\emph{J. Global Optim.}, 4:301-328, 1994.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpCliqueNumber}}
  \code{\link{qpGraphDensity}}
}
\examples{
require(mvtnorm)

nVar <- 50  ## number of variables
maxCon <- 5 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

A <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=Sigma)

## the higher the q the less complex the qp-graph

nrr.estimates <- qpNrr(X, q=1, verbose=FALSE)

qpClique(nrr.estimates, plot=FALSE)$complexity

nrr.estimates <- qpNrr(X, q=5, verbose=FALSE)

qpClique(nrr.estimates, plot=FALSE)$complexity
}
\keyword{models}
\keyword{multivariate}
