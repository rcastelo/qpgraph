\name{qpNrr}
\alias{qpNrr}
\alias{qpNrr,ExpressionSet-method}
\alias{qpNrr,data.frame-method}
\alias{qpNrr,matrix-method}

\title{
Non-rejection rate estimation
}
\description{
Estimates non-rejection rates for every pair of variables.
}
\usage{
\S4method{qpNrr}{ExpressionSet}(X, q=1, nTests=100, alpha=0.05, pairup.i=NULL,
                                pairup.j=NULL, verbose=TRUE, identicalQs=TRUE,
                                R.code.only=FALSE, clusterSize=1)
\S4method{qpNrr}{data.frame}(X, q=1, nTests=100, alpha=0.05, pairup.i=NULL,
                             pairup.j=NULL, long.dim.are.variables=TRUE, verbose=TRUE,
                             identicalQs=TRUE, R.code.only=FALSE, clusterSize=1)
\S4method{qpNrr}{matrix}(X, q=1, nTests=100, alpha=0.05, pairup.i=NULL,
                         pairup.j=NULL, long.dim.are.variables=TRUE, verbose=TRUE,
                         identicalQs=TRUE, R.code.only=FALSE, clusterSize=1)
}
\arguments{
  \item{X}{data set from where to estimate the non-rejection rates.
       It can be an ExpressionSet object, a data frame or a matrix.}
  \item{q}{partial-correlation order to be employed.}
  \item{nTests}{number of tests to perform for each pair for variables.}
  \item{alpha}{significance level of each test.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when data are in a data frame or in a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
  \item{verbose}{show progress on the calculations.}
  \item{identicalQs}{use identical conditioning subsets for every pair of vertices
       (default), otherwise sample a new collection of \code{nTests} subsets for
       each pair of vertices.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
  \item{clusterSize}{size of the cluster of processors to employ if we wish to
       speed-up the calculations by performing them in parallel. A value of 1
       (default) implies a single-processor execution. The use of a cluster of
       processors requires having previously loaded the packages \code{Rmpi}
       and \code{snow}.}
}
\details{
Note that the possible values of \code{q} should be in the range 1 to
\code{min(p,n-3)}, where \code{p} is the number of variables and \code{n}
the number of observations. The computational cost increases linearly with
\code{q} and quadratically in \code{p}. When setting \code{identicalQs} to
\code{FALSE} the computational cost may increase between 2 times and one
order of magnitude (depending on \code{p} and \code{q}) while asymptotically
the estimation of the non-rejection rate converges to the same value.
}
\value{
A \code{\link{dspMatrix-class}} symmetric matrix of estimated non-rejection
rates with the diagonal set to \code{NA} values. When using the arguments
\code{pairup.i} and \code{pairup.j}, those cells outside the constraint pairs
will get also a \code{NA} value.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpAvgNrr}}
  \code{\link{qpEdgeNrr}}
  \code{\link{qpHist}}
  \code{\link{qpGraphDensity}}
  \code{\link{qpClique}}
}
\examples{
require(mvtnorm)

nVar <- 50  ## number of variables
maxCon <- 5 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

A <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

nrr.estimates <- qpNrr(X, q=5, verbose=FALSE)

## distribution of non-rejection rates for the present edges
summary(nrr.estimates[upper.tri(nrr.estimates) & A])

## distribution of non-rejection rates for the missing edges
summary(nrr.estimates[upper.tri(nrr.estimates) & !A])

\dontrun{
library(Rmpi)
library(snow)

## only for moderate and large numbers of variables the
## use of a cluster of processors speeds up the calculations

nVar <- 500
A <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

system.time(nrr.estimates <- qpNrr(X, q=5, verbose=FALSE))
system.time(nrr.estimates <- qpNrr(X, q=5, verbose=FALSE, clusterSize=8))
}
}
\keyword{models}
\keyword{multivariate}
