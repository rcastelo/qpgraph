\name{qpAvgNrr}
\alias{qpAvgNrr}
\alias{qpAvgNrr,ExpressionSet-method}
\alias{qpAvgNrr,data.frame-method}
\alias{qpAvgNrr,matrix-method}

\title{
Average non-rejection rate estimation
}
\description{
Estimates average non-rejection rates for every pair of variables.
}
\usage{
\S4method{qpAvgNrr}{ExpressionSet}(X, qOrders=4, nTests=100, alpha=0.05,
                                   pairup.i=NULL, pairup.j=NULL,
                                   type=c("arith.mean"), verbose=TRUE,
                                   identicalQs=TRUE, R.code.only=FALSE,
                                   clusterSize=1)
\S4method{qpAvgNrr}{data.frame}(X, qOrders=4, nTests=100, alpha=0.05,
                                pairup.i=NULL, pairup.j=NULL,
                                long.dim.are.variables=TRUE,
                                type=c("arith.mean"), verbose=TRUE,
                                identicalQs=TRUE, R.code.only=FALSE,
                                clusterSize=1)
\S4method{qpAvgNrr}{matrix}(X, qOrders=4, nTests=100, alpha=0.05,
                            pairup.i=NULL, pairup.j=NULL,
                            long.dim.are.variables=TRUE,
                            type=c("arith.mean"), verbose=TRUE,
                            identicalQs=TRUE, R.code.only=FALSE,
                            clusterSize=1)
}
\arguments{
  \item{X}{data set from where to estimate the average non-rejection rates.
       It can be an ExpressionSet object, a data frame or a matrix.}
  \item{qOrders}{either a number of partial-correlation orders or a vector of
       vector of particular orders to be employed in the calculation.}
  \item{nTests}{number of tests to perform for each pair for variables.}
  \item{alpha}{significance level of each test.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when the data is a data frame or a matrix, the longer dimension
       is the one defining the random variables; if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
  \item{type}{type of average. By now only the arithmetic mean is available.}
  \item{verbose}{show progress on the calculations.}
  \item{identicalQs}{use identical conditioning subsets for every pair of vertices
       (default), otherwise sample a new collection of \code{nTests} subsets for
       each pair of vertices.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
  \item{clusterSize}{size of the cluster of processors to employ if we wish to
       speed-up the calculations by performing them in parallel. A value of 1
       (default) implies a single-processor execution. The use of a cluster of
       processors requires having previously loaded the packages \code{snow}
       and \code{rlecuyer}.}
}
\details{
Note that when specifying a vector of particular orders \code{q}, these values
should be in the range 1 to \code{min(p,n-3)}, where \code{p} is the number of
variables and \code{n} the number of observations. The computational cost
increases linearly within each \code{q} value and quadratically in \code{p}.
When setting \code{identicalQs} to \code{FALSE} the computational cost may
increase between 2 times and one order of magnitude (depending on \code{p} and
\code{q}) while asymptotically the estimation of the non-rejection rate
converges to the same value.
}
\value{
A \code{\link{dspMatrix-class}} symmetric matrix of estimated average
non-rejection rates with the diagonal set to \code{NA} values. When using the
arguments \code{pairup.i} and \code{pairup.j}, those cells outside the
constraint pairs will get also a \code{NA} value.
}
\references{
Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
networks from microarray data with qp-graphs. \emph{J. Comp. Biol.},
16(2):213-227, 2009.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
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

avgnrr.estimates <- qpAvgNrr(X, verbose=FALSE)

## distribution of average non-rejection rates for the present edges
summary(avgnrr.estimates[upper.tri(avgnrr.estimates) & A])

## distribution of average non-rejection rates for the missing edges
summary(avgnrr.estimates[upper.tri(avgnrr.estimates) & !A])

\dontrun{
library(snow)
library(rlecuyer)

## only for moderate and large numbers of variables the
## use of a cluster of processors speeds up the calculations

nVar <- 500
A <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

system.time(avgnrr.estimates <- qpAvgNrr(X, q=5, verbose=TRUE))
system.time(avgnrr.estimates <- qpAvgNrr(X, q=5, verbose=TRUE, clusterSize=8))
}
}
\keyword{models}
\keyword{multivariate}
