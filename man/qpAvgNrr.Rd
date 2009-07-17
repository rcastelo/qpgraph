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
\S4method{qpAvgNrr}{ExpressionSet}(data, qOrders=4, nTests=100, alpha=0.05,
                                   pairup.i=NULL, pairup.j=NULL,
                                   type=c("arith.mean"), verbose=TRUE,
                                   identicalQs=FALSE, R.code.only=FALSE)
\S4method{qpAvgNrr}{data.frame}(data, qOrders=4, nTests=100, alpha=0.05,
                                pairup.i=NULL, pairup.j=NULL,
                                long.dim.are.variables=TRUE,
                                type=c("arith.mean"), verbose=TRUE,
                                identicalQs=FALSE, R.code.only=FALSE)
\S4method{qpAvgNrr}{matrix}(data, qOrders=4, nTests=100, alpha=0.05,
                            pairup.i=NULL, pairup.j=NULL,
                            long.dim.are.variables=TRUE,
                            type=c("arith.mean"), verbose=TRUE,
                            identicalQs=FALSE, R.code.only=FALSE)
}
\arguments{
  \item{data}{data set from where to estimate the average non-rejection rates.
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
  \item{identicalQs}{experimental feature.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
Note that when specifying a vector of particular orders \code{q}, these values
should be in the range 1 to \code{min(p,n-3)}, where \code{p} is the number of
variables and \code{n} the number of observations. The computational cost
increases linearly within each \code{q} value and quadratically in \code{p}.
}
\value{
A symmetric matrix of estimated average non-rejection rates.
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
nVar <- 50 # number of variables
maxCon <- 5  # maximum connectivity per variable
nObs <- 30 # number of observations to simulate

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

X <- qpSampleMvnorm(K, nObs)

avgnrr.estimates <- qpAvgNrr(X, verbose=FALSE)

summary(avgnrr.estimates[upper.tri(avgnrr.estimates) & I])

summary(avgnrr.estimates[upper.tri(avgnrr.estimates) & !I])

}
\keyword{models}
\keyword{multivariate}
