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
\S4method{qpNrr}{ExpressionSet}(data, q=1, nTests=100, alpha=0.05, pairup.i=NULL,
                                pairup.j=NULL, verbose=TRUE, identicalQs=FALSE,
                                R.code.only=FALSE)
\S4method{qpNrr}{data.frame}(data, q=1, nTests=100, alpha=0.05, pairup.i=NULL,
                             pairup.j=NULL, long.dim.are.variables=TRUE,
                             verbose=TRUE, identicalQs=FALSE, R.code.only=FALSE)
\S4method{qpNrr}{matrix}(data, q=1, nTests=100, alpha=0.05, pairup.i=NULL,
                         pairup.j=NULL, long.dim.are.variables=TRUE,
                         verbose=TRUE, identicalQs=FALSE, R.code.only=FALSE)
}
\arguments{
  \item{data}{data set from where to estimate the non-rejection rates.
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
  \item{identicalQs}{experimental feature.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
Note that the possible values of \code{q} should be in the range 1 to
\code{min(p,n-3)}, where \code{p} is the number of variables and \code{n}
the number of observations. The computational cost increases linearly with
\code{q} and quadratically in \code{p}.
}
\value{
A symmetric matrix of estimated non-rejection rates.
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
nVar <- 50 # number of variables
maxCon <- 5  # maximum connectivity per variable
nObs <- 30 # number of observations to simulate

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

X <- qpSampleMvnorm(K, nObs)

nrr.estimates <- qpNrr(X, q=5, verbose=FALSE)

summary(nrr.estimates[upper.tri(nrr.estimates) & I])

summary(nrr.estimates[upper.tri(nrr.estimates) & !I])

}
\keyword{models}
\keyword{multivariate}
