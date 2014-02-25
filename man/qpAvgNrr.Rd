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
\S4method{qpAvgNrr}{ExpressionSet}(X, qOrders=4, I=NULL, restrict.Q=NULL,
                                   fix.Q=NULL, nTests=100, alpha=0.05,
                                   pairup.i=NULL, pairup.j=NULL, type=c("arith.mean"),
                                   verbose=TRUE, identicalQs=TRUE,
                                   exact.test=TRUE, use=c("complete.obs", "em"),
                                   tol=0.01, R.code.only=FALSE, clusterSize=1,
                                   estimateTime=FALSE, nAdj2estimateTime=10)
\S4method{qpAvgNrr}{data.frame}(X, qOrders=4, I=NULL, restrict.Q=NULL,
                                fix.Q=NULL, nTests=100, alpha=0.05, pairup.i=NULL,
                                pairup.j=NULL, long.dim.are.variables=TRUE,
                                type=c("arith.mean"), verbose=TRUE,
                                identicalQs=TRUE, exact.test=TRUE,
                                use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE,
                                clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10)
\S4method{qpAvgNrr}{matrix}(X, qOrders=4, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                            nTests=100, alpha=0.05, pairup.i=NULL,
                            pairup.j=NULL, long.dim.are.variables=TRUE,
                            type=c("arith.mean"), verbose=TRUE,
                            identicalQs=TRUE, exact.test=TRUE,
                            use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE,
                            clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10)
}
\arguments{
  \item{X}{data set from where to estimate the average non-rejection rates.
       It can be an ExpressionSet object, a data frame or a matrix.}
  \item{qOrders}{either a number of partial-correlation orders or a vector of
       vector of particular orders to be employed in the calculation.}
  \item{I}{indexes or names of the variables in \code{X} that are discrete.
       When \code{X} is an \code{ExpressionSet} then \code{I} may contain
       only names of the phenotypic variables in \code{X}. See details below
       regarding this argument.}
  \item{restrict.Q}{indexes or names of the variables in \code{X} that
       restrict the sample space of conditioning subsets Q.}
  \item{fix.Q}{indexes or names of the variables in \code{X} that should be
       fixed within every conditioning conditioning subsets Q.}
  \item{nTests}{number of tests to perform for each pair for variables.}
  \item{alpha}{significance level of each test.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{long.dim.are.variables}{logical; if \code{TRUE} it is assumed
       that when the data is a data frame or a matrix, the longer dimension
       is the one defining the random variables; if \code{FALSE}, then random
       variables are assumed to be at the columns of the data frame or matrix.}
  \item{type}{type of average. By now only the arithmetic mean is available.}
  \item{verbose}{show progress on the calculations.}
  \item{identicalQs}{use identical conditioning subsets for every pair of vertices
       (default), otherwise sample a new collection of \code{nTests} subsets for
       each pair of vertices.}
  \item{exact.test}{logical; if \code{FALSE} an asymptotic conditional independence
       test is employed with mixed (i.e., continuous and discrete) data;
       if \code{TRUE} (default) then an exact conditional independence test with
       mixed data is employed.}
  \item{use}{a character string defining the way in which calculations are done in the
       presence of missing values. It can be either \code{"complete.obs"} (default)
       or \code{"em"}.}
  \item{tol}{maximum tolerance controlling the convergence of the EM algorithm employed
             when the argument \code{use="em"}.}
  \item{R.code.only}{logical; if \code{FALSE} then the faster C implementation is used
       (default); if \code{TRUE} then only R code is executed.}
  \item{clusterSize}{size of the cluster of processors to employ if we wish to
       speed-up the calculations by performing them in parallel. A value of 1
       (default) implies a single-processor execution. The use of a cluster of
       processors requires having previously loaded the packages \code{snow}
       and \code{rlecuyer}.}
  \item{estimateTime}{logical; if \code{TRUE} then the time for carrying out the
       calculations with the given parameters is estimated by calculating for a
       limited number of adjacencies, specified by \code{nAdj2estimateTime}, and
       extrapolating the elapsed time; if \code{FALSE} (default) calculations are
       performed normally till they finish.}
  \item{nAdj2estimateTime}{number of adjacencies to employ when estimating the
       time of calculations (\code{estimateTime=TRUE}). By default this has a
       default value of 10 adjacencies and larger values should provide more
       accurate estimates. This might be relevant when using a cluster facility.}
}
\details{
Note that when specifying a vector of particular orders \code{q}, these values
should be in the range 1 to \code{min(p, n-3)}, where \code{p} is the number of
variables and \code{n} the number of observations. The computational cost
increases linearly within each \code{q} value and quadratically in \code{p}.
When setting \code{identicalQs} to \code{FALSE} the computational cost may
increase between 2 times and one order of magnitude (depending on \code{p} and
\code{q}) while asymptotically the estimation of the non-rejection rate
converges to the same value.

When \code{I} is set different to \code{NULL} then mixed graphical model theory
is employed and, concretely, it is assumed that the data comes from an homogeneous
conditional Gaussian distribution. In this setting further restrictions to the
maximum value of \code{q} apply, concretely, it cannot be smaller than
\code{p} plus the number of levels of the discrete variables involved in the
marginal distributions employed by the algorithm. By default, with
\code{exact.test=TRUE}, an exact test for conditional independence is employed,
otherwise an asymptotic one will be used. Full details on these features can
be found in Tur, Roverato and Castelo (2014).
}
\value{
A \code{\link{dspMatrix-class}} symmetric matrix of estimated average
non-rejection rates with the diagonal set to \code{NA} values. When using the
arguments \code{pairup.i} and \code{pairup.j}, those cells outside the
constraint pairs will get also a \code{NA} value.

Note, however, that when \code{estimateTime=TRUE}, then instead of the matrix
of estimated average non-rejection rates, a vector specifying the estimated
number of days, hours, minutes and seconds for completion of the calculations
is returned.
}
\references{
Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
networks from microarray data with qp-graphs. \emph{J. Comp. Biol.},
16(2):213-227, 2009.

Tur, I., Roverato, A. and Castelo, R. Mapping eQTL networks with mixed graphical models.
\emph{Submitted}, \url{http://arxiv.org/abs/1402.4547}, 2014.
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
maxCon <- 3 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

A <- qpRndGraph(p=nVar, d=maxCon)
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
maxCon <- 3
A <- qpRndGraph(p=nVar, d=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

system.time(avgnrr.estimates <- qpAvgNrr(X, q=10, verbose=TRUE))
system.time(avgnrr.estimates <- qpAvgNrr(X, q=10, verbose=TRUE, clusterSize=4))
}
}
\keyword{models}
\keyword{multivariate}
