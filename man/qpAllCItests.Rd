\name{qpAllCItests}
\alias{qpAllCItests}
\alias{qpAllCItests,matrix-method}

\title{
Tests of conditional independence
}
\description{
Performs a test of conditional independence for every pair of variables.
}
\usage{
\S4method{qpAllCItests}{matrix}(X, I=NULL, Q=NULL, pairup.i=NULL, pairup.j=NULL,
                                long.dim.are.variables=TRUE, exact.test=TRUE,
                                use=c("complete.obs", "em"), tol=0.01,
                                return.type=c("p.value", "statn", "all"), verbose=TRUE,
                                R.code.only=FALSE, clusterSize=1, estimateTime=FALSE,
                                nAdj2estimateTime=10)
}
\arguments{
  \item{X}{data set from where to estimate the non-rejection rates.
       It can be an ExpressionSet object, a data frame or a matrix.}
  \item{I}{indexes or names of the variables in \code{X} that are discrete.
       See details below regarding this argument.}
  \item{Q}{indexes or names of the variables in \code{X} forming the conditioning set.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{long.dim.are.variables}{logical; if \code{TRUE} it is assumed
       that when data are in a data frame or in a matrix, the longer dimension
       is the one defining the random variables (default); if \code{FALSE}, then
       random variables are assumed to be at the columns of the data frame or matrix.}
  \item{exact.test}{logical; if \code{FALSE} an asymptotic conditional independence
       test is employed with mixed (i.e., continuous and discrete) data;
       if \code{TRUE} (default) then an exact conditional independence test with
       mixed data is employed. See details below regarding this argument.}
  \item{use}{a character string defining the way in which calculations are done in the
       presence of missing values. It can be either \code{"complete.obs"} (default)
       or \code{"em"}.}
  \item{tol}{maximum tolerance controlling the convergence of the EM algorithm employed
             when the argument \code{use="em"}.}
  \item{return.type}{type of value returned by this function. By default \code{"p.value"}
       indicates that a list containing a matrix of p-values from all performed conditional
       independence (CI) tests will be returned. If \code{return.type="statn"} then a list
       containing the matrix of the statistics and the sample sizes on each CI test, will be
       returned. If \code{return.type="all"} then all previous matrices of values will be
       returned within a list.}
  \item{verbose}{show progress on the calculations.}
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
When \code{I} is set different to \code{NULL} then mixed graphical model theory
is employed and, concretely, it is assumed that the data comes from an homogeneous
conditional Gaussian distribution. By default, with \code{exact.test=TRUE}, an
exact test for conditional independence is employed, otherwise an asymptotic one
will be used. Full details on these features can be found in Tur, Roverato and Castelo (2014).
}
\value{
A list with three entries called \code{p.value}, \code{statistic} and \code{n}
corresponding to a \code{\link{dspMatrix-class}} symmetric matrix of p-values for the null
hypothesis of coindtional independence with the diagonal set to \code{NA} values,
an analogous matrix of the statistics of each test and of the sample sizes, respectively.
These returned values, however, depend on the setting of argument \code{return.type} which,
by default, enables only returning the matrix of p-values.
If arguments \code{pairup.i} and \code{pairup.j} are employed, those cells outside
the constrained pairs will get also a \code{NA} value.

Note, however, that when \code{estimateTime=TRUE}, then instead of the matrix
of estimated non-rejection rates, a vector specifying the estimated number of
days, hours, minutes and seconds for completion of the calculations is returned.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Tur, I., Roverato, A. and Castelo, R. Mapping eQTL networks with mixed graphical models.
\emph{Submitted}, \url{http://arxiv.org/abs/1402.4547}, 2014.
}
\author{R. Castelo, A. Roverato and I. Tur}
\seealso{
  \code{\link{qpCItest}}
}
\examples{
library(mvtnorm)

nVar <- 50  ## number of variables
maxCon <- 3 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

A <- qpRndGraph(p=nVar, d=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

alltests <- qpAllCItests(X, verbose=FALSE)

## distribution of p-values for the present edges
summary(alltests$p.value[upper.tri(alltests$p.value) & A])

## distribution of p-values for the missing edges
summary(alltests$p.value[upper.tri(alltests$p.value) & !A])

}
\keyword{models}
\keyword{multivariate}
