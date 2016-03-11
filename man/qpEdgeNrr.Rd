\name{qpEdgeNrr}
\alias{qpEdgeNrr}
\alias{qpEdgeNrr,ExpressionSet-method}
\alias{qpEdgeNrr,data.frame-method}
\alias{qpEdgeNrr,matrix-method}
\alias{qpEdgeNrr,SsdMatrix-method}
\alias{qpEdgeCor,UGgmm-method}
\alias{qpEdgeCor,matrix-method}
\alias{qpEdgeCor}
\alias{qpPathWeight,character-method}
\alias{qpPathWeight}

\title{
Non-rejection rate estimation for a pair of variables
}
\description{
Estimates the non-rejection rate for one pair of variables.
}
\usage{
\S4method{qpEdgeNrr}{ExpressionSet}(X, i=1, j=2, q=1, restrict.Q=NULL, fix.Q=NULL,
                                    nTests=100, alpha=0.05, exact.test=TRUE,
                                    use=c("complete.obs", "em"), tol=0.01,
                                    R.code.only=FALSE)
\S4method{qpEdgeNrr}{data.frame}(X, i=1, j=2, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                                 nTests=100, alpha=0.05, long.dim.are.variables=TRUE,
                                 exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                                 R.code.only=FALSE)
\S4method{qpEdgeNrr}{matrix}(X, i=1, j=2, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL,
                             nTests=100, alpha=0.05, long.dim.are.variables=TRUE,
                             exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                             R.code.only=FALSE)
\S4method{qpEdgeNrr}{SsdMatrix}(X, i=1, j=2, q=1, restrict.Q=NULL, fix.Q=NULL,
                                nTests=100, alpha=0.05, R.code.only=FALSE)
}
\arguments{
  \item{X}{data set from where the non-rejection rate should be estimated. It
       can be either an \code{ExpressionSet} object a
       data frame, a matrix or an \code{\link{SsdMatrix-class}} object. In the
       latter case, the input matrix should correspond to a sample covariance matrix
       of data from which we want to estimate the non-rejection rate for a pair of
       variables. The function \code{\link{qpCov}()} can be used to estimate such
       matrices.}
  \item{i}{index or name of one of the two variables in \code{X} to test.}
  \item{j}{index or name of the other variable in \code{X} to test.}
  \item{q}{order of the conditioning subsets employed in the calculation.}
  \item{I}{indexes or names of the variables in \code{X} that are discrete
       when \code{X} is a matrix or a data frame.}
  \item{restrict.Q}{indexes or names of the variables in \code{X} that
       restrict the sample space of conditioning subsets Q.}
  \item{fix.Q}{indexes or names of the variables in \code{X} that should be
       fixed within every conditioning conditioning subsets Q.}
  \item{nTests}{number of tests to perform for each pair for variables.}
  \item{alpha}{significance level of each test.}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when data are in a data frame or in a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
  \item{exact.test}{logical; if \code{FALSE} an asymptotic conditional independence
       test is employed with mixed (i.e., continuous and discrete) data;
       if \code{TRUE} (default) then an exact conditional independence test with
       mixed data is employed.See details below regarding this argument.}
  \item{use}{a character string defining the way in which calculations are done in the
       presence of missing values. It can be either \code{"complete.obs"} (default)
       or \code{"em"}.}
  \item{tol}{maximum tolerance controlling the convergence of the EM algorithm employed
             when the argument \code{use="em"}.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
The estimation of the non-rejection rate for a pair of variables is calculated
as the fraction of tests that accept the null hypothesis of conditional independence
given a set of randomly sampled q-order conditionals.

Note that the possible values of \code{q} should be in the range 1 to
\code{min(p,n-3)}, where \code{p} is the number of variables and \code{n}
the number of observations. The computational cost increases linearly with
\code{q}.

When \code{I} is set different to \code{NULL} then mixed graphical model theory
is employed and, concretely, it is assumed that the data comes from an homogeneous
conditional Gaussian distribution. In this setting further restrictions to the
maximum value of \code{q} apply, concretely, it cannot be smaller than
\code{p} plus the number of levels of the discrete variables involved in the
marginal distributions employed by the algorithm. By default, with
\code{exact.test=TRUE}, an exact test for conditional independence is employed,
otherwise an asymptotic one will be used. Full details on these features can
be found in Tur, Roverato and Castelo (2014).

The argument \code{I} specifying what variables are discrete actually applies only
when \code{X} is a matrix object since in the other cases data types are specified
for each data columns or slot.
}
\value{
An estimate of the non-rejection rate for the particular given pair of
variables.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Tur, I., Roverato, A. and Castelo, R. Mapping eQTL networks with mixed graphical models.
\emph{Submitted}, \url{http://arxiv.org/abs/1402.4547}, 2014.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
  \code{\link{qpAvgNrr}}
  \code{\link{qpHist}}
  \code{\link{qpGraphDensity}}
  \code{\link{qpClique}}
  \code{\link{qpCov}}
}
\examples{
require(mvtnorm)

nObs <- 100 ## number of observations to simulate

## the following adjacency matrix describes an undirected graph
## where vertex 3 is conditionally independent of 4 given 1 AND 2
A <- matrix(c(FALSE,  TRUE,  TRUE,  TRUE,
              TRUE,  FALSE,  TRUE,  TRUE,
              TRUE,   TRUE, FALSE, FALSE,
              TRUE,   TRUE, FALSE, FALSE), nrow=4, ncol=4, byrow=TRUE)
Sigma <- qpG2Sigma(A, rho=0.5)

X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

qpEdgeNrr(X, i=3, j=4, q=1, long.dim.are.variables=FALSE)

qpEdgeNrr(X, i=3, j=4, q=2, long.dim.are.variables=FALSE)
}
\keyword{models}
\keyword{multivariate}
