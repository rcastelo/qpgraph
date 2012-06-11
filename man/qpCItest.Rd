\name{qpCItest}
\alias{qpCItest}
\alias{qpCItest,smlSet-method}
\alias{qpCItest,ExpressionSet-method}
\alias{qpCItest,data.frame-method}
\alias{qpCItest,matrix-method}

\title{
Conditional independence test
}
\description{
Performs a conditional independence test between two variables given
a conditioning set.
}
\usage{
\S4method{qpCItest}{smlSet}(X, i=1, j=2, Q=c(), exact.test=TRUE, use=c("complete.obs", "em"),
                            tol=0.01, R.code.only=FALSE)
\S4method{qpCItest}{ExpressionSet}(X, i=1, j=2, Q=c(), exact.test=TRUE, use=c("complete.obs", "em"),
                                   tol=0.01, R.code.only=FALSE)
\S4method{qpCItest}{data.frame}(X, i=1, j=2, Q=c(), I=NULL, long.dim.are.variables=TRUE,
                                exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE)
\S4method{qpCItest}{matrix}(X, i=1, j=2, Q=c(), I=NULL, n=NULL, long.dim.are.variables=TRUE,
                            exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE)
}
\arguments{
  \item{X}{data set where the test should be performed. It can be either an
       \code{ExpressionSet} object, a data frame, or a matrix. If it is a matrix
       and the matrix is squared then this function assumes the matrix corresponds
       to the sample covariance matrix of the data and the sample size parameter
       \code{n} should be provided.}
  \item{i}{index or name of one of the two variables in \code{X} to test.}
  \item{j}{index or name of the other variable in \code{X} to test.}
  \item{Q}{indexes or names of the variables in \code{X} forming the conditioning set.}
  \item{I}{indexes or names of the variables in \code{X} that are discrete. See details
       below regarding this argument.}
  \item{n}{number of observations in the data set. Only necessary when the
       sample covariance matrix is provided through the \code{X} parameter.}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when data are in a data frame or in a matrix, the longer dimension
       is the one defining the random variables (default); if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
  \item{exact.test}{logical; if \code{FALSE} an asymptotic likelihood ratio test of
       conditional independence test is employed with mixed (i.e., continuous and discrete)
       data; if \code{TRUE} (default) then an exact likelihood ratio test of conditional
       independence with mixed data is employed. See details below regarding this argument.}
  \item{use}{a character string defining the way in which calculations are done in the
       presence of missing values. It can be either \code{"complete.obs"} (default)
       or \code{"em"}.}
  \item{tol}{maximum tolerance controlling the convergence of the EM algorithm employed
             when the argument \code{use="em"}.}
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
}
\details{
When variables in \code{i, j} and \code{Q} are continuous this function performs a
conditional independence test using a t test for zero partial regression coefficient
(Lauritzen, 1996, pg. 150). Note that the size of possible \code{Q} sets should be in
the range 1 to \code{min(p,n-3)}, where \code{p} is the number of variables and \code{n}
the number of observations. The computational cost increases linearly with
the number of variables in \code{Q}.

When variables in \code{i, j} and \code{Q} are continuous and discrete (mixed data),
indicated with the \code{I} argument when \code{X} is a matrix, then mixed graphical
model theory (Lauritzen and Wermuth, 1989) is employed and, concretely, it is assumed
that data come from an homogeneous conditional Gaussian distribution. By default, with
\code{exact.test=TRUE}, an exact likelihood ratio test for conditional independence is
performed (Lauritzen, 1996, pg. 192-194; Tur and Castelo, 2011), otherwise an asymptotic
one is used.

In this setting further restrictions to the maximum value of \code{q} apply, concretely,
it cannot be smaller than \code{p} plus the number of levels of the discrete variables
involved in the marginal distributions employed by the algorithm.
}
\value{
A list with two members, the value of the statistic and its corresponding
P-value of rejecting the null hypothesis of conditional independence.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Lauritzen, S.L. \emph{Graphical models}. Oxford University Press, 1996.

Lauritzen, S.L and Wermuth, N. Graphical Models for associations between variables,
some of which are qualitative and some quantitative. \emph{Ann. Stat.}, 17(1):31-57, 1989.

Tur, I. and Castelo, R. Learning mixed graphical models from data with p larger than n,
In \emph{Proc. 27th Conference on Uncertainty in Artificial Intelligence},
F.G. Cozman and A. Pfeffer eds., pp. 689-697, AUAI Press, ISBN 978-0-9749039-7-2,
Barcelona, 2011.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
  \code{\link{qpEdgeNrr}}
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

qpCItest(X, i=3, j=4, Q=1, long.dim.are.variables=FALSE)

qpCItest(X, i=3, j=4, Q=c(1,2), long.dim.are.variables=FALSE)
}
\keyword{models}
\keyword{multivariate}
