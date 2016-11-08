\name{qpCItest}
\alias{qpCItest}
\alias{qpCItest,ExpressionSet-method}
\alias{qpCItest,cross-method}
\alias{qpCItest,data.frame-method}
\alias{qpCItest,matrix-method}
\alias{qpCItest,SsdMatrix-method}

\title{
Conditional independence test
}
\description{
Performs a conditional independence test between two variables given
a conditioning set.
}
\usage{
\S4method{qpCItest}{ExpressionSet}(X, i=1, j=2, Q=c(), exact.test=TRUE, use=c("complete.obs", "em"),
                                   tol=0.01, R.code.only=FALSE)
\S4method{qpCItest}{cross}(X, i=1, j=2, Q=c(), exact.test=TRUE, use=c("complete.obs", "em"),
                           tol=0.01, R.code.only=FALSE)
\S4method{qpCItest}{data.frame}(X, i=1, j=2, Q=c(), I=NULL, long.dim.are.variables=TRUE,
                                exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE)
\S4method{qpCItest}{matrix}(X, i=1, j=2, Q=c(), I=NULL, long.dim.are.variables=TRUE,
                            exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE)
\S4method{qpCItest}{SsdMatrix}(X, i=1, j=2, Q=c(), R.code.only=FALSE)
}
\arguments{
  \item{X}{data set where the test should be performed. It can be either
       an \code{ExpressionSet} object, a \code{qtl::cross} object,
       a data frame, a matrix or an \code{\link{SsdMatrix-class}} object. In the latter case,
       the input matrix should correspond to a sample covariance matrix of data on which we
       want to test for conditional independence. The function \code{\link{qpCov}()}
       can be used to estimate such matrices.}
  \item{i}{index or name of one of the two variables in \code{X} to test.}
  \item{j}{index or name of the other variable in \code{X} to test.}
  \item{Q}{indexes or names of the variables in \code{X} forming the conditioning set.}
  \item{I}{indexes or names of the variables in \code{X} that are discrete. See details
       below regarding this argument.}
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
When variables in \code{i, j} and \code{Q} are continuous and \code{I=NULL}, this function
performs a conditional independence test using a t-test for zero partial regression coefficient
(Lauritzen, 1996, pg. 150). Note that the size of possible \code{Q} sets should be in
the range 1 to \code{min(p,n-3)}, where \code{p} is the number of variables and \code{n}
the number of observations. The computational cost increases linearly with
the number of variables in \code{Q}.

When variables in \code{i, j} and \code{Q} are continuous and discrete (mixed data),
indicated with the \code{I} argument when \code{X} is a matrix, then mixed graphical
model theory (Lauritzen and Wermuth, 1989) is employed and, concretely, it is assumed
that data come from an homogeneous conditional Gaussian distribution. By default, with
\code{exact.test=TRUE}, an exact likelihood ratio test for conditional independence is
performed (Lauritzen, 1996, pg. 192-194; Tur, Roverato and Castelo, 2014), otherwise an
asymptotic one is used.

In this setting further restrictions to the maximum value of \code{q} apply, concretely,
it cannot be smaller than \code{p} plus the number of levels of the discrete variables
involved in the marginal distributions employed by the algorithm.
}
\value{
A list with class \code{"htest"} containing the following components:
\item{statistic}{in case of pure continuous data and \code{I=NULL}, the t-statistic for
  zero partial regression coefficient; when \code{I!=NULL}, the value \code{Lambda} of
  the likelihood ratio if \code{exact.test=TRUE} and \code{-n log Lambda} otherwise.}
\item{parameter}{in case of pure continuous data and \code{I=NULL}, the degrees of freedom
  for the t-statistic (\code{n-q-2}); when \code{I!=NULL}, the degrees of freedom for
  \code{-n log Lambda} of a chi-square distribution under the null hypothesis if
  \code{exact.test=FALSE} and the \code{(a, b)} parameters of a beta distribution under
  the null if \code{exact.test=TRUE}.}
\item{p.value}{the p-value for the test.}
\item{estimate}{in case of pure continuous data (\code{I=NULL}), the estimated partial
  regression coefficient. In case of mixed continuous and discrete data with \code{I!=NULL},
  the estimated partial eta-squared: the fraction of variance from \code{i} or \code{j}
  explained by the other tested variable after excluding the variance explained by the
  variables in \code{Q}. If one of the tested variables \code{i} or \code{j} is discrete,
  then the partial eta-squared is calculated on the tested continuous variable. If both,
  \code{i} and \code{j} are continuous, then the partial eta-squared is calculated on
  variable \code{i}.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of conditional independence test was
  performed.}
\item{data.name}{a character string giving the name(s) of the random variables involved in
  the conditional independence test.}
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Lauritzen, S.L. \emph{Graphical models}. Oxford University Press, 1996.

Lauritzen, S.L and Wermuth, N. Graphical Models for associations between variables,
some of which are qualitative and some quantitative. \emph{Ann. Stat.}, 17(1):31-57, 1989.

Tur, I., Roverato, A. and Castelo, R. Mapping eQTL networks with mixed graphical Markov models.
\emph{Genetics}, 198:1377-1393, 2014.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpCov}}
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
