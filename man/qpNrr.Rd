\name{qpNrr}
\alias{qpNrr}
\alias{qpNrr,ExpressionSet-method}
\alias{qpNrr,cross-method}
\alias{qpNrr,data.frame-method}
\alias{qpNrr,matrix-method}

\title{
Non-rejection rate estimation
}
\description{
Estimates non-rejection rates for every pair of variables.
}
\usage{
\S4method{qpNrr}{ExpressionSet}(X, q=1, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                                alpha=0.05, pairup.i=NULL, pairup.j=NULL,
                                verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                                use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE,
                                clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10)
\S4method{qpNrr}{cross}(X, q=1, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                         alpha=0.05, pairup.i=NULL, pairup.j=NULL, verbose=TRUE,
                         identicalQs=TRUE, exact.test=TRUE, use=c("complete.obs", "em"),
                         tol=0.01, R.code.only=FALSE, clusterSize=1, estimateTime=FALSE,
                         nAdj2estimateTime=10)
\S4method{qpNrr}{data.frame}(X, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                             alpha=0.05, pairup.i=NULL, pairup.j=NULL,
                             long.dim.are.variables=TRUE, verbose=TRUE,
                             identicalQs=TRUE, exact.test=TRUE, use=c("complete.obs", "em"),
                             tol=0.01, R.code.only=FALSE, clusterSize=1,
                             estimateTime=FALSE, nAdj2estimateTime=10)
\S4method{qpNrr}{matrix}(X, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=100,
                         alpha=0.05, pairup.i=NULL, pairup.j=NULL,
                         long.dim.are.variables=TRUE, verbose=TRUE, identicalQs=TRUE,
                         exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                         R.code.only=FALSE, clusterSize=1, estimateTime=FALSE,
                         nAdj2estimateTime=10)
}
\arguments{
  \item{X}{data set from where to estimate the non-rejection rates.
       It can be an \code{ExpressionSet} object, a \code{qtl/cross} object,
       a \code{data.frame} object or a \code{matrix} object.}
  \item{q}{partial-correlation order to be employed.}
  \item{I}{indexes or names of the variables in \code{X} that are discrete.
       See details below regarding this argument.}
  \item{restrict.Q}{indexes or names of the variables in \code{X} that
       restrict the sample space of conditioning subsets Q.}
  \item{fix.Q}{indexes or names of the variables in \code{X} that should be
       fixed within every conditioning conditioning subsets Q.}
  \item{nTests}{number of tests to perform for each pair for variables.}
  \item{alpha}{significance level of each test.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{long.dim.are.variables}{logical; if \code{TRUE} it is assumed
       that when data are in a data frame or in a matrix, the longer dimension
       is the one defining the random variables (default); if \code{FALSE}, then
       random variables are assumed to be at the columns of the data frame or matrix.}
  \item{verbose}{show progress on the calculations.}
  \item{identicalQs}{use identical conditioning subsets for every pair of vertices
       (default), otherwise sample a new collection of \code{nTests} subsets for
       each pair of vertices.}
  \item{exact.test}{logical; if \code{FALSE} an asymptotic conditional independence
       test is employed with mixed (i.e., continuous and discrete) data;
       if \code{TRUE} (default) then an exact conditional independence test with
       mixed data is employed. See details below regarding this argument.}
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
Note that for pure continuous data the possible values of \code{q} should be in the
range 1 to \code{min(p, n-3)}, where \code{p} is the number of variables and
\code{n} the number of observations. The computational cost increases linearly
with \code{q} and quadratically in \code{p}. When setting \code{identicalQs}
to \code{FALSE} the computational cost may increase between 2 times and one
order of magnitude (depending on \code{p} and \code{q}) while asymptotically
the estimation of the non-rejection rate converges to the same value. Full
details on the calculation of the non-rejection rate can be found in
Castelo and Roverato (2006).

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

In the case that \code{X} is a \code{qtl/cross} object, the default \code{NULL}
values in arguments \code{pairup.i} and \code{pairup.j} actually imply pairing
all markers and phenotypes with numerical phenotypes only (including integer phenotypes).
Likewise, the default argument \code{restrict.Q=NULL} implies setting \code{restrict.Q}
to all numeric phenotypes. Setting these arguments to values other than \code{NULL}
allows the user to use those particular values being set.
}
\value{
A \code{\link{dspMatrix-class}} symmetric matrix of estimated non-rejection
rates with the diagonal set to \code{NA} values. If arguments \code{pairup.i}
and \code{pairup.j} are employed, those cells outside the constrained pairs
will get also a \code{NA} value.

Note, however, that when \code{estimateTime=TRUE}, then instead of the matrix
of estimated non-rejection rates, a vector specifying the estimated number of
days, hours, minutes and seconds for completion of the calculations is returned.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Tur, I., Roverato, A. and Castelo, R. Mapping eQTL networks with mixed graphical Markov models.
\emph{Genetics}, 198(4):1377-1393, 2014.
}
\author{R. Castelo, A. Roverato and I. Tur}
\seealso{
  \code{\link{qpAvgNrr}}
  \code{\link{qpEdgeNrr}}
  \code{\link{qpHist}}
  \code{\link{qpGraphDensity}}
  \code{\link{qpClique}}
}
\examples{
nVar <- 50  ## number of variables
maxCon <- 3 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

## simulate an undirected Gaussian graphical model
## determined by some random undirected d-regular graph
model <- rUGgmm(dRegularGraphParam(p=nVar, d=maxCon), rho=0.5)
model

## simulate data from this model
X <- rmvnorm(nObs, model)
dim(X)

## estimate non-rejection rates with q=3
nrr.estimates <- qpNrr(X, q=3, verbose=FALSE)

## create an adjacency matrix of the undirected graph
## determining the undirected Gaussian graphical model
A <- as(model$g, "matrix") == 1

## distribution of non-rejection rates for the present edges
summary(nrr.estimates[upper.tri(nrr.estimates) & A])

## distribution of non-rejection rates for the missing edges
summary(nrr.estimates[upper.tri(nrr.estimates) & !A])

\dontrun{
## using R code only this would take much more time
qpNrr(X, q=3, R.code.only=TRUE, estimateTime=TRUE)

## only for moderate and large numbers of variables the
## use of a cluster of processors speeds up the calculations

library(snow)
library(rlecuyer)

nVar <- 500
maxCon <- 3
model <- rUGgmm(dRegularGraphParam(p=nVar, d=maxCon), rho=0.5)
X <- rmvnorm(nObs, model)

system.time(nrr.estimates <- qpNrr(X, q=10, verbose=TRUE))
system.time(nrr.estimates <- qpNrr(X, q=10, verbose=TRUE, clusterSize=4))
}
}
\keyword{models}
\keyword{multivariate}
