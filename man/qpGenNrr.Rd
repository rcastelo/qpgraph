\name{qpGenNrr}
\alias{qpGenNrr}
\alias{qpGenNrr,ExpressionSet-method}
\alias{qpGenNrr,data.frame-method}
\alias{qpGenNrr,matrix-method}
\alias{qpGenNrr,list-method}

\title{
Generalized non-rejection rate estimation
}
\description{
Estimates generalized non-rejection rates for every pair of variables from
two or more data sets.
}
\usage{
\S4method{qpGenNrr}{ExpressionSet}(X, datasetIdx=1, qOrders=NULL, I=NULL, restrict.Q=NULL,
                                   fix.Q=NULL, return.all=FALSE, nTests=100, alpha=0.05,
                                   pairup.i=NULL, pairup.j=NULL, verbose=TRUE, identicalQs=TRUE,
                                   exact.test=TRUE, use=c("complete.obs", "em"), tol=0.01,
                                   R.code.only=FALSE, clusterSize=1, estimateTime=FALSE,
                                   nAdj2estimateTime=10)
\S4method{qpGenNrr}{data.frame}(X, datasetIdx=1, qOrders=NULL, I=NULL, restrict.Q=NULL,
                                fix.Q=NULL, return.all=FALSE, nTests=100, alpha=0.05,
                                pairup.i=NULL, pairup.j=NULL, long.dim.are.variables=TRUE,
                                verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                                use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE,
                                clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10)
\S4method{qpGenNrr}{matrix}(X, datasetIdx=1, qOrders=NULL, I=NULL, restrict.Q=NULL,
                            fix.Q=NULL, return.all=FALSE, nTests=100, alpha=0.05,
                            pairup.i=NULL, pairup.j=NULL, long.dim.are.variables=TRUE,
                            verbose=TRUE, identicalQs=TRUE, exact.test=TRUE,
                            use=c("complete.obs", "em"), tol=0.01, R.code.only=FALSE,
                            clusterSize=1, estimateTime=FALSE, nAdj2estimateTime=10)
}
\arguments{
  \item{X}{data set from where to estimate the average non-rejection rates.
       It can be an ExpressionSet object, a data frame or a matrix.}
  \item{datasetIdx}{either a single number, or a character string, indicating
       the column in the phenotypic data of the \code{ExpressionSet} object,
       or in the input matrix or data frame, containing the indexes to the
       data sets. Alternatively, it can be a vector of these indexes with
       as many positions as samples.}
  \item{qOrders}{either a NULL value (default) indicating that a default guess
       on the q-order will be employed for each data set or a vector of particular
       orders with one for each data set. The default guess corresponds to
       the floor of the median value among the valid q orders of the data set.}
  \item{I}{indexes or names of the variables in \code{X} that are discrete.
       When \code{X} is an \code{ExpressionSet} then \code{I} may contain
       only names of the phenotypic variables in \code{X}. See details below
       regarding this argument.}
  \item{restrict.Q}{indexes or names of the variables in \code{X} that
       restrict the sample space of conditioning subsets Q.}
  \item{fix.Q}{indexes or names of the variables in \code{X} that should be
       fixed within every conditioning conditioning subsets Q.}
  \item{return.all}{logical; if TRUE all intervining non-rejection rates will be
       return in a matrix per dataset within a list; FALSE (default) if only
       generalized non-rejection rates should be returned.}
  \item{nTests}{number of tests to perform for each pair for variables.}
  \item{alpha}{significance level of each test.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{long.dim.are.variables}{logical; if TRUE it is assumed
       that when the data is a data frame or a matrix, the longer dimension
       is the one defining the random variables; if FALSE, then random
       variables are assumed to be at the columns of the data frame or matrix.}
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
  \item{R.code.only}{logical; if FALSE then the faster C implementation is used
       (default); if TRUE then only R code is executed.}
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
should be in the range 1 to \code{min(p,n-3)}, where \code{p} is the number of
variables and \code{n} the number of observations for the corresponding data set.
The computational cost increases linearly within each \code{q} value and
quadratically in \code{p}.  When setting \code{identicalQs} to \code{FALSE} the
computational cost may increase between 2 times and one order of magnitude
(depending on \code{p} and \code{q}) while asymptotically the estimation of the
non-rejection rate converges to the same value.

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
A list containing the following two or more entries: a first one with name
\code{genNrr} with a \code{\link{dspMatrix-class}} symmetric matrix of estimated
generalized non-rejection rates with the diagonal set to \code{NA} values. When
using the arguments \code{pairup.i} and \code{pairup.j}, those cells outside the
constraint pairs will get also a \code{NA} value; a second one with name
\code{qOrders} with the q-orders employed in the calculation for each data set;
if \code{return.all=TRUE} then there will be one additional entry for each data
set containing the matrix of the non-rejection rates estimated from that data
set with the corresponding q-order, using the indexing value of the data set as
entry name.

Note, however, that when \code{estimateTime=TRUE}, then instead of the list with
matrices of estimated (generalized) non-rejection rates, a vector specifying the
estimated number of days, hours, minutes and seconds for completion of the
calculations is returned.
}
\references{
Castelo, R. and Roverato, A. Reverse engineering molecular regulatory
networks from microarray data with qp-graphs. \emph{J. Comp. Biol.},
16(2):213-227, 2009.

Tur, I., Roverato, A. and Castelo, R. Mapping eQTL networks with mixed graphical Markov models.
\emph{Genetics}, 198:1377-1393, 2014.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
  \code{\link{qpAvgNrr}}
  \code{\link{qpEdgeNrr}}
  \code{\link{qpHist}}
  \code{\link{qpGraphDensity}}
  \code{\link{qpClique}}
}
\examples{
nVar <- 50  ## number of variables
maxCon <- 5 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

## simulate two independent Gaussian graphical models determined
## by two undirected d-regular graphs
model1 <- rUGgmm(dRegularGraphParam(p=nVar, d=maxCon), rho=0.5)
model2 <- rUGgmm(dRegularGraphParam(p=nVar, d=maxCon), rho=0.5)

## simulate two independent data sets from the previous graphical models
X1 <- rmvnorm(nObs, model1)
dim(X1)
X2 <- rmvnorm(nObs, model2)
dim(X2)

## estimate generalized non-rejection rates from the joint data
nrr.estimates <- qpGenNrr(rbind(X1, X2),
                          datasetIdx=rep(1:2, each=nObs),
                          qOrders=c("1"=5, "2"=5),
                          long.dim.are.variables=FALSE, verbose=FALSE)

## create adjacency matrices from the undirected graphs
## determining the two Gaussian graphical models
A1 <- as(model1$g, "matrix") == 1
A2 <- as(model2$g, "matrix") == 1

## distribution of generalized non-rejection rates for the common present edges
summary(nrr.estimates$genNrr[upper.tri(nrr.estimates$genNrr) & A1 & A2])

## distribution of generalized non-rejection rates for the present edges specific to A1
summary(nrr.estimates$genNrr[upper.tri(nrr.estimates$genNrr) & A1 & !A2])

## distribution of generalized non-rejection rates for the present edges specific to A2
summary(nrr.estimates$genNrr[upper.tri(nrr.estimates$genNrr) & !A1 & A2])

## distribution of generalized non-rejection rates for the common missing edges
summary(nrr.estimates$genNrr[upper.tri(nrr.estimates$genNrr) & !A1 & !A2])

## compare with the average non-rejection rate on the pooled data set
avgnrr.estimates <- qpNrr(rbind(X1, X2), q=5, long.dim.are.variables=FALSE, verbose=FALSE)

## distribution of average non-rejection rates for the common present edges
summary(avgnrr.estimates[upper.tri(avgnrr.estimates) & A1 & A2])

## distribution of average non-rejection rates for the present edges specific to A1
summary(avgnrr.estimates[upper.tri(avgnrr.estimates) & A1 & !A2])

## distribution of average non-rejection rates for the present edges specific to A2
summary(avgnrr.estimates[upper.tri(avgnrr.estimates) & !A1 & A2])

## distribution of average non-rejection rates for the common missing edges
summary(avgnrr.estimates[upper.tri(avgnrr.estimates) & !A1 & !A2])
}
\keyword{models}
\keyword{multivariate}
