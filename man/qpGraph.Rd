\name{qpGraph}
\alias{qpGraph}

\title{
The qp-graph
}
\description{
Obtains a qp-graph from a matrix of non-rejection rates
}
\usage{
qpGraph(nrrMatrix, threshold=NULL, topPairs=NULL, pairup.i=NULL, pairup.j=NULL,
        return.type=c("adjacency.matrix", "edge.list", "graphNEL", "graphAM"))
}
\arguments{
  \item{nrrMatrix}{matrix of non-rejection rates.}
  \item{threshold}{threshold on the non-rejection rate above which pairs of
       variables are assumed to be disconnected in the resulting qp-graph.}
  \item{topPairs}{number of edges from the top of the ranking, defined by the
       non-rejection rates in \code{nrrMatrix}, to use to form the resulting
       qp-graph. This parameter is incompatible with a value different from
       \code{NULL} in \code{threshold}.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{return.type}{type of data structure on which the resulting undirected
       graph should be returned. Either a logical adjacency matrix with cells
       set to TRUE when the two indexing variables are connected in the qp-graph
       (default), or a list of edges in a matrix where each row corresponds to
       one edge and the two columns contain the two vertices defining each edge,
       or a \code{graphNEL-class} object, or a \code{graphAM-class} object.}
}
\details{
This function requires the \code{graph} package when \code{return.type=graphNEL}
or \code{return.type=graphAM}.
}
\value{
The resulting qp-graph as either an adjacency matrix, a \code{graphNEL} object or
a \code{graphAM} object, depending on the value of the \code{return.type}
parameter. Note that when some gold-standard graph is available for comparison,
a value for the parameter \code{threshold} can be found by calculating a
precision-recall curve with \code{qpPrecisionRecall} with respect to this
gold-standard, and then using \code{qpPRscoreThreshold}. Parameters
\code{threshold} and \code{topPairs} are mutually exclusive, that is, when
we specify with \code{topPairs=n} that we want a qp-graph with \code{n} edges
then \code{threshold} cannot be used.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
  \code{\link{qpAvgNrr}}
  \code{\link{qpEdgeNrr}}
  \code{\link{qpAnyGraph}}
  \code{\link{qpGraphDensity}}
  \code{\link{qpClique}}
  \code{\link{qpPrecisionRecall}}
  \code{\link{qpPRscoreThreshold}}
}
\examples{
require(mvtnorm)

nVar <- 50  ## number of variables
maxCon <- 5 ## maximum connectivity per variable
nObs <- 30  ## number of observations to simulate

set.seed(123)

A <- qpRndGraph(p=nVar, d=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=as.matrix(Sigma))

## estimate non-rejection rates
nrr.estimates <- qpNrr(X, q=5, verbose=FALSE)

## the higher the threshold
g <- qpGraph(nrr.estimates, threshold=0.9)

## the denser the qp-graph
(sum(g)/2) / (nVar*(nVar-1)/2)

## the lower the threshold
g <- qpGraph(nrr.estimates, threshold=0.5)

## the sparser the qp-graph
(sum(g)/2) / (nVar*(nVar-1)/2)

}
\keyword{models}
\keyword{multivariate}
