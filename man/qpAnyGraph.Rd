\name{qpAnyGraph}
\alias{qpAnyGraph}

\title{
A graph
}
\description{
Obtains an undirected graph from a matrix of pairwise measurements
}
\usage{
qpAnyGraph(measurementsMatrix, threshold=NULL, remove=c("below", "above"),
           topPairs=NULL, decreasing=TRUE, pairup.i=NULL, pairup.j=NULL,
           return.type=c("adjacency.matrix", "edge.list", "graphNEL", "graphAM"))
}
\arguments{
  \item{measurementsMatrix}{matrix of pairwise measurements.}
  \item{threshold}{threshold on the measurements below or above which pairs of
       variables are assumed to be disconnected in the resulting graph.}
  \item{remove}{direction of the removal with the threshold. It should be
       either \code{"below"} (default) or \code{"above"}.}
  \item{topPairs}{number of edges from the top of the ranking, defined by the
       pairwise measurements in \code{measurementsMatrix}, to use to form the
       resulting graph. This parameter is incompatible with a value different
       from \code{NULL} in \code{threshold}.}
  \item{decreasing}{logical, only applies when topPairs is set; if \code{TRUE}
       then the ranking is made in decreasing order; if \code{FALSE} then is made
       in increasing order.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}}
  \item{return.type}{type of data structure on which the resulting undirected
       graph should be returned. Either a logical adjacency matrix with cells
       set to TRUE when the two indexing variables are connected in the graph
       (default), or a list of edges in a matrix where each row corresponds to
       one edge and the two columns contain the two vertices defining each edge,
       or a \code{graphNEL-class} object, or a \code{graphAM-class} object.}
}
\details{
This function requires the \code{graph} package when \code{return.type=graphNEL}
or \code{return.type=graphAM}.
}
\value{
The resulting undirected graph as either an adjacency matrix, a \code{graphNEL}
object or a \code{graphAM} object, depending on the value of the \code{return.type}
parameter. Note that when some gold-standard graph is available for comparison,
a value for the parameter \code{threshold} can be found by calculating a
precision-recall curve with \code{qpPrecisionRecall} with respect to this
gold-standard, and then using \code{qpPRscoreThreshold}. Parameters
\code{threshold} and \code{topPairs} are mutually exclusive, that is, when
we specify with \code{topPairs=n} that we want a graph with \code{n} edges
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
  \code{\link{qpGraph}}
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

A <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
Sigma <- qpG2Sigma(A, rho=0.5)
X <- rmvnorm(nObs, sigma=Sigma)

## estimate Pearson correlations
pcc.estimates <- qpPCC(X)

## the higher the threshold
g <- qpAnyGraph(abs(pcc.estimates$R), threshold=0.9,
                remove="below")

## the sparser the qp-graph
(sum(g)/2) / (nVar*(nVar-1)/2)

## the lower the threshold
g <- qpAnyGraph(abs(pcc.estimates$R), threshold=0.5,
                remove="below")

# the denser the graph
(sum(g)/2) / (nVar*(nVar-1)/2)
}
\keyword{models}
\keyword{multivariate}
