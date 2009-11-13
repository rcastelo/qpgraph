\name{qpRndGraph}
\alias{qpRndGraph}

\title{
Random undirected graphs with maximum connectivity degree
}
\description{
Builds a random undirected graph with a bounded maximum connectivity degree
(boundary) on every vertex.
}
\usage{
qpRndGraph(n.vtx, n.bd)
}
\arguments{
  \item{n.vtx}{number of vertices.}
  \item{n.bd}{maximum boundary for every vertex.}
}
\details{
This is a very simple function to generate random undirected graphs where
we impose a maximum order of correlation between disconnected vertices
when using it to sample multivariate normal data reflecting the conditional
independencies encoded in this graph. Note that the maximum order of
correlation between two disconnected vertices is bounded by the minimum
degree of connectivity of the two vertices.

The algorithm employed is not designed to ensure a uniform probability
distribution on the set of graphs with the given maximum boundary that may be
sampled with positive probability.
}
\value{
The adjacency matrix of the resulting graph.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpNrr}}
}
\examples{
nVar <- 50  ## number of vertices
maxCon <- 5 ## maximum connectivity per vertex

A <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)

summary(rowSums(A))
}
\keyword{models}
\keyword{multivariate}
