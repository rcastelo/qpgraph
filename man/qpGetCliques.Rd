\name{qpGetCliques}
\alias{qpGetCliques}

\title{
Clique list
}
\description{
Finds the set of (maximal) cliques of a given undirected graph.
}
\usage{
qpGetCliques(g, clqspervtx=FALSE, verbose=TRUE)
}
\arguments{
  \item{g}{either a \code{graphNEL} object or an incidence matrix of the given
       undirected graph.}
  \item{clqspervtx}{logical; if TRUE then the resulting list returned by the
       function includes additionally p entries at the beginning (p=number of
       variables) each corresponding to a vertex in the graph and containing the
       indices of the cliques where that vertex belongs to; if FALSE these
       additional entries are not included (default).}
  \item{verbose}{show progress on calculations.}
}
\details{
To find the list of all (maximal) cliques in an undirected graph is an NP-complete
problem which means that its computational cost is bounded by an exponential
running time (Karp, 1972). For this reason, this is an extremely time and memory
consuming computation for large dense graphs. The current implementation uses C
code from the GNU GPL Cliquer library by Niskanen and Ostergard (2003).

}
\value{
A list of maximal cliques. When \code{clqspervtx=TRUE} the first p entries
(p=number of variables) contain, each of them, the indices of the cliques where
that particular vertex belongs to.
}
\references{
Castelo, R. and Roverato, A. A robust procedure for
Gaussian graphical model search from microarray data with p larger than n.
\emph{J. Mach. Learn. Res.}, 7:2621-2650, 2006.

Niskanen, S. Ostergard, P. Cliquer User's Guide, Version 1.0.
Communications Laboratory, Helsinki University of Technology, Espoo, Finland,
Tech. Rep. T48, 2003. (\url{http://users.tkk.fi/~pat/cliquer.html})

Karp, R.M. Reducibility among combinatorial problems.
In R.E. Miller and J.W. Thatcher (eds.): \emph{Complexity of Computer
Computations}, 85-103, New York: Plenum, 1972.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpCliqueNumber}}
  \code{\link{qpIPF}}
}
\examples{
nVertices <- 50 # number of vertices
maxCon <- 5  # maximum connectivity per vertex

I <- qpRndGraph(n.vtx=nVertices, n.bd=maxCon)

clqs <- qpGetCliques(I, verbose=FALSE)

summary(unlist(lapply(clqs, length)))

maxCon <- 10  # maximum connectivity per vertex

I <- qpRndGraph(n.vtx=nVertices, n.bd=maxCon)

clqs <- qpGetCliques(I, verbose=FALSE)

summary(unlist(lapply(clqs, length)))

}
\keyword{models}
\keyword{multivariate}
