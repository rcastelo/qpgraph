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
  \item{g}{either a \code{graphNEL} object or an adjacency matrix of the given
       undirected graph.}
  \item{clqspervtx}{logical; if TRUE then the resulting list returned by the
       function includes additionally p entries at the beginning (p=number of
       variables) each corresponding to a vertex in the graph and containing the
       indices of the cliques where that vertex belongs to; if FALSE these
       additional entries are not included (default).}
  \item{verbose}{show progress on calculations.}
}
\details{
To find the list of all (maximal) cliques in an undirected graph is an NP-hard
problem which means that its computational cost is bounded by an exponential
running time (Garey and Johnson, 1979). For this reason, this is an extremely
time and memory consuming computation for large dense graphs. The current
implementation uses C code from the GNU GPL Cliquer library by Niskanen and
Ostergard (2003).
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

Garey, M.R. and Johnson D.S. \emph{Computers and intractability: a guide to the
theory of NP-completeness}. W.H. Freeman, San Francisco, 1979.

Niskanen, S. Ostergard, P. Cliquer User's Guide, Version 1.0.
Communications Laboratory, Helsinki University of Technology, Espoo, Finland,
Tech. Rep. T48, 2003. (\url{http://users.tkk.fi/~pat/cliquer.html})
}
\author{R. Castelo}
\seealso{
  \code{\link{qpCliqueNumber}}
  \code{\link{qpIPF}}
}
\examples{
require(graph)

set.seed(123)
nVar <- 50
g1 <- randomEGraph(V=as.character(1:nVar), p=0.3)
clqs1 <- qpGetCliques(g1, verbose=FALSE)

length(clqs1)

summary(sapply(clqs1, length))

g2 <- randomEGraph(V=as.character(1:nVar), p=0.7)
clqs2 <- qpGetCliques(g2, verbose=FALSE)

length(clqs2)

clqs2 <- qpGetCliques(g2, verbose=FALSE)

summary(sapply(clqs2, length))
}
\keyword{models}
\keyword{multivariate}
