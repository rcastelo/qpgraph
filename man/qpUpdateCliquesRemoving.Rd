\name{qpUpdateCliquesRemoving}
\alias{qpUpdateGetCliquesRemoving}

\title{
Update clique list when removing one edge
}
\description{
Updates the set of (maximal) cliques of a given undirected graph when removing one edge.
}
\usage{
qpUpdateCliquesRemoving <- function(g, clqlst, v, w, verbose=TRUE)
}
\arguments{
  \item{g}{either a \code{graphNEL} object or an adjacency matrix of the given
       undirected graph.}
  \item{clqlst}{list of cliques of the graph encoded in g. this list should
                start on element n+1 (for n vertices) while between elements
                1 to n there should be references to the cliques to which each
                of the 1 to n vertices belong to (i.e., the output of
                \code{\link{qpGetCliques}}) with parameter \code{clqspervtx=TRUE}.}
  \item{v}{vertex of the edge being removed.}
  \item{w}{vertex of the edge being removed.}
  \item{verbose}{show progress on calculations.}
}
\details{
To find the list of all (maximal) cliques in an undirected graph is an NP-hard
problem which means that its computational cost is bounded by an exponential
running time (Garey and Johnson, 1979). For this reason, this is an extremely
time and memory consuming computation for large dense graphs. If we spend the time
to obtain one such list of cliques and we remove one edge of the graph with this
function we may be able to update the set of maximal cliques instead of having
to generate it again entirely with \code{\link{qpGetCliques}} but it requires that
in the first call to \code{\link{qpGetCliques}} we set \code{clqspervtx=TRUE}.
It calls a C implementation of the algorithm from Stix (2004).
}
\value{
A list of maximal cliques. When \code{clqspervtx=TRUE} the first p entries
(p=number of variables) contain, each of them, the indices of the cliques where
that particular vertex belongs to.
}
\references{
Garey, M.R. and Johnson D.S. \emph{Computers and intractability: a guide to the
theory of NP-completeness}. W.H. Freeman, San Francisco, 1979.

Stix, V. Finding all maximal cliques in dynamic graphs
\emph{Comput. Optimization and Appl.}, 27:173-186, 2004.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpCliqueNumber}}
  \code{\link{qpGetCliques}}
  \code{\link{qpIPF}}
}
\examples{
require(graph)

## this example should be updated
set.seed(123)
nVar <- 50
g1 <- randomEGraph(V=as.character(1:nVar), p=0.3)
clqs1 <- qpGetCliques(g1, verbose=FALSE)

length(clqs1)

summary(sapply(clqs1, length))
}
\keyword{models}
\keyword{multivariate}
