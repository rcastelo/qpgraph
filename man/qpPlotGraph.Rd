\name{qpPlotGraph}
\alias{qpPlotGraph}

\title{
Plots a graph
}
\description{
Plots a graph using the \code{Rgraphviz} library
}
\usage{
qpPlotGraph(g, vertexSubset=nodes(g), boundary=FALSE, minimumSizeConnComp=2,
            pairup.i=NULL, pairup.j=NULL, annotation=NULL)
}
\arguments{
  \item{g}{graph to plot provided as a \code{graphNEL-class} object.}
  \item{vertexSubset}{subset of vertices that define the induced subgraph to be plotted.}
  \item{boundary}{flag set to \code{TRUE} when we wish that the subset specified in \code{vertexSubset} also includes the vertices connected to them; \code{FALSE} otherwise.}
  \item{minimumSizeConnComp}{minimum size of the connected components to be plotted.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}.}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}.}
  \item{annotation}{name of an annotation package to transform gene identifiers into gene symbols when vertices correspond to genes.}
}
\details{
}
\value{
The plotted graph is invisibly returned as a \code{graphNEL-class} object.
}
\references{
}
\author{R. Castelo}
\seealso{
  \code{\link{qpGraph}}
}
\examples{
## missing
}
\keyword{models}
\keyword{multivariate}
