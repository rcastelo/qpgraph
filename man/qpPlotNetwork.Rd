\name{qpPlotNetwork}
\alias{qpPlotNetwork}

\title{
Plots a graph
}
\description{
Plots a graph using the \code{Rgraphviz} library
}
\usage{
qpPlotNetwork(g, vertexSubset=graph::nodes(g), boundary=FALSE,
              minimumSizeConnComp=2, pairup.i=NULL, pairup.j=NULL,
              highlight=NULL, annotation=NULL, layout=c("twopi", "dot", "neato", "circo", "fdp"))
}
\arguments{
  \item{g}{graph to plot provided as a \code{graphNEL-class} object.}
  \item{vertexSubset}{subset of vertices that define the induced subgraph to be plotted.}
  \item{boundary}{flag set to \code{TRUE} when we wish that the subset specified in \code{vertexSubset} also includes the vertices connected to them; \code{FALSE} otherwise.}
  \item{minimumSizeConnComp}{minimum size of the connected components to be plotted.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}.}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}.}
  \item{highlight}{subset of vertices to highlight by setting the color font to red.}
  \item{annotation}{name of an annotation package to transform gene identifiers into gene symbols when vertices correspond to genes.}
  \item{layout}{layout argument for the Rgraphviz library that plots the network. Possible values are \code{twopi} (default), \code{dot}, \code{neato}, \code{circo}, \code{fdp}.}
}
\details{
This function acts as a wrapper for the functionality provided by the \code{Rgraphviz}
package to plot graphs in R. It should help to plot networks obtained with methods from
the\code{qpgraph} package.
}
\value{
The plotted graph is invisibly returned as a \code{graphNEL-class} object.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpGraph}}
  \code{\link{qpAnyGraph}}
}
\examples{
\dontrun{
require(Rgraphviz)

rndassociations <- qpUnifRndAssociation(10)
g <- qpAnyGraph(abs(rndassociations), threshold=0.7, remove="below")
qpPlotNetwork(g) ## this does not work at the moment and should be fixed
}
}
\keyword{models}
\keyword{multivariate}
