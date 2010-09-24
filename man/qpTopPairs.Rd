\name{qpTopPairs}
\alias{qpTopPairs}

\title{
Report pairs of variables
}
\description{
Report a top number of pairs of variables according to either an association measure and/or
occurring in a given reference graph.
}
\usage{
qpTopPairs(measurementsMatrix=NULL, refGraph=NULL, n=6L, file=NULL,
           decreasing=FALSE, pairup.i=NULL, pairup.j=NULL,
           annotation=NULL, fcOutput=NULL, fcOutput.na.rm=FALSE,
           digits=2)
}
\arguments{
  \item{measurementsMatrix}{matrix containing the measure of association between
       all pairs of variables.}
  \item{refGraph}{a reference graph containing the pairs that should be reported
       and provided either as an adjacency matrix, a \code{graphNEL-class} object
       or a \code{graphAM-class} object.}
  \item{n}{number of pairs to report, 6 by default, use \code{Inf} for reporting all of them.}
  \item{file}{file name to dump the pairs information as tab-separated column text.}
  \item{decreasing}{logical; if TRUE then the measurements are employed to be
       ordered in decreasing order; if FALSE then in increasing order.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}.}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}.}
  \item{annotation}{name of an annotation package to transform gene identifiers into gene symbols when variables correspond to genes.}
  \item{fcOutput}{output of \code{\link{qpFunctionalCoherence}}.}
  \item{fcOutput.na.rm}{flag set to TRUE when pairs with \code{NA} values from \code{fcOutput} should not be reported; \code{FALSE} (default) otherwise.}
  \item{digits}{number of decimal digits reported in the association measure and functional coherence values.}
}
\details{
The \code{measurementsMatrix} should be symmetric and may have also contain
\code{NA} values which will not be taken into account. That is an alternative
way to restricting the variable pairs with the parameters \code{pairup.i} and
\code{pairup.j}. The same holds for \code{refGraph}. One of these two, should
be specified.
}
\value{
The ranking of pairs is invisibly returned.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpGraph}}
  \code{\link{qpPrecisionRecall}}
  \code{\link{qpFunctionalCoherence}}
}
\examples{
qpTopPairs(matrix(runif(100), nrow=10, dimnames=list(1:10,1:10)))
}
\keyword{models}
\keyword{multivariate}
