\name{qpPrecisionRecall}
\alias{qpPrecisionRecall}

\title{
Calculation of precision-recall curves
}
\description{
Calculates the precision-recall curve (see Fawcett, 2006) for a given measure of
association between all pairs of variables in a matrix.
}
\usage{
qpPrecisionRecall(measurementsMatrix, refGraph, decreasing=TRUE, pairup.i=NULL,
                  pairup.j=NULL, recallSteps=seq(0, 1, by=0.1))
}
\arguments{
  \item{measurementsMatrix}{matrix containing the measure of association between
       all pairs of variables.}
  \item{refGraph}{a reference graph from which to calculate the precision-recall
       curve provided either as an adjacency matrix, a two-column matrix of edges,
       a \code{graphNEL-class} object or a \code{graphAM-class} object.}
  \item{decreasing}{logical; if TRUE then the measurements are ordered in
       decreasing order; if FALSE then in increasing order.}
  \item{pairup.i}{subset of vertices to pair up with subset \code{pairup.j}.}
  \item{pairup.j}{subset of vertices to pair up with subset \code{pairup.i}.}
  \item{recallSteps}{steps of the recall on which to calculate precision.}
}
\details{
The \code{measurementsMatrix} should be symmetric and may have also contain
\code{NA} values which will not be taken into account. That is an alternative
way to restricting the variable pairs with the parameters \code{pairup.i} and
\code{pairup.j}.
}
\value{
A matrix where rows correspond to recall steps and columns correspond,
respetively, to the actual recall, the precision, the number of true positives at
that recall rate and the threshold score that yields that recall rate.
}
\references{
Fawcett, T. An introduction to ROC analysis.
\emph{Pattern Recogn. Lett.}, 27:861-874, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpPRscoreThreshold}}
  \code{\link{qpGraph}}
  \code{\link{qpAvgNrr}}
  \code{\link{qpPCC}}
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

## estimate Pearson correlation coefficients
pcc.estimates <- qpPCC(X)

## calculate area under the precision-recall curve
## for both sets of estimated values of association
nrr.prerec <- qpPrecisionRecall(nrr.estimates, refGraph=A, decreasing=FALSE,
                                recallSteps=seq(0, 1, 0.1))
f <- approxfun(nrr.prerec[, c("Recall", "Precision")])
integrate(f, 0, 1)$value

pcc.prerec <- qpPrecisionRecall(abs(pcc.estimates$R), refGraph=A,
                                recallSteps=seq(0, 1, 0.1))
f <- approxfun(pcc.prerec[, c("Recall", "Precision")])
integrate(f, 0, 1)$value
}
\keyword{models}
\keyword{multivariate}
