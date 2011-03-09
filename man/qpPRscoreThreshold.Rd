\name{qpPRscoreThreshold}
\alias{qpPRscoreThreshold}

\title{
Calculation of scores thresholds attaining nominal precision or recall levels
}
\description{
Calculates the score threshold at a given precision or recall level from a given
precision-recall curve.
}
\usage{
qpPRscoreThreshold(preRecFun, level, recall.level=TRUE, max.score=9999999)
}
\arguments{
  \item{preRecFun}{precision-recall function (output from
       \code{\link{qpPrecisionRecall}}).}
  \item{level}{recall or precision level.}
  \item{recall.level}{logical; if TRUE then it is assumed that the value given in
       the level parameter corresponds to a desired level of recall; if FALSE
       then it is assumed a desired level of precision.}
  \item{max.score}{maximum score given by the method that produced the
       precision-recall function to an association.}
}
\value{
The score threshold at which a given level of precision or recall is attained by
the given precision-recall function. For levels that do not form part of the
given function their score is calculated by linear interpolation and for this
reason is important to carefully specify a proper value for the \code{max.score}
parameter.
}
\references{
Fawcett, T. An introduction to ROC analysis.
\emph{Pattern Recogn. Lett.}, 27:861-874, 2006.
}
\author{R. Castelo and A. Roverato}
\seealso{
  \code{\link{qpPrecisionRecall}}
  \code{\link{qpGraph}}
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

nrr.estimates <- qpNrr(X, q=1, verbose=FALSE)

nrr.prerec <- qpPrecisionRecall(nrr.estimates, A, decreasing=FALSE,
                                recallSteps=seq(0, 1, by=0.1))

qpPRscoreThreshold(nrr.prerec, level=0.5, recall.level=TRUE, max.score=0)

qpPRscoreThreshold(nrr.prerec, level=0.5, recall.level=FALSE, max.score=0)
}
\keyword{models}
\keyword{multivariate}
