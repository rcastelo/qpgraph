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
qpPrecisionRecall(measurementsMatrix, refI, decreasing=TRUE, pairup.i=NULL,
                  pairup.j=NULL, recallSteps=c(seq(0,0.1,0.005),seq(0.2,1.0,0.1)))
}
\arguments{
  \item{measurementsMatrix}{matrix containing the measure of association between
       all pairs of variables.}
  \item{refI}{incidence matrix of reference from which to calculate the
       precision-recall curve.}
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
nVar <- 50 # number of variables
maxCon <- 5  # maximum connectivity per variable
nObs <- 30 # number of observations to simulate

I <- qpRndGraph(n.vtx=nVar, n.bd=maxCon)
K <- qpI2K(I)

X <- qpSampleMvnorm(K, nObs)

nrr.estimates <- qpNrr(X, q=5, verbose=FALSE) # NRR !

pcc.estimates <- qpPCC(X) # PCC !

# calculate area under the precision-recall curve

nrr.prerec <- qpPrecisionRecall(nrr.estimates, refI=K!=0, decreasing=FALSE,
                                recallSteps=seq(0, 1, 0.1))
f <- approxfun(nrr.prerec[, c("Recall", "Precision")])
integrate(f, 0, 1)$value

pcc.prerec <- qpPrecisionRecall(abs(pcc.estimates$R), refI=K!=0,
                                recallSteps=seq(0, 1, 0.1))
f <- approxfun(pcc.prerec[, c("Recall", "Precision")])
integrate(f, 0, 1)$value
}
\keyword{models}
\keyword{multivariate}
