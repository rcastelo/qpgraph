\name{qpPlotMap}
\alias{qpPlotMap}

\title{
Plots a map of associated pairs
}
\description{
Plots a map of associated pairs defined by adjusted p-values
}
\usage{
qpPlotMap(p.valueMatrix, markerPos, genePos, chrLen,
          p.value=0.05, adjust.method="holm",
          xlab="Ordered Markers", ylab="Ordered Genes",
          main="", ...)
}
\arguments{
  \item{p.valueMatrix}{squared symmetric matrix with raw p-values for all pairs.}
  \item{markerPos}{two-column matrix containing chromosome and position of each genetic marker.}
  \item{genePos}{two-column matrix containing chromosome and position of each gene.}
  \item{chrLen}{named vector with chromosome lengths. Vector names should correspond to chromosome names, which are displayed in the axes of the plot. This vector should be ordered following the same convention for chromosomes in arguments \code{markerPos} and \code{genePos}.}
  \item{p.value}{adjusted p-value cutoff.}
  \item{adjust.method}{method employed to adjust the raw p-values. It is passed in a call to \code{p.adjust()} in its \code{method} argument.}
  \item{xlab}{label for the x-axis.}
  \item{ylab}{label for the y-axis.}
  \item{main}{main title of the plot, set to the empty string by default.}
  \item{...}{further arguments passed to the \code{plot()} function.}
}
\details{
This function plots a map of present associations, typically between genetic markers and gene expression profiles (i.e., eQTL associations), according to the chromosomal locations of both the genetic markers and the genes. The input argument \code{p.valueMatrix} should contain the raw p-values of these associations. Present associations are selected by a cutoff given in the \code{p.value} argument applied to the adjusted p-values.

The input raw p-values can be obtained with the function \code{\link{qpAllCItests}}.
}
\value{
The selected present associations are invisibly returned.
}
\author{R. Castelo}
\seealso{
  \code{\link{qpAllCItests}}
}
\examples{
## generate uniformly random p-values for synthetic associations
## between m genetic markers and g genes into a symmetric matrix
m <- 100
g <- 100
p <- m + g
markerids <- paste0("m", 1:m)
geneids <- paste0("g", 1:g)
rndpvalues <- matrix(0, nrow=p, ncol=p,
                     dimnames=list(c(markerids, geneids), c(markerids, geneids)))
rndpvalues[1:m,(m+1):p] <- runif(m*g)

## put significant cis associations
rndpvalues[cbind(1:m, (m+1):p)] <- rnorm(m, mean=1e-4, sd=1e-2)^2

## put one hotspot locus with significant, but somehat weaker, trans associations
hotspotmarker <- sample(1:m, size=1)
rndpvalues[cbind(hotspotmarker, (m+1):p)] <- rnorm(g, mean=1e-2, sd=1e-2)^2

## make matrix symmetric
rndpvalues <- rndpvalues + t(rndpvalues)
stopifnot(isSymmetric(rndpvalues))
rndpvalues[1:m, 1:m] <- rndpvalues[(m+1):p,(m+1):p] <- NA

## create chromosomal map
chrlen <- c("chr1"=1000)
posmarkers <- matrix(c(rep(1, m), seq(1, chrlen, length.out=m)), nrow=m)
posgenes <- matrix(c(rep(1, g), seq(1, chrlen, length.out=g)), nrow=g)
rownames(posmarkers) <- paste0("m", 1:m)
rownames(posgenes) <- paste0("g", 1:g)

qpPlotMap(rndpvalues, posmarkers, posgenes, chrlen, cex=3)
}
\keyword{models}
\keyword{multivariate}
