\name{qpgraph-defunct}
\alias{qpgraph-defunct}
\title{Defunct Functions in Package \pkg{qpgraph}}
\description{
  These functions are provided for compatibility with older versions of
  \pkg{qpgraph} only, and may disappear as soon as the next release.
}
\usage{
qpRndHMGM(nDiscrete=1, nContinuous=3, d=2, mixedIntStrength=5, rho=0.5, G=NULL)
qpSampleFromHMGM(n=10, hmgm=qpRndHMGM())
}
\details{
  \code{qpRndHMGM()} has been replaced by \code{rHMgmm()}.

  \code{qpSampleFromHMGM()} has been replaced by \code{rcmvnorm()}.

  The original help page for these functions is available as
  \code{help("oldName-defunct")} (note the quotes).
}
\seealso{
  \code{\link{Defunct}}
}
\keyword{internal}
\keyword{misc}
